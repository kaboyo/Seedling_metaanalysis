import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import zipfile
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


RANK_COLUMNS = ["Phyla", "Family", "Genus"]
UNASSIGNED_BY_RANK = {
    "Phyla": "unassigned_phyla",
    "Family": "unassigned_family",
    "Genus": "unassigned_genus",
}


def load_config(config_path: Path) -> tuple[dict, Path]:
    config_path = Path(config_path).resolve()
    with config_path.open("r", encoding="utf-8") as handle:
        return json.load(handle), config_path.parent


def resolve_path(base_dir: Path, configured_path: str) -> Path:
    path = Path(configured_path)
    return path if path.is_absolute() else (base_dir / path).resolve()


#: Decoding order for the metadata table. Datesets/metadata.csv is MAC ROMAN, not
#: latin-1: the distinguishing byte is 0x96, which Mac Roman decodes to "n-tilde"
#: (Acuna) while latin-1 maps it to an unassigned control character and cp1252 to an
#: en-dash. Order matters and cannot be inferred, because latin-1 and mac_roman both
#: decode *every* byte without raising - whichever is tried first always wins. latin-1
#: first silently mangled 33 author names, so mac_roman now precedes it and latin-1
#: remains last as a never-fails backstop.
METADATA_ENCODINGS = ("utf-8", "mac_roman", "latin-1")


def load_metadata(path: Path, encodings: tuple[str, ...] = METADATA_ENCODINGS) -> pd.DataFrame:
    last_error = None
    for encoding in encodings:
        try:
            frame = pd.read_csv(path, encoding=encoding)
            break
        except UnicodeDecodeError as exc:
            last_error = exc
    else:
        raise last_error

    unnamed = [column for column in frame.columns if str(column).startswith("Unnamed:") or column == ""]
    if unnamed:
        frame = frame.drop(columns=unnamed)

    if "run_accession" not in frame.columns:
        raise ValueError("metadata file must contain 'run_accession'")

    frame["run_accession"] = frame["run_accession"].astype(str).str.strip()
    frame = frame[frame["run_accession"] != ""].copy()
    return frame


def extract_timepoint(metadata: pd.DataFrame) -> pd.DataFrame:
    """Extract Timepoint from sample_description for PRJNA649549 (Denmark) samples."""
    timepoints: list[str] = []
    for _, row in metadata.iterrows():
        study = str(row.get("study_alias", ""))
        desc = str(row.get("sample_description", ""))
        if study != "PRJNA649549":
            timepoints.append("")
            continue
        m = re.match(r"^seed:(T0)\.Seed:", desc)
        if m:
            timepoints.append(m.group(1))
            continue
        m = re.match(r"^seed:(\d+)\s*DAS:", desc)
        if m:
            timepoints.append(f"{m.group(1)}_DAS")
            continue
        m = re.match(r"^soil:(\d+)DBS:", desc)
        if m:
            timepoints.append(f"{m.group(1)}_DBS")
            continue
        timepoints.append("")
    metadata = metadata.copy()
    metadata["Timepoint"] = timepoints
    return metadata


def load_taxonomy(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    if "OTUID" not in frame.columns:
        raise ValueError("taxonomy file must contain 'OTUID'")
    frame["OTUID"] = frame["OTUID"].astype(str)
    return frame.set_index("OTUID")


def extract_newick_text(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in {".nexus", ".nex"}:
        text = path.read_text(encoding="utf-8")
        match = re.search(r"tree\s+[^=]+=\s*(?:\[[^\]]*\]\s*)?(.+?);", text, flags=re.IGNORECASE | re.DOTALL)
        if not match:
            raise ValueError(f"could not find tree statement in {path}")
        return match.group(1).strip() + ";"

    if suffix in {".nwk", ".newick", ".tree"}:
        return path.read_text(encoding="utf-8").strip()

    if suffix == ".qza":
        with zipfile.ZipFile(path) as archive:
            candidates = [name for name in archive.namelist() if name.endswith("tree.nwk")]
            if not candidates:
                raise ValueError(f"could not find tree.nwk in {path}")
            return archive.read(candidates[0]).decode("utf-8").strip()

    raise ValueError(f"unsupported tree format: {path.suffix}")


def ensure_tree_nexus(tree_path: Path, output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if tree_path.suffix.lower() in {".nexus", ".nex"}:
        if tree_path.resolve() != output_path.resolve():
            output_path.write_text(tree_path.read_text(encoding="utf-8"), encoding="utf-8")
        return output_path

    newick = extract_newick_text(tree_path)
    nexus = "#NEXUS\nBegin trees;\nTree tree1 = [&R] " + newick + "\nEnd;\n"
    output_path.write_text(nexus, encoding="utf-8")
    return output_path


def extract_tree_tips(path: Path) -> set[str]:
    newick = extract_newick_text(path)
    quoted = re.findall(r"(?:(?<=\()|(?<=,))\s*'([^']+)'(?=\s*:)", newick)
    scrubbed = re.sub(r"'[^']+'", "", newick)
    unquoted = re.findall(r"(?:(?<=\()|(?<=,))\s*([^'(),:;\s]+)(?=\s*:)", scrubbed)
    return set(quoted) | set(unquoted)


def build_taxonomy_mask(
    feature_ids: list[str],
    taxonomy: pd.DataFrame,
    filters: dict,
) -> np.ndarray:
    """Return a boolean mask (len == len(feature_ids)) marking features that pass
    taxonomy filters: keep only allowed kingdoms, exclude chloroplasts (Order)
    and mitochondria (Family).

    Config keys used from *filters*:
      keep_kingdom   – list[str], e.g. ["Bacteria"]  (empty/absent = keep all)
      exclude_order  – list[str], e.g. ["Chloroplast"]
      exclude_family – list[str], e.g. ["Mitochondria"]
    """
    n = len(feature_ids)
    mask = np.ones(n, dtype=bool)

    keep_kingdoms: list[str] = filters.get("keep_kingdom", [])
    exclude_orders: list[str] = filters.get("exclude_order", [])
    exclude_families: list[str] = filters.get("exclude_family", [])

    if not (keep_kingdoms or exclude_orders or exclude_families):
        return mask  # nothing to filter

    # Build sets for fast lookup (case-sensitive – matches taxonomy file)
    keep_kingdoms_set = set(keep_kingdoms) if keep_kingdoms else None
    exclude_orders_set = set(exclude_orders)
    exclude_families_set = set(exclude_families)

    for idx, fid in enumerate(feature_ids):
        if fid not in taxonomy.index:
            # Feature not in taxonomy table → exclude if kingdom filter is active
            if keep_kingdoms_set:
                mask[idx] = False
            continue

        row = taxonomy.loc[fid]

        # Kingdom filter
        if keep_kingdoms_set and row.get("Kingdom", "") not in keep_kingdoms_set:
            mask[idx] = False
            continue

        # Order filter (chloroplasts)
        if exclude_orders_set and row.get("Order", "") in exclude_orders_set:
            mask[idx] = False
            continue

        # Family filter (mitochondria)
        if exclude_families_set and row.get("Family", "") in exclude_families_set:
            mask[idx] = False
            continue

    return mask


def parse_count(value: str) -> int:
    value = value.strip()
    if not value:
        return 0
    return int(float(value))


def shannon_index(counts: np.ndarray) -> float:
    total = counts.sum()
    if total <= 0:
        return 0.0
    proportions = counts[counts > 0] / total
    return float(-(proportions * np.log(proportions)).sum())


def simpson_index(counts: np.ndarray) -> float:
    total = counts.sum()
    if total <= 0:
        return 0.0
    proportions = counts[counts > 0] / total
    return float(1.0 - np.square(proportions).sum())


def safe_pielou(shannon: float, richness: int) -> float:
    if richness <= 1:
        return 0.0
    return float(shannon / np.log(richness))


def first_pass_counts(
    counts_path: Path,
    max_samples: int | None = None,
    taxonomy_mask: np.ndarray | None = None,
) -> dict:
    """Stream the ASV count table and compute per-sample / per-feature summaries.

    Parameters
    ----------
    taxonomy_mask : optional boolean array (len == n_features).
        When provided, only features where mask is True contribute to
        per-sample alpha metrics (library_size, richness, shannon, simpson).
        Per-feature prevalence/totals are still recorded for ALL features so
        that `select_feature_indices` can combine the taxonomy mask with
        prevalence/abundance filters later.
    """
    with counts_path.open("r", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        feature_ids = [item.strip() for item in header[1:]]
        n_features = len(feature_ids)
        feature_prevalence = np.zeros(n_features, dtype=np.int64)
        feature_totals = np.zeros(n_features, dtype=np.int64)
        sample_ids: list[str] = []
        library_sizes: list[int] = []
        richness: list[int] = []
        shannon: list[float] = []
        simpson: list[float] = []
        blank_cells = 0
        zero_rows = 0

        # Pre-compute mask indices for alpha-metric calculation
        if taxonomy_mask is not None:
            taxon_idx = np.flatnonzero(taxonomy_mask)
        else:
            taxon_idx = None

        for row_index, row in enumerate(reader):
            if max_samples is not None and row_index >= max_samples:
                break
            sample_ids.append(row[0].strip())
            counts = np.fromiter(
                (parse_count(v) for v in row[1:]), dtype=np.int64, count=n_features
            )
            blank_cells += sum(1 for value in row[1:] if value.strip() == "")

            # Per-feature stats: ALL features (unfiltered)
            nonzero_all = counts > 0
            feature_prevalence += nonzero_all.astype(np.int64)
            feature_totals += counts

            # Per-sample alpha metrics: taxonomy-filtered features only
            if taxon_idx is not None:
                filtered_counts = counts[taxon_idx]
            else:
                filtered_counts = counts

            nonzero = filtered_counts > 0
            sample_total = int(filtered_counts.sum())
            sample_richness = int(nonzero.sum())
            library_sizes.append(sample_total)
            richness.append(sample_richness)
            shannon_value = shannon_index(filtered_counts)
            shannon.append(shannon_value)
            simpson.append(simpson_index(filtered_counts))
            if sample_total == 0:
                zero_rows += 1

    return {
        "feature_ids": feature_ids,
        "feature_prevalence": feature_prevalence,
        "feature_totals": feature_totals,
        "sample_ids": sample_ids,
        "library_sizes": np.array(library_sizes, dtype=np.int64),
        "richness": np.array(richness, dtype=np.int64),
        "shannon": np.array(shannon, dtype=float),
        "simpson": np.array(simpson, dtype=float),
        "blank_cells": blank_cells,
        "zero_rows": zero_rows,
    }


def select_feature_indices(
    first_pass: dict,
    filters: dict,
    taxonomy_mask: np.ndarray | None = None,
) -> np.ndarray:
    """Select feature indices that pass prevalence, abundance, AND taxonomy filters."""
    prevalence = first_pass["feature_prevalence"]
    totals = first_pass["feature_totals"]
    keep_mask = (
        (prevalence >= int(filters["min_feature_prevalence"]))
        & (totals >= int(filters["min_feature_total_count"]))
    )

    # Apply taxonomy mask (Bacteria only, no chloroplast/mitochondria)
    if taxonomy_mask is not None:
        keep_mask &= taxonomy_mask

    retained = np.flatnonzero(keep_mask)
    max_features = filters.get("max_retained_features")
    if max_features and len(retained) > int(max_features):
        score = np.lexsort((totals[retained], prevalence[retained]))
        retained = retained[score[-int(max_features):]]
        retained = np.sort(retained)
    return retained


def second_pass_collect(
    counts_path: Path,
    sample_ids: list[str],
    feature_ids: list[str],
    selected_indices: np.ndarray,
    taxonomy: pd.DataFrame,
    taxonomy_mask: np.ndarray | None = None,
) -> dict:
    """Build dense matrices and rank-collapsed abundance tables.

    Parameters
    ----------
    taxonomy_mask : optional boolean array (len == n_features).
        When provided, rank-collapsed tables are built using only features
        that pass the taxonomy filter (Bacteria, no chloroplast/mitochondria).
    """
    feature_to_rank: dict[str, list[str]] = {}
    for rank in RANK_COLUMNS:
        default = UNASSIGNED_BY_RANK[rank]
        feature_to_rank[rank] = [taxonomy[rank].get(feature_id, default) for feature_id in feature_ids]

    selected_feature_ids = [feature_ids[idx] for idx in selected_indices]
    sample_count = len(sample_ids)
    selected_matrix = np.zeros((sample_count, len(selected_indices)), dtype=np.float32)
    rank_tables: dict[str, list[dict[str, int]]] = {rank: [] for rank in RANK_COLUMNS}
    unassigned_fraction = np.zeros(sample_count, dtype=float)

    selected_lookup = {feature_index: target_index for target_index, feature_index in enumerate(selected_indices.tolist())}

    with counts_path.open("r", newline="") as handle:
        reader = csv.reader(handle)
        next(reader)
        for row_index, row in enumerate(reader):
            if row_index >= sample_count:
                break
            counts = np.fromiter((parse_count(v) for v in row[1:]), dtype=np.int64, count=len(feature_ids))

            if len(selected_indices):
                selected_matrix[row_index, :] = counts[selected_indices]

            nonzero_indices = np.flatnonzero(counts)
            rank_totals: dict[str, dict[str, int]] = {rank: {} for rank in RANK_COLUMNS}
            unassigned = 0
            filtered_total = 0
            for feature_index in nonzero_indices:
                # Skip features that didn't pass taxonomy filter
                if taxonomy_mask is not None and not taxonomy_mask[feature_index]:
                    continue
                count = int(counts[feature_index])
                filtered_total += count
                for rank in RANK_COLUMNS:
                    taxon = feature_to_rank[rank][feature_index]
                    rank_totals[rank][taxon] = rank_totals[rank].get(taxon, 0) + count
                if feature_to_rank["Genus"][feature_index] == UNASSIGNED_BY_RANK["Genus"]:
                    unassigned += count

            for rank in RANK_COLUMNS:
                rank_tables[rank].append(rank_totals[rank])
            unassigned_fraction[row_index] = float(unassigned / filtered_total) if filtered_total > 0 else 0.0

    rank_frames = {}
    relative_rank_frames = {}
    for rank, rows in rank_tables.items():
        frame = pd.DataFrame(rows, index=sample_ids).fillna(0.0).astype(float)
        frame.index.name = "run_accession"
        rank_frames[rank] = frame
        relative = frame.div(frame.sum(axis=1).replace(0, 1), axis=0)
        relative_rank_frames[rank] = relative

    return {
        "selected_matrix": selected_matrix,
        "selected_feature_ids": selected_feature_ids,
        "rank_frames": rank_frames,
        "relative_rank_frames": relative_rank_frames,
        "unassigned_fraction": unassigned_fraction,
    }


def bray_curtis_distance(matrix: np.ndarray, block_size: int = 32) -> np.ndarray:
    matrix = matrix.astype(np.float32, copy=False)
    n_samples = matrix.shape[0]
    row_sums = matrix.sum(axis=1)
    distances = np.zeros((n_samples, n_samples), dtype=np.float32)
    for i in range(0, n_samples, block_size):
        x = matrix[i : i + block_size]
        sx = row_sums[i : i + block_size][:, None]
        for j in range(i, n_samples, block_size):
            y = matrix[j : j + block_size]
            sy = row_sums[j : j + block_size][None, :]
            mins = np.minimum(x[:, None, :], y[None, :, :]).sum(axis=2)
            denom = sx + sy
            ratio = np.divide(2.0 * mins, denom, out=np.zeros_like(mins, dtype=np.float32), where=denom > 0)
            block = np.where(denom > 0, 1.0 - ratio, 0.0)
            distances[i : i + block.shape[0], j : j + block.shape[1]] = block
            if i != j:
                distances[j : j + block.shape[1], i : i + block.shape[0]] = block.T
    return distances


def jaccard_distance(matrix: np.ndarray, block_size: int = 128) -> np.ndarray:
    binary = (matrix > 0).astype(np.int16)
    n_samples = binary.shape[0]
    row_sums = binary.sum(axis=1)
    distances = np.zeros((n_samples, n_samples), dtype=np.float32)
    for i in range(0, n_samples, block_size):
        x = binary[i : i + block_size]
        sx = row_sums[i : i + block_size][:, None]
        for j in range(i, n_samples, block_size):
            y = binary[j : j + block_size]
            sy = row_sums[j : j + block_size][None, :]
            intersection = x @ y.T
            union = sx + sy - intersection
            ratio = np.divide(intersection, union, out=np.zeros_like(intersection, dtype=np.float32), where=union > 0)
            block = np.where(union > 0, 1.0 - ratio, 0.0)
            distances[i : i + block.shape[0], j : j + block.shape[1]] = block
            if i != j:
                distances[j : j + block.shape[1], i : i + block.shape[0]] = block.T
    return distances


def pcoa(distance_matrix: np.ndarray, axes: int = 2) -> tuple[pd.DataFrame, np.ndarray]:
    n_samples = distance_matrix.shape[0]
    squared = np.square(distance_matrix.astype(np.float64))
    centering = np.eye(n_samples) - np.ones((n_samples, n_samples)) / n_samples
    b_matrix = -0.5 * centering @ squared @ centering
    eigenvalues, eigenvectors = np.linalg.eigh(b_matrix)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    positive = eigenvalues > 0
    positive_values = eigenvalues[positive]
    positive_vectors = eigenvectors[:, positive]
    coords = positive_vectors[:, :axes] * np.sqrt(positive_values[:axes])
    columns = [f"PC{i + 1}" for i in range(coords.shape[1])]
    explained = positive_values[:axes] / positive_values.sum() if positive_values.size else np.array([])
    return pd.DataFrame(coords, columns=columns), explained


def permanova(distance_matrix: np.ndarray, groups: list[str], permutations: int = 99, seed: int = 0) -> dict:
    groups = np.array(groups)
    unique_groups, counts = np.unique(groups, return_counts=True)
    if len(unique_groups) < 2 or np.any(counts < 2):
        return {"f_stat": np.nan, "p_value": np.nan, "groups": len(unique_groups), "status": "insufficient_groups"}

    n = len(groups)
    g = len(unique_groups)
    centered = np.eye(n) - np.ones((n, n)) / n
    a = -0.5 * np.square(distance_matrix.astype(np.float64))
    gram = centered @ a @ centered
    total_ss = np.trace(gram)

    def _within_sum(labels: np.ndarray) -> float:
        total = 0.0
        for value in np.unique(labels):
            idx = np.flatnonzero(labels == value)
            sub = gram[np.ix_(idx, idx)]
            total += np.trace(sub) - sub.sum() / len(idx)
        return float(total)

    within_ss = _within_sum(groups)
    among_ss = total_ss - within_ss
    df_between = g - 1
    df_within = n - g
    if df_within <= 0 or within_ss <= 0:
        return {"f_stat": np.nan, "p_value": np.nan, "groups": len(unique_groups), "status": "degenerate"}
    observed_f = (among_ss / df_between) / (within_ss / df_within)

    rng = np.random.default_rng(seed)
    permuted = 0
    for _ in range(permutations):
        shuffled = rng.permutation(groups)
        perm_within = _within_sum(shuffled)
        perm_among = total_ss - perm_within
        perm_f = (perm_among / df_between) / (perm_within / df_within) if perm_within > 0 else np.inf
        if perm_f >= observed_f:
            permuted += 1

    p_value = (permuted + 1) / (permutations + 1)
    return {"f_stat": float(observed_f), "p_value": float(p_value), "groups": len(unique_groups), "status": "ok"}


def permdisp(distance_matrix: np.ndarray, groups: list[str], permutations: int = 99, seed: int = 0) -> dict:
    groups = np.array(groups)
    unique_groups, counts = np.unique(groups, return_counts=True)
    if len(unique_groups) < 2 or np.any(counts < 2):
        return {"f_stat": np.nan, "p_value": np.nan, "groups": len(unique_groups), "status": "insufficient_groups"}

    coords, _ = pcoa(distance_matrix, axes=min(10, distance_matrix.shape[0] - 1))
    coord_values = coords.to_numpy()
    distances_to_centroid = np.zeros(len(groups), dtype=float)
    for value in unique_groups:
        idx = np.flatnonzero(groups == value)
        centroid = coord_values[idx].mean(axis=0)
        distances_to_centroid[idx] = np.linalg.norm(coord_values[idx] - centroid, axis=1)

    def _f_stat(labels: np.ndarray) -> float:
        grand_mean = distances_to_centroid.mean()
        ss_between = 0.0
        ss_within = 0.0
        for value in np.unique(labels):
            idx = np.flatnonzero(labels == value)
            group_mean = distances_to_centroid[idx].mean()
            ss_between += len(idx) * ((group_mean - grand_mean) ** 2)
            ss_within += ((distances_to_centroid[idx] - group_mean) ** 2).sum()
        df_between = len(np.unique(labels)) - 1
        df_within = len(labels) - len(np.unique(labels))
        if df_between <= 0 or df_within <= 0 or ss_within <= 0:
            return np.nan
        return float((ss_between / df_between) / (ss_within / df_within))

    observed_f = _f_stat(groups)
    rng = np.random.default_rng(seed)
    permuted = 0
    for _ in range(permutations):
        shuffled = rng.permutation(groups)
        perm_f = _f_stat(shuffled)
        if np.isfinite(perm_f) and perm_f >= observed_f:
            permuted += 1
    p_value = (permuted + 1) / (permutations + 1)
    return {"f_stat": float(observed_f), "p_value": float(p_value), "groups": len(unique_groups), "status": "ok"}


def kruskal_wallis_tests(
    alpha: pd.DataFrame,
    group_column: str,
    metrics: list[str] | None = None,
) -> pd.DataFrame:
    """Run Kruskal-Wallis H-test for each alpha diversity metric across groups."""
    if metrics is None:
        metrics = ["shannon", "observed_features", "simpson", "pielou_evenness"]
    metrics = [m for m in metrics if m in alpha.columns]
    if group_column not in alpha.columns:
        return pd.DataFrame()

    subset = alpha.dropna(subset=[group_column])
    groups = subset[group_column].unique()
    if len(groups) < 2:
        return pd.DataFrame()

    rows: list[dict] = []
    for metric in metrics:
        group_data = [
            subset.loc[subset[group_column] == g, metric].dropna().values for g in groups
        ]
        group_data = [g for g in group_data if len(g) >= 1]
        if len(group_data) < 2:
            rows.append({"metric": metric, "H_statistic": np.nan, "p_value": np.nan, "n_groups": len(groups)})
            continue
        h_stat, p_val = stats.kruskal(*group_data)
        rows.append({"metric": metric, "H_statistic": float(h_stat), "p_value": float(p_val), "n_groups": len(groups)})
    return pd.DataFrame(rows)


def dunns_posthoc(
    alpha: pd.DataFrame,
    group_column: str,
    metric: str = "shannon",
) -> pd.DataFrame:
    """Dunn's post-hoc test with Bonferroni correction for pairwise group comparisons."""
    if group_column not in alpha.columns or metric not in alpha.columns:
        return pd.DataFrame()

    subset = alpha[[group_column, metric]].dropna()
    groups = sorted(subset[group_column].unique())
    if len(groups) < 2:
        return pd.DataFrame()

    n_total = len(subset)
    ranked = stats.rankdata(subset[metric].values)
    subset = subset.copy()
    subset["_rank"] = ranked

    group_stats: dict[str, tuple[float, int]] = {}
    for g in groups:
        mask = subset[group_column] == g
        group_stats[g] = (float(subset.loc[mask, "_rank"].mean()), int(mask.sum()))

    tie_counts = np.unique(ranked, return_counts=True)[1]
    tie_correction = 1.0 - np.sum(tie_counts**3 - tie_counts) / (n_total**3 - n_total)

    rows: list[dict] = []
    pairs = list(combinations(groups, 2))
    n_comparisons = len(pairs)

    for g1, g2 in pairs:
        mean_r1, n1 = group_stats[g1]
        mean_r2, n2 = group_stats[g2]
        se = math.sqrt(
            tie_correction * (n_total * (n_total + 1) / 12.0) * (1.0 / n1 + 1.0 / n2)
        )
        if se <= 0:
            rows.append({"group1": g1, "group2": g2, "z_statistic": np.nan, "p_value": np.nan, "p_adjusted": np.nan})
            continue
        z = (mean_r1 - mean_r2) / se
        p_val = 2.0 * stats.norm.sf(abs(z))
        p_adj = min(p_val * n_comparisons, 1.0)
        rows.append({
            "group1": str(g1), "group2": str(g2),
            "z_statistic": float(z), "p_value": float(p_val), "p_adjusted": float(p_adj),
        })
    return pd.DataFrame(rows)


def nmds(distance_matrix: np.ndarray, dimensions: int = 2, max_iter: int = 300, seed: int = 0) -> pd.DataFrame:
    """Non-metric multidimensional scaling via SMACOF algorithm."""
    n = distance_matrix.shape[0]
    rng = np.random.default_rng(seed)
    coords = rng.standard_normal((n, dimensions)).astype(np.float64)

    target = distance_matrix.astype(np.float64)
    np.fill_diagonal(target, 0.0)

    for _ in range(max_iter):
        current_dist = np.sqrt(
            np.sum((coords[:, None, :] - coords[None, :, :]) ** 2, axis=2)
        )
        np.fill_diagonal(current_dist, 1.0)

        # Isotonic regression (pool adjacent violators on flattened upper triangle)
        triu_idx = np.triu_indices(n, k=1)
        target_flat = target[triu_idx]
        dist_flat = current_dist[triu_idx]
        order = np.argsort(target_flat)
        dist_ordered = dist_flat[order]

        # Simple isotonic regression
        isotonic = dist_ordered.copy()
        i = 0
        while i < len(isotonic) - 1:
            if isotonic[i] > isotonic[i + 1]:
                j = i + 1
                while j < len(isotonic) and isotonic[j] < isotonic[i]:
                    j += 1
                isotonic[i:j] = isotonic[i:j].mean()
            i += 1

        disparities = np.zeros_like(current_dist)
        reverse_order = np.argsort(order)
        disparities[triu_idx] = isotonic[reverse_order]
        disparities = disparities + disparities.T

        ratio = np.divide(
            disparities, current_dist,
            out=np.zeros_like(current_dist),
            where=current_dist > 0,
        )
        new_coords = np.zeros_like(coords)
        for idx in range(n):
            weights = ratio[idx]
            new_coords[idx] = coords[idx] + (weights[:, None] * (coords - coords[idx])).sum(axis=0) / max(weights.sum(), 1e-10)
            new_coords[idx] = new_coords[idx] - new_coords[idx].mean()

        shift = np.sqrt(np.sum((new_coords - coords) ** 2))
        coords = new_coords
        if shift < 1e-6:
            break

    # Compute stress
    final_dist = np.sqrt(np.sum((coords[:, None, :] - coords[None, :, :]) ** 2, axis=2))
    stress = np.sqrt(
        np.sum((target[triu_idx] - final_dist[triu_idx]) ** 2)
        / max(np.sum(final_dist[triu_idx] ** 2), 1e-10)
    )

    columns = [f"NMDS{i + 1}" for i in range(dimensions)]
    result = pd.DataFrame(coords, columns=columns)
    result.attrs["stress"] = float(stress)
    return result


def differential_abundance(
    abundance: pd.DataFrame,
    metadata: pd.DataFrame,
    group_column: str,
    group1: str,
    group2: str,
    min_prevalence: float = 0.1,
) -> pd.DataFrame:
    """CLR-based differential abundance testing with Welch's t-test and BH correction."""
    merged = abundance.join(
        metadata.set_index("run_accession")[[group_column]], how="inner"
    )
    subset = merged[merged[group_column].isin([group1, group2])].copy()
    if subset.empty:
        return pd.DataFrame()

    taxa_cols = [c for c in abundance.columns if c != group_column]

    # Filter by prevalence in at least one group
    keep_taxa: list[str] = []
    for taxon in taxa_cols:
        prev1 = (subset.loc[subset[group_column] == group1, taxon] > 0).mean()
        prev2 = (subset.loc[subset[group_column] == group2, taxon] > 0).mean()
        if max(prev1, prev2) >= min_prevalence:
            keep_taxa.append(taxon)

    if not keep_taxa:
        return pd.DataFrame()

    # CLR transform (add pseudocount)
    counts = subset[keep_taxa].values.astype(np.float64) + 1e-6
    log_counts = np.log(counts)
    geometric_mean = log_counts.mean(axis=1, keepdims=True)
    clr = log_counts - geometric_mean

    clr_df = pd.DataFrame(clr, columns=keep_taxa, index=subset.index)
    clr_df[group_column] = subset[group_column].values

    rows: list[dict] = []
    g1_mask = clr_df[group_column] == group1
    g2_mask = clr_df[group_column] == group2

    for taxon in keep_taxa:
        vals1 = clr_df.loc[g1_mask, taxon].values
        vals2 = clr_df.loc[g2_mask, taxon].values
        mean1 = vals1.mean()
        mean2 = vals2.mean()
        t_stat, p_val = stats.ttest_ind(vals1, vals2, equal_var=False)
        log2fc = float((mean1 - mean2) / np.log(2))
        rows.append({
            "taxon": taxon,
            "mean_clr_group1": float(mean1),
            "mean_clr_group2": float(mean2),
            "log2_fold_change": log2fc,
            "t_statistic": float(t_stat),
            "p_value": float(p_val),
        })

    result = pd.DataFrame(rows)
    if result.empty:
        return result

    # Benjamini-Hochberg correction
    p_values = result["p_value"].values
    n_tests = len(p_values)
    ranked = stats.rankdata(p_values)
    adjusted = np.minimum(p_values * n_tests / ranked, 1.0)
    # Enforce monotonicity
    order = np.argsort(p_values)[::-1]
    adjusted_sorted = adjusted[order]
    for i in range(1, len(adjusted_sorted)):
        adjusted_sorted[i] = min(adjusted_sorted[i], adjusted_sorted[i - 1])
    adjusted[order] = adjusted_sorted
    result["p_adjusted"] = adjusted
    result = result.sort_values("p_value").reset_index(drop=True)
    return result


def core_microbiome(
    abundance: pd.DataFrame,
    metadata: pd.DataFrame,
    group_column: str,
    prevalence_threshold: float = 0.5,
    abundance_threshold: float = 0.001,
) -> pd.DataFrame:
    """Identify core taxa present above thresholds in each group."""
    merged = abundance.join(
        metadata.set_index("run_accession")[[group_column]], how="inner"
    )
    groups = sorted(merged[group_column].dropna().unique())
    taxa_cols = [c for c in abundance.columns if c != group_column]

    rows: list[dict] = []
    for group in groups:
        group_data = merged[merged[group_column] == group][taxa_cols]
        n_samples = len(group_data)
        if n_samples == 0:
            continue
        for taxon in taxa_cols:
            prevalence = (group_data[taxon] > 0).sum() / n_samples
            mean_abundance = group_data[taxon].mean()
            if prevalence >= prevalence_threshold and mean_abundance >= abundance_threshold:
                rows.append({
                    "group": str(group),
                    "taxon": taxon,
                    "prevalence": float(prevalence),
                    "mean_relative_abundance": float(mean_abundance),
                })
    return pd.DataFrame(rows)


def indicator_species(
    abundance: pd.DataFrame,
    metadata: pd.DataFrame,
    group_column: str,
    permutations: int = 999,
    seed: int = 0,
) -> pd.DataFrame:
    """Indicator Value (IndVal) analysis to identify group-specific taxa."""
    merged = abundance.join(
        metadata.set_index("run_accession")[[group_column]], how="inner"
    )
    groups = sorted(merged[group_column].dropna().unique())
    taxa_cols = [c for c in abundance.columns if c != group_column]

    if len(groups) < 2:
        return pd.DataFrame()

    group_labels = merged[group_column].values
    abundance_matrix = merged[taxa_cols].values

    rows: list[dict] = []
    rng = np.random.default_rng(seed)

    for t_idx, taxon in enumerate(taxa_cols):
        vals = abundance_matrix[:, t_idx]
        best_indval = 0.0
        best_group = groups[0]

        for group in groups:
            mask = group_labels == group
            n_group = mask.sum()
            if n_group == 0:
                continue
            # Specificity: mean abundance in group / sum of mean abundances across groups
            group_mean = vals[mask].mean()
            total_mean = sum(vals[group_labels == g].mean() for g in groups)
            specificity = group_mean / total_mean if total_mean > 0 else 0.0
            # Fidelity: proportion of group samples where taxon is present
            fidelity = (vals[mask] > 0).sum() / n_group
            indval = specificity * fidelity
            if indval > best_indval:
                best_indval = indval
                best_group = group

        # Permutation test
        count_ge = 0
        for _ in range(permutations):
            perm_labels = rng.permutation(group_labels)
            for group in groups:
                mask = perm_labels == group
                n_group = mask.sum()
                if n_group == 0:
                    continue
                group_mean = vals[mask].mean()
                total_mean = sum(vals[perm_labels == g].mean() for g in groups)
                specificity = group_mean / total_mean if total_mean > 0 else 0.0
                fidelity = (vals[mask] > 0).sum() / n_group
                if specificity * fidelity >= best_indval:
                    count_ge += 1
                    break

        p_val = (count_ge + 1) / (permutations + 1)
        rows.append({
            "taxon": taxon,
            "indicator_group": str(best_group),
            "indval": float(best_indval),
            "p_value": float(p_val),
        })

    result = pd.DataFrame(rows).sort_values("indval", ascending=False).reset_index(drop=True)
    return result


def hedges_g(group1: np.ndarray, group2: np.ndarray) -> tuple[float, float]:
    """Compute Hedges' g effect size and its standard error."""
    n1, n2 = len(group1), len(group2)
    if n1 < 2 or n2 < 2:
        return np.nan, np.nan
    mean1, mean2 = group1.mean(), group2.mean()
    var1, var2 = group1.var(ddof=1), group2.var(ddof=1)
    pooled_sd = math.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    if pooled_sd <= 0:
        return np.nan, np.nan
    d = (mean1 - mean2) / pooled_sd
    # Small-sample correction factor
    correction = 1.0 - 3.0 / (4.0 * (n1 + n2 - 2) - 1.0)
    g = d * correction
    se = math.sqrt((n1 + n2) / (n1 * n2) + g**2 / (2.0 * (n1 + n2)))
    return float(g), float(se)


def random_effects_meta(
    effects: np.ndarray,
    se: np.ndarray,
) -> dict:
    """DerSimonian-Laird random-effects meta-analysis."""
    valid = np.isfinite(effects) & np.isfinite(se) & (se > 0)
    effects = effects[valid]
    se = se[valid]
    k = len(effects)
    if k == 0:
        return {"pooled_effect": np.nan, "pooled_se": np.nan, "ci_lower": np.nan,
                "ci_upper": np.nan, "Q": np.nan, "I_squared": np.nan, "tau_squared": np.nan, "k": 0}

    weights_fe = 1.0 / se**2
    pooled_fe = np.sum(weights_fe * effects) / np.sum(weights_fe)
    Q = np.sum(weights_fe * (effects - pooled_fe) ** 2)
    df = k - 1
    C = np.sum(weights_fe) - np.sum(weights_fe**2) / np.sum(weights_fe)

    tau_sq = max(0.0, (Q - df) / C) if C > 0 and df > 0 else 0.0
    I_sq = max(0.0, (Q - df) / Q * 100.0) if Q > 0 and df > 0 else 0.0

    weights_re = 1.0 / (se**2 + tau_sq)
    pooled_re = float(np.sum(weights_re * effects) / np.sum(weights_re))
    pooled_se_re = float(math.sqrt(1.0 / np.sum(weights_re)))
    ci_lower = pooled_re - 1.96 * pooled_se_re
    ci_upper = pooled_re + 1.96 * pooled_se_re

    p_value = float(2.0 * stats.norm.sf(abs(pooled_re / pooled_se_re))) if pooled_se_re > 0 else np.nan

    return {
        "pooled_effect": pooled_re,
        "pooled_se": pooled_se_re,
        "ci_lower": float(ci_lower),
        "ci_upper": float(ci_upper),
        "z_value": float(pooled_re / pooled_se_re) if pooled_se_re > 0 else np.nan,
        "p_value": p_value,
        "Q": float(Q),
        "I_squared": float(I_sq),
        "tau_squared": float(tau_sq),
        "k": int(k),
    }


def per_study_effect_sizes(
    alpha: pd.DataFrame,
    group_column: str,
    study_column: str,
    metric: str = "shannon",
    group1: str = "Seeds",
    group2: str = "Seedlings",
) -> pd.DataFrame:
    """Compute Hedges' g for each study comparing two groups on an alpha metric."""
    if any(c not in alpha.columns for c in [group_column, study_column, metric]):
        return pd.DataFrame()

    subset = alpha[alpha[group_column].isin([group1, group2])].dropna(subset=[metric])
    studies = sorted(subset[study_column].dropna().unique())

    rows: list[dict] = []
    for study in studies:
        study_data = subset[subset[study_column] == study]
        g1_vals = study_data.loc[study_data[group_column] == group1, metric].values
        g2_vals = study_data.loc[study_data[group_column] == group2, metric].values
        if len(g1_vals) < 2 or len(g2_vals) < 2:
            continue
        g, se = hedges_g(g1_vals, g2_vals)
        rows.append({
            "study": str(study),
            "n_group1": len(g1_vals),
            "n_group2": len(g2_vals),
            "mean_group1": float(g1_vals.mean()),
            "mean_group2": float(g2_vals.mean()),
            "hedges_g": g,
            "se": se,
            "ci_lower": g - 1.96 * se if np.isfinite(se) else np.nan,
            "ci_upper": g + 1.96 * se if np.isfinite(se) else np.nan,
        })
    return pd.DataFrame(rows)


def top_taxa_frame(relative_frame: pd.DataFrame, top_n: int, exclude_taxa: set[str] | None = None) -> pd.DataFrame:
    exclude_taxa = exclude_taxa or set()
    usable = [column for column in relative_frame.columns if column not in exclude_taxa]
    if not usable:
        return pd.DataFrame(index=relative_frame.index)
    totals = relative_frame[usable].sum(axis=0).sort_values(ascending=False)
    return relative_frame[totals.head(top_n).index]


def write_distance_matrix(distance_matrix: np.ndarray, sample_ids: list[str], output_path: Path) -> None:
    frame = pd.DataFrame(distance_matrix, index=sample_ids, columns=sample_ids)
    frame.index.name = "run_accession"
    frame.to_csv(output_path)


def plotting_libs():
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    import matplotlib.pyplot as plt
    import seaborn as sns

    return plt, sns


def plot_library_sizes(sample_qc: pd.DataFrame, output_path: Path) -> None:
    plt, sns = plotting_libs()
    plt.figure(figsize=(9, 5))
    sns.histplot(sample_qc["library_size"], bins=30)
    plt.xlabel("Library Size")
    plt.ylabel("Samples")
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def plot_alpha_boxplots(alpha: pd.DataFrame, group_column: str, output_path: Path) -> None:
    if group_column not in alpha.columns:
        return
    plt, sns = plotting_libs()
    subset = alpha[[group_column, "shannon", "observed_features"]].dropna()
    if subset.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    sns.boxplot(data=subset, x=group_column, y="shannon", ax=axes[0])
    axes[0].tick_params(axis="x", rotation=30)
    sns.boxplot(data=subset, x=group_column, y="observed_features", ax=axes[1])
    axes[1].tick_params(axis="x", rotation=30)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_pcoa(coords: pd.DataFrame, metadata: pd.DataFrame, group_column: str, explained: np.ndarray, output_path: Path) -> None:
    plt, sns = plotting_libs()
    merged = coords.join(metadata.set_index("run_accession")[[group_column]], how="left")
    plt.figure(figsize=(7, 6))
    sns.scatterplot(data=merged, x="PC1", y="PC2", hue=group_column, s=50)
    x_label = "PC1"
    y_label = "PC2"
    if explained.size >= 2:
        x_label = f"PC1 ({explained[0] * 100:.1f}%)"
        y_label = f"PC2 ({explained[1] * 100:.1f}%)"
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def plot_unassigned_fraction(sample_qc: pd.DataFrame, group_column: str, output_path: Path) -> None:
    if group_column not in sample_qc.columns:
        return
    plt, sns = plotting_libs()
    plt.figure(figsize=(9, 5))
    sns.boxplot(data=sample_qc, x=group_column, y="unassigned_fraction")
    plt.xticks(rotation=30)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def plot_taxa_stack(relative_frame: pd.DataFrame, metadata: pd.DataFrame, group_column: str, output_path: Path) -> None:
    if relative_frame.empty or group_column not in metadata.columns:
        return
    plt, _ = plotting_libs()
    merged = relative_frame.join(metadata.set_index("run_accession")[[group_column]], how="left")
    grouped = merged.groupby(group_column, dropna=False).mean(numeric_only=True)
    grouped.plot(kind="bar", stacked=True, figsize=(10, 6), colormap="tab20")
    plt.ylabel("Mean Relative Abundance")
    plt.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def plot_alpha_significance(
    alpha: pd.DataFrame,
    group_column: str,
    kw_results: pd.DataFrame,
    dunn_results: pd.DataFrame,
    output_path: Path,
) -> None:
    """Alpha diversity boxplots with significance brackets from Dunn's test."""
    if group_column not in alpha.columns:
        return
    plt, sns = plotting_libs()
    metrics = ["shannon", "observed_features", "simpson", "pielou_evenness"]
    metrics = [m for m in metrics if m in alpha.columns]
    n_metrics = len(metrics)
    if n_metrics == 0:
        return

    fig, axes = plt.subplots(1, n_metrics, figsize=(4 * n_metrics, 6))
    if n_metrics == 1:
        axes = [axes]

    for ax, metric in zip(axes, metrics):
        subset = alpha[[group_column, metric]].dropna()
        sns.boxplot(data=subset, x=group_column, y=metric, ax=ax, palette="Set2")
        sns.stripplot(data=subset, x=group_column, y=metric, ax=ax, color="0.3", size=2, alpha=0.4)
        ax.set_title(metric.replace("_", " ").title())
        ax.tick_params(axis="x", rotation=30)

        # Add significance brackets for significant pairs
        if not dunn_results.empty:
            metric_dunns = dunn_results if metric == "shannon" else pd.DataFrame()
            sig_pairs = metric_dunns[metric_dunns["p_adjusted"] < 0.05] if not metric_dunns.empty else pd.DataFrame()
            if not sig_pairs.empty:
                groups_order = [t.get_text() for t in ax.get_xticklabels()]
                y_max = subset[metric].max()
                y_range = subset[metric].max() - subset[metric].min()
                for idx, row in sig_pairs.head(3).iterrows():
                    if row["group1"] in groups_order and row["group2"] in groups_order:
                        x1 = groups_order.index(row["group1"])
                        x2 = groups_order.index(row["group2"])
                        y = y_max + y_range * 0.05 * (idx + 1)
                        ax.plot([x1, x1, x2, x2], [y, y + y_range * 0.01, y + y_range * 0.01, y], color="k", lw=1)
                        p_text = f"p={row['p_adjusted']:.3f}" if row["p_adjusted"] >= 0.001 else "p<0.001"
                        ax.text((x1 + x2) / 2, y + y_range * 0.01, p_text, ha="center", va="bottom", fontsize=7)

        # Add KW p-value to title
        if not kw_results.empty:
            kw_row = kw_results[kw_results["metric"] == metric]
            if not kw_row.empty:
                p = kw_row.iloc[0]["p_value"]
                p_text = f"p={p:.3f}" if p >= 0.001 else "p<0.001"
                ax.set_title(f"{metric.replace('_', ' ').title()}\nKW {p_text}")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_pcoa_ellipses(
    coords: pd.DataFrame,
    metadata: pd.DataFrame,
    group_column: str,
    explained: np.ndarray,
    output_path: Path,
) -> None:
    """PCoA scatter with 95% confidence ellipses per group."""
    plt, sns = plotting_libs()
    from matplotlib.patches import Ellipse

    merged = coords.join(metadata.set_index("run_accession")[[group_column]], how="left")
    fig, ax = plt.subplots(figsize=(8, 7))
    palette = sns.color_palette("Set2", n_colors=merged[group_column].nunique())
    groups = sorted(merged[group_column].dropna().unique())

    for idx, group in enumerate(groups):
        group_data = merged[merged[group_column] == group]
        color = palette[idx % len(palette)]
        ax.scatter(group_data["PC1"], group_data["PC2"], color=color, s=40, alpha=0.6, label=group, edgecolors="white", linewidth=0.5)

        if len(group_data) >= 3:
            mean_x, mean_y = group_data["PC1"].mean(), group_data["PC2"].mean()
            cov = np.cov(group_data["PC1"].values, group_data["PC2"].values)
            eigenvalues, eigenvectors = np.linalg.eigh(cov)
            order = eigenvalues.argsort()[::-1]
            eigenvalues = eigenvalues[order]
            eigenvectors = eigenvectors[:, order]
            angle = np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))
            width, height = 2 * np.sqrt(eigenvalues * 5.991)  # 95% CI with chi2(2)
            ellipse = Ellipse(
                xy=(mean_x, mean_y), width=width, height=height, angle=angle,
                facecolor=color, alpha=0.15, edgecolor=color, linewidth=1.5,
            )
            ax.add_patch(ellipse)

    x_label = "PC1"
    y_label = "PC2"
    if explained.size >= 2:
        x_label = f"PC1 ({explained[0] * 100:.1f}%)"
        y_label = f"PC2 ({explained[1] * 100:.1f}%)"
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_heatmap(
    abundance: pd.DataFrame,
    metadata: pd.DataFrame,
    group_column: str,
    top_n: int,
    output_path: Path,
) -> None:
    """Heatmap of top taxa grouped by condition."""
    plt, sns = plotting_libs()
    merged = abundance.join(metadata.set_index("run_accession")[[group_column]], how="inner")
    grouped = merged.groupby(group_column, dropna=False).mean(numeric_only=True)
    top_cols = grouped.sum().sort_values(ascending=False).head(top_n).index
    data = grouped[top_cols]

    fig, ax = plt.subplots(figsize=(max(8, top_n * 0.5), max(4, len(data) * 0.8)))
    sns.heatmap(data, cmap="YlOrRd", annot=top_n <= 20, fmt=".3f", ax=ax, linewidths=0.5)
    ax.set_ylabel(group_column)
    ax.set_xlabel("Taxon")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_volcano(da_results: pd.DataFrame, output_path: Path) -> None:
    """Volcano plot for differential abundance results."""
    if da_results.empty:
        return
    plt, _ = plotting_libs()
    fig, ax = plt.subplots(figsize=(8, 6))

    neg_log10_p = -np.log10(da_results["p_adjusted"].clip(lower=1e-300))
    lfc = da_results["log2_fold_change"]

    sig_mask = (da_results["p_adjusted"] < 0.05) & (abs(lfc) > 1)
    ax.scatter(lfc[~sig_mask], neg_log10_p[~sig_mask], c="grey", alpha=0.5, s=20, label="NS")
    ax.scatter(lfc[sig_mask], neg_log10_p[sig_mask], c="red", alpha=0.7, s=30, label="Significant")

    ax.axhline(-np.log10(0.05), color="grey", linestyle="--", linewidth=0.8)
    ax.axvline(-1, color="grey", linestyle="--", linewidth=0.8)
    ax.axvline(1, color="grey", linestyle="--", linewidth=0.8)

    # Label top significant taxa
    top_sig = da_results[sig_mask].nsmallest(10, "p_adjusted")
    for _, row in top_sig.iterrows():
        ax.annotate(
            row["taxon"], (row["log2_fold_change"], -np.log10(max(row["p_adjusted"], 1e-300))),
            fontsize=6, alpha=0.8, ha="center",
        )

    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10(adjusted p-value)")
    ax.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_forest(effect_sizes: pd.DataFrame, meta_result: dict, output_path: Path) -> None:
    """Forest plot for meta-analysis results."""
    if effect_sizes.empty:
        return
    plt, _ = plotting_libs()
    n_studies = len(effect_sizes)
    fig, ax = plt.subplots(figsize=(8, max(4, n_studies * 0.6 + 2)))

    y_positions = list(range(n_studies, 0, -1))

    for i, (_, row) in enumerate(effect_sizes.iterrows()):
        y = y_positions[i]
        ax.plot([row["ci_lower"], row["ci_upper"]], [y, y], color="black", linewidth=1)
        ax.plot(row["hedges_g"], y, "s", color="steelblue", markersize=8)
        ax.text(
            ax.get_xlim()[0] if i == 0 else effect_sizes["ci_lower"].min() - 0.5,
            y, f"{row['study']} (n={row['n_group1']}+{row['n_group2']})",
            va="center", ha="right", fontsize=8,
        )

    # Add pooled estimate diamond
    pooled = meta_result.get("pooled_effect", np.nan)
    ci_lo = meta_result.get("ci_lower", np.nan)
    ci_hi = meta_result.get("ci_upper", np.nan)
    if np.isfinite(pooled):
        y = 0
        diamond_height = 0.3
        ax.fill(
            [ci_lo, pooled, ci_hi, pooled],
            [y, y + diamond_height, y, y - diamond_height],
            color="red", alpha=0.6,
        )
        ax.text(
            effect_sizes["ci_lower"].min() - 0.5, y, "Pooled RE",
            va="center", ha="right", fontsize=8, fontweight="bold",
        )

    ax.axvline(0, color="grey", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Hedges' g")
    ax.set_yticks([])

    i_sq = meta_result.get("I_squared", np.nan)
    p_val = meta_result.get("p_value", np.nan)
    title = "Forest Plot: Seeds vs Seedlings (Shannon)"
    if np.isfinite(i_sq):
        title += f"\nI²={i_sq:.1f}%"
    if np.isfinite(p_val):
        title += f", p={'<0.001' if p_val < 0.001 else f'{p_val:.3f}'}"
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close(fig)


# ── R bridge functions ───────────────────────────────────────────────────────


def export_filtered_feature_ids(feature_ids: list[str], output_path: Path) -> Path:
    """Write selected feature IDs to CSV for R to consume."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"OTUID": feature_ids}).to_csv(output_path, index=False)
    return output_path


def run_phyloseq_analysis(
    config_path: Path,
    filtered_features_path: Path,
    output_dir: Path,
) -> dict:
    """Run phyloseq_analysis.R as a subprocess, returning status dict."""
    if shutil.which("Rscript") is None:
        print("Rscript not found on PATH — skipping R phylogenetic analysis")
        return {"status": "skipped", "reason": "rscript_not_found"}

    r_script = Path(__file__).resolve().parent / "phyloseq_analysis.R"
    if not r_script.exists():
        print(f"R script not found at {r_script} — skipping")
        return {"status": "skipped", "reason": "r_script_missing"}

    cmd = [
        "Rscript", "--vanilla", str(r_script),
        "--config", str(config_path),
    ]
    if filtered_features_path.exists():
        cmd.extend(["--filtered-features", str(filtered_features_path)])

    print(f"Running: {' '.join(cmd)}")
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600,
            cwd=str(Path(__file__).resolve().parent),
        )
    except subprocess.TimeoutExpired:
        return {"status": "failed", "reason": "timeout"}

    if result.stderr:
        for line in result.stderr.strip().splitlines()[-10:]:
            print(f"  R: {line}")

    # Parse JSON status from the last stdout line containing __R_STATUS__
    r_status: dict = {}
    for line in reversed(result.stdout.splitlines()):
        if "__R_STATUS__" in line:
            json_str = line.split("__R_STATUS__", 1)[1].strip()
            try:
                r_status = json.loads(json_str)
            except json.JSONDecodeError:
                pass
            break

    if result.returncode != 0 and not r_status:
        return {"status": "failed", "reason": "nonzero_exit", "returncode": result.returncode}

    r_status["status"] = "completed"
    return r_status


def integrate_r_outputs(output_dir: Path, alpha: pd.DataFrame) -> pd.DataFrame:
    """Merge Faith's PD from R output into the alpha diversity DataFrame."""
    faith_path = output_dir / "diversity" / "faith_pd.csv"
    if not faith_path.exists():
        return alpha

    faith = pd.read_csv(faith_path)
    if "sample_id" not in faith.columns or "faith_pd" not in faith.columns:
        return alpha

    faith_lookup = faith.set_index("sample_id")["faith_pd"]
    alpha = alpha.copy()
    alpha["faith_pd"] = alpha["run_accession"].map(faith_lookup)
    return alpha


def run_pipeline(config_path: Path) -> None:
    config, base_dir = load_config(config_path)
    counts_path = resolve_path(base_dir, config["inputs"]["counts"])
    taxonomy_path = resolve_path(base_dir, config["inputs"]["taxonomy"])
    metadata_path = resolve_path(base_dir, config["inputs"]["metadata"])
    tree_path = resolve_path(base_dir, config["inputs"]["tree"])
    output_dir = resolve_path(base_dir, config["outputs"]["directory"])
    generate_plots = bool(config["outputs"].get("generate_plots", True))
    qc_dir = output_dir / "qc"
    diversity_dir = output_dir / "diversity"
    taxonomy_dir = output_dir / "taxonomy"
    derived_dir = output_dir / "derived"
    da_dir = output_dir / "differential_abundance"
    meta_dir = output_dir / "meta_analysis"
    figures_dir = output_dir / "figures"
    for path in (qc_dir, diversity_dir, taxonomy_dir, derived_dir, da_dir, meta_dir, figures_dir):
        path.mkdir(parents=True, exist_ok=True)

    ensured_nexus = ensure_tree_nexus(tree_path, derived_dir / "tree.nexus")
    tree_tips = extract_tree_tips(ensured_nexus)

    metadata = load_metadata(metadata_path)
    metadata = extract_timepoint(metadata)
    taxonomy = load_taxonomy(taxonomy_path)

    # ── Taxonomy-based feature mask ──────────────────────────────────
    # Read just the header to get feature IDs for mask construction
    with counts_path.open("r", newline="") as handle:
        header_reader = csv.reader(handle)
        header_row = next(header_reader)
    all_feature_ids = [item.strip() for item in header_row[1:]]

    taxonomy_mask = build_taxonomy_mask(all_feature_ids, taxonomy, config["filters"])
    n_total_features = len(all_feature_ids)
    n_pass_taxonomy = int(taxonomy_mask.sum())
    n_excluded_kingdom = 0
    n_excluded_chloroplast = 0
    n_excluded_mitochondria = 0

    keep_kingdoms = set(config["filters"].get("keep_kingdom", []))
    exclude_orders = set(config["filters"].get("exclude_order", []))
    exclude_families = set(config["filters"].get("exclude_family", []))

    for fid in all_feature_ids:
        if fid not in taxonomy.index:
            if keep_kingdoms:
                n_excluded_kingdom += 1
            continue
        row = taxonomy.loc[fid]
        kingdom = row.get("Kingdom", "")
        if keep_kingdoms and kingdom not in keep_kingdoms:
            n_excluded_kingdom += 1
        elif exclude_orders and row.get("Order", "") in exclude_orders:
            n_excluded_chloroplast += 1
        elif exclude_families and row.get("Family", "") in exclude_families:
            n_excluded_mitochondria += 1

    print(f"Taxonomy filter: {n_total_features} total ASVs")
    print(f"  Kept (Bacteria, non-chloroplast, non-mitochondria): {n_pass_taxonomy}")
    print(f"  Excluded by Kingdom filter: {n_excluded_kingdom}")
    print(f"  Excluded chloroplasts (Order): {n_excluded_chloroplast}")
    print(f"  Excluded mitochondria (Family): {n_excluded_mitochondria}")

    first_pass = first_pass_counts(
        counts_path,
        max_samples=config["analysis"].get("max_samples"),
        taxonomy_mask=taxonomy_mask,
    )
    sample_ids = first_pass["sample_ids"]
    metadata_ids = set(metadata["run_accession"])
    count_ids = set(sample_ids)

    selected_indices = select_feature_indices(first_pass, config["filters"], taxonomy_mask)
    second_pass = second_pass_collect(
        counts_path=counts_path,
        sample_ids=sample_ids,
        feature_ids=first_pass["feature_ids"],
        selected_indices=selected_indices,
        taxonomy=taxonomy,
        taxonomy_mask=taxonomy_mask,
    )

    # ── QC summary: raw metrics from all taxonomy-filtered features, ALL samples ──
    sample_qc = pd.DataFrame(
        {
            "run_accession": sample_ids,
            "library_size": first_pass["library_sizes"],
            "observed_features": first_pass["richness"],
            "shannon": first_pass["shannon"],
            "simpson": first_pass["simpson"],
            "pielou_evenness": [safe_pielou(h, r) for h, r in zip(first_pass["shannon"], first_pass["richness"])],
            "unassigned_fraction": second_pass["unassigned_fraction"],
        }
    )
    min_lib = int(config["filters"]["min_library_size"])
    sample_qc["passes_library_filter"] = sample_qc["library_size"] >= min_lib
    sample_qc = sample_qc.merge(metadata, on="run_accession", how="left")
    sample_qc.to_csv(qc_dir / "sample_qc_summary.csv", index=False)

    feature_filter_summary = pd.DataFrame(
        [
            {
                "input_features": n_total_features,
                "pass_taxonomy_filter": n_pass_taxonomy,
                "excluded_non_bacteria": n_excluded_kingdom,
                "excluded_chloroplast": n_excluded_chloroplast,
                "excluded_mitochondria": n_excluded_mitochondria,
                "retained_features_alpha": int(len(selected_indices)),
                "min_feature_prevalence": int(config["filters"]["min_feature_prevalence"]),
                "min_feature_total_count": int(config["filters"]["min_feature_total_count"]),
                "max_retained_features": int(config["filters"]["max_retained_features"]),
                "min_library_size": min_lib,
                "blank_cells_in_counts": int(first_pass["blank_cells"]),
                "zero_library_samples": int(first_pass["zero_rows"]),
                "metadata_rows": int(len(metadata)),
                "count_table_rows": int(len(sample_ids)),
                "metadata_only_ids": int(len(metadata_ids - count_ids)),
                "count_only_ids": int(len(count_ids - metadata_ids)),
                "tree_tip_count": int(len(tree_tips)),
                "tree_selected_feature_overlap": int(len(tree_tips & set(second_pass["selected_feature_ids"]))),
                "unassigned_taxa_rows": int((taxonomy["Kingdom"] == "Unassigned").sum()),
            }
        ]
    )
    feature_filter_summary.to_csv(qc_dir / "feature_filter_summary.csv", index=False)

    # ── Alpha diversity: feature-filtered metrics, sample-filtered by min_library_size ──
    # Compute alpha from the prevalence/abundance-filtered feature matrix
    sel_mat = second_pass["selected_matrix"]
    n_samp = len(sample_ids)
    alpha_observed = np.array([(sel_mat[i] > 0).sum() for i in range(n_samp)], dtype=int)
    alpha_shannon_vals = np.array([shannon_index(sel_mat[i]) for i in range(n_samp)])
    alpha_simpson_vals = np.array([simpson_index(sel_mat[i]) for i in range(n_samp)])
    alpha_pielou_vals = np.array(
        [safe_pielou(h, int(r)) for h, r in zip(alpha_shannon_vals, alpha_observed)]
    )

    alpha_all = pd.DataFrame(
        {
            "run_accession": sample_ids,
            "library_size": first_pass["library_sizes"],
            "observed_features": alpha_observed,
            "shannon": alpha_shannon_vals,
            "simpson": alpha_simpson_vals,
            "pielou_evenness": alpha_pielou_vals,
            "unassigned_fraction": second_pass["unassigned_fraction"],
        }
    )

    # Apply min_library_size threshold (only affects alpha diversity output)
    alpha = alpha_all[alpha_all["library_size"] >= min_lib].copy()
    meta_cols_alpha = [c for c in ["Developmental_Stage", "study_alias", "Scientific_Name", "Geography", "Timepoint"] if c in metadata.columns]
    if meta_cols_alpha:
        meta_lookup = metadata.set_index("run_accession")[meta_cols_alpha]
        alpha = alpha.join(meta_lookup, on="run_accession")
    alpha.to_csv(diversity_dir / "alpha_diversity.csv", index=False)

    n_alpha_samples = len(alpha)
    n_excluded_by_depth = n_samp - n_alpha_samples
    print(f"Alpha diversity: {len(selected_indices)} features, "
          f"{n_alpha_samples} samples (excluded {n_excluded_by_depth} below {min_lib} reads)")

    genus_relative = second_pass["relative_rank_frames"]["Genus"]
    genus_relative.to_csv(taxonomy_dir / "genus_relative_abundance.csv")
    second_pass["relative_rank_frames"]["Phyla"].to_csv(taxonomy_dir / "phylum_relative_abundance.csv")
    second_pass["relative_rank_frames"]["Family"].to_csv(taxonomy_dir / "family_relative_abundance.csv")

    beta_rank = config["analysis"].get("beta_rank", "Genus")
    beta_source = second_pass["relative_rank_frames"][beta_rank]
    beta_frame = top_taxa_frame(
        beta_source,
        top_n=int(config["filters"]["beta_top_taxa"]),
        exclude_taxa={UNASSIGNED_BY_RANK.get(beta_rank, "unassigned")},
    )
    if beta_frame.empty:
        beta_frame = beta_source.iloc[:, : min(10, beta_source.shape[1])]
    beta_matrix = beta_frame.to_numpy(dtype=np.float32)

    bray = bray_curtis_distance(beta_matrix)
    jaccard = jaccard_distance(beta_matrix)
    write_distance_matrix(bray, sample_ids, diversity_dir / "braycurtis_distance.csv")
    write_distance_matrix(jaccard, sample_ids, diversity_dir / "jaccard_distance.csv")

    coords, explained = pcoa(bray, axes=int(config["analysis"].get("pcoa_axes", 2)))
    coords.insert(0, "run_accession", sample_ids)
    coords.to_csv(diversity_dir / "braycurtis_pcoa.csv", index=False)

    # NMDS ordination
    nmds_coords = nmds(bray, dimensions=2)
    nmds_coords.insert(0, "run_accession", sample_ids)
    nmds_coords.to_csv(diversity_dir / "nmds_coordinates.csv", index=False)

    grouping_columns = [config["analysis"]["primary_group"]] + list(config["analysis"].get("secondary_groups", []))
    permanova_rows = []
    permdisp_rows = []
    merged_meta = metadata.set_index("run_accession").reindex(sample_ids)
    permutations = int(config["analysis"].get("permutations", 99))
    for column in grouping_columns:
        if column not in merged_meta.columns:
            continue
        labels = merged_meta[column].fillna("missing").astype(str).tolist()
        permanova_result = permanova(bray, labels, permutations=permutations)
        permanova_rows.append({"metric": "braycurtis", "group_column": column, **permanova_result})
        permdisp_result = permdisp(bray, labels, permutations=permutations)
        permdisp_rows.append({"metric": "braycurtis", "group_column": column, **permdisp_result})

    pd.DataFrame(permanova_rows).to_csv(diversity_dir / "permanova_results.csv", index=False)
    pd.DataFrame(permdisp_rows).to_csv(diversity_dir / "permdisp_results.csv", index=False)

    # Kruskal-Wallis and Dunn's post-hoc tests
    primary_group = config["analysis"]["primary_group"]
    kw_results = kruskal_wallis_tests(alpha, primary_group)
    kw_results.to_csv(diversity_dir / "kruskal_wallis_results.csv", index=False)

    dunn_results = pd.DataFrame()
    for metric in ["shannon", "observed_features", "simpson", "pielou_evenness"]:
        if metric in alpha.columns:
            metric_dunns = dunns_posthoc(alpha, primary_group, metric=metric)
            if not metric_dunns.empty:
                metric_dunns.insert(0, "metric", metric)
                dunn_results = pd.concat([dunn_results, metric_dunns], ignore_index=True)
    dunn_results.to_csv(diversity_dir / "dunns_posthoc_results.csv", index=False)

    # Differential abundance analysis
    da_config = config.get("differential_abundance", {})
    da_groups = da_config.get("groups", ["Seeds", "Seedlings"])
    da_min_prevalence = float(da_config.get("min_prevalence", 0.1))
    da_results = pd.DataFrame()
    if len(da_groups) >= 2 and primary_group in alpha.columns:
        da_results = differential_abundance(
            genus_relative, metadata, primary_group,
            group1=da_groups[0], group2=da_groups[1],
            min_prevalence=da_min_prevalence,
        )
        da_results.to_csv(da_dir / "da_results.csv", index=False)

    # Core microbiome analysis
    core_threshold = float(config.get("taxonomy_analysis", {}).get("core_prevalence_threshold", 0.5))
    core_abundance_threshold = float(config.get("taxonomy_analysis", {}).get("core_abundance_threshold", 0.001))
    core_results = core_microbiome(
        genus_relative, metadata, primary_group,
        prevalence_threshold=core_threshold,
        abundance_threshold=core_abundance_threshold,
    )
    core_results.to_csv(taxonomy_dir / "core_microbiome.csv", index=False)

    # Indicator species analysis
    indval_perms = int(config.get("taxonomy_analysis", {}).get("indicator_permutations", 99))
    indval_results = indicator_species(
        genus_relative, metadata, primary_group,
        permutations=indval_perms,
    )
    indval_results.to_csv(taxonomy_dir / "indicator_species.csv", index=False)

    # Meta-analysis: per-study effect sizes and random-effects pooling
    study_column = "study_alias"
    meta_effect_sizes = pd.DataFrame()
    meta_pooled = {}
    if study_column in alpha.columns and primary_group in alpha.columns and len(da_groups) >= 2:
        meta_effect_sizes = per_study_effect_sizes(
            alpha, primary_group, study_column,
            metric="shannon", group1=da_groups[0], group2=da_groups[1],
        )
        meta_effect_sizes.to_csv(meta_dir / "per_study_effect_sizes.csv", index=False)

        if len(meta_effect_sizes) >= 2:
            meta_pooled = random_effects_meta(
                meta_effect_sizes["hedges_g"].values,
                meta_effect_sizes["se"].values,
            )
            pd.DataFrame([meta_pooled]).to_csv(meta_dir / "pooled_meta_analysis.csv", index=False)

    # Plotting
    if generate_plots:
        plot_library_sizes(sample_qc, qc_dir / "library_size_histogram.png")
        plot_alpha_boxplots(sample_qc, primary_group, diversity_dir / "alpha_by_primary_group.png")
        plot_pcoa(coords.set_index("run_accession"), metadata, primary_group, explained, diversity_dir / "braycurtis_pcoa.png")
        plot_unassigned_fraction(sample_qc, primary_group, qc_dir / "unassigned_fraction_by_primary_group.png")
        plot_taxa_stack(
            top_taxa_frame(
                second_pass["relative_rank_frames"]["Genus"],
                top_n=10,
                exclude_taxa={UNASSIGNED_BY_RANK["Genus"]},
            ),
            metadata,
            primary_group,
            taxonomy_dir / "top_genus_by_primary_group.png",
        )

        # Enhanced visualizations
        plot_alpha_significance(
            alpha, primary_group, kw_results,
            dunn_results[dunn_results["metric"] == "shannon"] if not dunn_results.empty else pd.DataFrame(),
            figures_dir / "alpha_significance.png",
        )
        plot_pcoa_ellipses(
            coords.set_index("run_accession"), metadata, primary_group, explained,
            figures_dir / "pcoa_ellipses.png",
        )
        plot_heatmap(
            genus_relative, metadata, primary_group,
            top_n=20, output_path=figures_dir / "genus_heatmap.png",
        )
        if not da_results.empty:
            plot_volcano(da_results, figures_dir / "volcano_plot.png")
        if not meta_effect_sizes.empty and meta_pooled:
            plot_forest(meta_effect_sizes, meta_pooled, figures_dir / "forest_plot.png")

    # ── R phylogenetic analysis bridge ──────────────────────────────────
    filtered_ids_path = qc_dir / "filtered_feature_ids.csv"
    export_filtered_feature_ids(second_pass["selected_feature_ids"], filtered_ids_path)

    r_status = run_phyloseq_analysis(
        config_path.resolve(), filtered_ids_path, output_dir,
    )
    print(f"R phyloseq status: {r_status.get('status', 'unknown')}")

    faith_pd_status = "skipped_missing_skbio"
    unifrac_status = "skipped_missing_skbio"

    if r_status.get("status") == "completed":
        faith_sub = r_status.get("faith_pd", {})
        unifrac_sub = r_status.get("unifrac", {})
        faith_pd_status = faith_sub.get("status", "unknown") if isinstance(faith_sub, dict) else str(faith_sub)
        unifrac_status = unifrac_sub.get("status", "unknown") if isinstance(unifrac_sub, dict) else str(unifrac_sub)

        # Merge Faith's PD into alpha diversity
        alpha = integrate_r_outputs(output_dir, alpha)
        alpha.to_csv(diversity_dir / "alpha_diversity.csv", index=False)

    phylo_summary = pd.DataFrame(
        [
            {
                "tree_path": str(ensured_nexus),
                "tree_tip_count": len(tree_tips),
                "selected_feature_overlap": len(
                    tree_tips & set(second_pass["selected_feature_ids"])
                ),
                "faith_pd_status": faith_pd_status,
                "unifrac_status": unifrac_status,
                "r_script_status": r_status.get("status", "unknown"),
            }
        ]
    )
    phylo_summary.to_csv(diversity_dir / "phylogenetic_metrics_status.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the seedling microbiome analysis workflow.")
    parser.add_argument("--config", default="analysis_config.json", help="Path to the analysis config JSON file.")
    args = parser.parse_args()
    run_pipeline(Path(args.config))


if __name__ == "__main__":
    main()
