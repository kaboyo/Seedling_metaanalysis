import csv
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


class MicrobiomeAnalysisTests(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.base = Path(self.tmpdir.name)
        self.datasets = self.base / "Datesets"
        self.datasets.mkdir()

    def tearDown(self):
        self.tmpdir.cleanup()

    def write_text(self, relative_path: str, content: str, encoding: str = "utf-8") -> Path:
        path = self.base / relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding=encoding)
        return path

    def _write_metadata(self, encoding: str) -> None:
        self.write_text(
            "Datesets/metadata.csv",
            (
                "run_accession,,Author,Geography,Cultivar,Developmental_Stage,"
                "Scientific_Name,sample_description,Title,study_alias,command_to_download_the_reads\n"
                "S1,,Acu\xf1a,Chile,c1,Seeds,Plant A,desc,title,Study1,cmd\n"
                "S2,,Acu\xf1a,Chile,c1,Seedlings,Plant A,desc,title,Study1,cmd\n"
            ),
            encoding=encoding,
        )

    def test_load_metadata_decodes_mac_roman_by_default(self):
        """The project's own Datesets/metadata.csv is Mac Roman, not latin-1.

        The giveaway byte is 0x96: Mac Roman decodes it to 'n-tilde', latin-1 maps it
        to an unassigned control character. latin-1 never raises, so trying it first
        silently mangled 33 author names (Acuna -> Acua) rather than failing loudly.
        """
        self._write_metadata("mac_roman")

        from microbiome_analysis import load_metadata

        metadata = load_metadata(self.base / "Datesets" / "metadata.csv")
        self.assertEqual(list(metadata["run_accession"]), ["S1", "S2"])
        self.assertIn("Acu\xf1a", set(metadata["Author"]))

    def test_load_metadata_honours_explicit_encoding_order(self):
        """latin-1 and mac_roman both decode every byte, so order decides the result.

        No single default can serve both, which is why the order is injectable: a
        genuinely latin-1 file is still readable by naming that order explicitly.
        """
        self._write_metadata("latin-1")

        from microbiome_analysis import load_metadata

        metadata = load_metadata(
            self.base / "Datesets" / "metadata.csv",
            encodings=("utf-8", "latin-1"),
        )
        self.assertEqual(list(metadata["run_accession"]), ["S1", "S2"])
        self.assertIn("Acu\xf1a", set(metadata["Author"]))

    def test_extract_tree_tips_from_nexus(self):
        self.write_text(
            "Datesets/tree.nexus",
            (
                "#NEXUS\n"
                "Begin trees;\n"
                "Tree tree1 = [&R] ('ASV1':0.1,('ASV2':0.2,'ASV3':0.3):0.4);\n"
                "End;\n"
            ),
        )

        from microbiome_analysis import extract_tree_tips

        tips = extract_tree_tips(self.base / "Datesets" / "tree.nexus")
        self.assertEqual(tips, {"ASV1", "ASV2", "ASV3"})

    def test_kruskal_wallis_tests(self):
        from microbiome_analysis import kruskal_wallis_tests

        alpha = pd.DataFrame({
            "Developmental_Stage": ["Seeds"] * 10 + ["Seedlings"] * 10 + ["soil"] * 10,
            "shannon": np.concatenate([
                np.random.default_rng(42).normal(3.0, 0.5, 10),
                np.random.default_rng(43).normal(4.0, 0.5, 10),
                np.random.default_rng(44).normal(2.5, 0.5, 10),
            ]),
            "observed_features": np.concatenate([
                np.random.default_rng(42).integers(50, 150, 10),
                np.random.default_rng(43).integers(100, 200, 10),
                np.random.default_rng(44).integers(30, 100, 10),
            ]),
        })
        result = kruskal_wallis_tests(alpha, "Developmental_Stage", metrics=["shannon", "observed_features"])
        self.assertEqual(len(result), 2)
        self.assertIn("H_statistic", result.columns)
        self.assertIn("p_value", result.columns)
        # With these well-separated groups, we expect significant results
        for _, row in result.iterrows():
            self.assertTrue(np.isfinite(row["H_statistic"]))
            self.assertGreater(row["H_statistic"], 0)
            self.assertGreaterEqual(row["p_value"], 0)
            self.assertLessEqual(row["p_value"], 1)

    def test_dunns_posthoc(self):
        from microbiome_analysis import dunns_posthoc

        alpha = pd.DataFrame({
            "group": ["A"] * 10 + ["B"] * 10 + ["C"] * 10,
            "shannon": np.concatenate([
                np.random.default_rng(1).normal(2.0, 0.3, 10),
                np.random.default_rng(2).normal(4.0, 0.3, 10),
                np.random.default_rng(3).normal(3.0, 0.3, 10),
            ]),
        })
        result = dunns_posthoc(alpha, "group", metric="shannon")
        self.assertEqual(len(result), 3)  # 3 pairwise comparisons
        self.assertIn("z_statistic", result.columns)
        self.assertIn("p_adjusted", result.columns)
        for _, row in result.iterrows():
            self.assertGreaterEqual(row["p_adjusted"], 0)
            self.assertLessEqual(row["p_adjusted"], 1)

    def test_differential_abundance(self):
        from microbiome_analysis import differential_abundance

        rng = np.random.default_rng(42)
        n = 20
        abundance = pd.DataFrame({
            "TaxonA": np.concatenate([rng.uniform(0.1, 0.5, n), rng.uniform(0.01, 0.05, n)]),
            "TaxonB": rng.uniform(0.01, 0.1, 2 * n),
            "TaxonC": np.concatenate([rng.uniform(0.0, 0.01, n), rng.uniform(0.2, 0.4, n)]),
        }, index=[f"S{i}" for i in range(2 * n)])
        abundance.index.name = "run_accession"

        metadata = pd.DataFrame({
            "run_accession": [f"S{i}" for i in range(2 * n)],
            "group": ["GroupA"] * n + ["GroupB"] * n,
        })

        result = differential_abundance(abundance, metadata, "group", "GroupA", "GroupB", min_prevalence=0.1)
        self.assertGreater(len(result), 0)
        self.assertIn("log2_fold_change", result.columns)
        self.assertIn("p_adjusted", result.columns)
        for _, row in result.iterrows():
            self.assertGreaterEqual(row["p_value"], 0)
            self.assertLessEqual(row["p_value"], 1)
            self.assertGreaterEqual(row["p_adjusted"], 0)

    def test_core_microbiome(self):
        from microbiome_analysis import core_microbiome

        abundance = pd.DataFrame({
            "CoreTaxon": [0.3, 0.25, 0.35, 0.28, 0.0, 0.0],
            "RareTaxon": [0.0, 0.01, 0.0, 0.0, 0.0, 0.02],
            "SharedTaxon": [0.1, 0.15, 0.12, 0.08, 0.1, 0.11],
        }, index=["S1", "S2", "S3", "S4", "S5", "S6"])
        abundance.index.name = "run_accession"

        metadata = pd.DataFrame({
            "run_accession": ["S1", "S2", "S3", "S4", "S5", "S6"],
            "group": ["A", "A", "A", "A", "B", "B"],
        })

        result = core_microbiome(abundance, metadata, "group", prevalence_threshold=0.5, abundance_threshold=0.001)
        self.assertGreater(len(result), 0)
        self.assertIn("taxon", result.columns)
        self.assertIn("prevalence", result.columns)
        # CoreTaxon should be core in group A (present in 4/4 samples above thresholds)
        core_a = result[(result["group"] == "A") & (result["taxon"] == "CoreTaxon")]
        self.assertEqual(len(core_a), 1)

    def test_indicator_species(self):
        from microbiome_analysis import indicator_species

        abundance = pd.DataFrame({
            "IndicatorA": [0.5, 0.4, 0.45, 0.0, 0.0, 0.0],
            "IndicatorB": [0.0, 0.0, 0.0, 0.3, 0.35, 0.4],
            "Generalist": [0.1, 0.12, 0.11, 0.09, 0.1, 0.08],
        }, index=["S1", "S2", "S3", "S4", "S5", "S6"])
        abundance.index.name = "run_accession"

        metadata = pd.DataFrame({
            "run_accession": ["S1", "S2", "S3", "S4", "S5", "S6"],
            "group": ["A", "A", "A", "B", "B", "B"],
        })

        result = indicator_species(abundance, metadata, "group", permutations=99, seed=42)
        self.assertEqual(len(result), 3)
        self.assertIn("indval", result.columns)
        # IndicatorA should be indicator for group A
        ind_a = result[result["taxon"] == "IndicatorA"]
        self.assertEqual(ind_a.iloc[0]["indicator_group"], "A")
        self.assertGreater(ind_a.iloc[0]["indval"], 0.5)

    def test_hedges_g(self):
        from microbiome_analysis import hedges_g

        rng = np.random.default_rng(42)
        group1 = rng.normal(5.0, 1.0, 30)
        group2 = rng.normal(3.0, 1.0, 30)
        g, se = hedges_g(group1, group2)
        self.assertTrue(np.isfinite(g))
        self.assertTrue(np.isfinite(se))
        self.assertGreater(g, 0)  # group1 mean > group2 mean
        self.assertGreater(se, 0)

    def test_hedges_g_small_sample(self):
        from microbiome_analysis import hedges_g

        g, se = hedges_g(np.array([1.0]), np.array([2.0, 3.0]))
        self.assertTrue(np.isnan(g))

    def test_random_effects_meta(self):
        from microbiome_analysis import random_effects_meta

        effects = np.array([0.5, 0.8, 0.3, 0.6, 0.4])
        se = np.array([0.2, 0.15, 0.25, 0.18, 0.22])
        result = random_effects_meta(effects, se)
        self.assertIn("pooled_effect", result)
        self.assertIn("I_squared", result)
        self.assertIn("tau_squared", result)
        self.assertEqual(result["k"], 5)
        self.assertTrue(np.isfinite(result["pooled_effect"]))
        self.assertGreaterEqual(result["I_squared"], 0)
        self.assertGreaterEqual(result["tau_squared"], 0)
        self.assertGreater(result["ci_upper"], result["ci_lower"])

    def test_nmds(self):
        from microbiome_analysis import nmds

        rng = np.random.default_rng(42)
        n = 10
        dist = rng.uniform(0, 1, (n, n))
        dist = (dist + dist.T) / 2
        np.fill_diagonal(dist, 0)

        result = nmds(dist, dimensions=2, max_iter=50, seed=0)
        self.assertEqual(result.shape, (n, 2))
        self.assertIn("NMDS1", result.columns)
        self.assertIn("NMDS2", result.columns)

    def test_per_study_effect_sizes(self):
        from microbiome_analysis import per_study_effect_sizes

        rng = np.random.default_rng(42)
        alpha = pd.DataFrame({
            "Developmental_Stage": (["Seeds"] * 10 + ["Seedlings"] * 10) * 3,
            "study_alias": ["Study1"] * 20 + ["Study2"] * 20 + ["Study3"] * 20,
            "shannon": np.concatenate([
                rng.normal(3.0, 0.5, 10), rng.normal(4.0, 0.5, 10),
                rng.normal(2.5, 0.3, 10), rng.normal(3.5, 0.3, 10),
                rng.normal(3.2, 0.4, 10), rng.normal(3.8, 0.4, 10),
            ]),
        })

        result = per_study_effect_sizes(
            alpha, "Developmental_Stage", "study_alias",
            metric="shannon", group1="Seeds", group2="Seedlings",
        )
        self.assertEqual(len(result), 3)
        self.assertIn("hedges_g", result.columns)
        self.assertIn("se", result.columns)
        for _, row in result.iterrows():
            self.assertTrue(np.isfinite(row["hedges_g"]))
            self.assertTrue(np.isfinite(row["se"]))

    def test_run_pipeline_writes_outputs_for_tiny_dataset(self):
        self.write_text(
            "Datesets/asv-counts.csv",
            (
                "id,ASV1,ASV2,ASV3,ASV4\n"
                "S1,10,0,5,0\n"
                "S2,0,20,1,0\n"
                "S3,3,1,0,7\n"
                "S4,8,2,3,1\n"
                "S5,0,15,0,5\n"
                "S6,5,0,8,2\n"
            ),
        )
        self.write_text(
            "Datesets/taxa.csv",
            (
                "OTUID,Kingdom,Phyla,Class,Order,Family,Genus,Species\n"
                "ASV1,Bacteria,Proteobacteria,C1,O1,F1,G1,sp1\n"
                "ASV2,Bacteria,Firmicutes,C2,O2,F2,G2,sp2\n"
                "ASV3,Unassigned,unassigned_phyla,unassigned_class,unassigned_order,unassigned_family,unassigned_genus,unassigned_species\n"
                "ASV4,Bacteria,Actinobacteriota,C3,O3,F3,G3,sp3\n"
            ),
        )
        self.write_text(
            "Datesets/metadata.csv",
            (
                "run_accession,,Author,Geography,Cultivar,Developmental_Stage,"
                "Scientific_Name,sample_description,Title,study_alias,command_to_download_the_reads\n"
                "S1,,A1,G1,C1,Seeds,PlantA,desc,title,Study1,cmd\n"
                "S2,,A1,G1,C1,Seedlings,PlantA,desc,title,Study1,cmd\n"
                "S3,,A2,G2,C2,Seeds,PlantB,desc,title,Study2,cmd\n"
                "S4,,A2,G2,C2,Seeds,PlantA,desc,title,Study1,cmd\n"
                "S5,,A2,G2,C2,Seedlings,PlantB,desc,title,Study2,cmd\n"
                "S6,,A1,G1,C1,Seeds,PlantA,desc,title,Study2,cmd\n"
            ),
        )
        self.write_text(
            "Datesets/tree.nexus",
            (
                "#NEXUS\n"
                "Begin trees;\n"
                "Tree tree1 = [&R] ('ASV1':0.1,('ASV2':0.2,('ASV3':0.3,'ASV4':0.1):0.2):0.2);\n"
                "End;\n"
            ),
        )
        config = {
            "inputs": {
                "counts": "Datesets/asv-counts.csv",
                "taxonomy": "Datesets/taxa.csv",
                "metadata": "Datesets/metadata.csv",
                "tree": "Datesets/tree.nexus",
            },
            "filters": {
                "min_library_size": 1,
                "min_feature_prevalence": 1,
                "min_feature_total_count": 1,
                "max_retained_features": 10,
                "beta_top_taxa": 10,
            },
            "analysis": {
                "primary_group": "Developmental_Stage",
                "secondary_groups": ["study_alias", "Scientific_Name", "Geography"],
                "permutations": 9,
                "max_samples": None,
                "beta_rank": "Genus",
                "pcoa_axes": 2,
            },
            "differential_abundance": {
                "groups": ["Seeds", "Seedlings"],
                "min_prevalence": 0.1,
            },
            "taxonomy_analysis": {
                "core_prevalence_threshold": 0.5,
                "core_abundance_threshold": 0.001,
                "indicator_permutations": 9,
            },
            "outputs": {
                "directory": "analysis_results",
                "generate_plots": False,
            },
        }
        config_path = self.base / "analysis_config.json"
        config_path.write_text(json.dumps(config), encoding="utf-8")

        from microbiome_analysis import run_pipeline

        run_pipeline(config_path)

        results = self.base / "analysis_results"

        # Original outputs
        self.assertTrue((results / "qc" / "sample_qc_summary.csv").exists())
        self.assertTrue((results / "diversity" / "alpha_diversity.csv").exists())
        self.assertTrue((results / "diversity" / "permanova_results.csv").exists())
        self.assertTrue((results / "taxonomy" / "genus_relative_abundance.csv").exists())

        # New diversity outputs
        self.assertTrue((results / "diversity" / "kruskal_wallis_results.csv").exists())
        self.assertTrue((results / "diversity" / "dunns_posthoc_results.csv").exists())
        self.assertTrue((results / "diversity" / "nmds_coordinates.csv").exists())

        # Differential abundance outputs
        self.assertTrue((results / "differential_abundance" / "da_results.csv").exists())

        # Taxonomy analysis outputs
        self.assertTrue((results / "taxonomy" / "core_microbiome.csv").exists())
        self.assertTrue((results / "taxonomy" / "indicator_species.csv").exists())

        # Meta-analysis outputs
        self.assertTrue((results / "meta_analysis" / "per_study_effect_sizes.csv").exists())

        # Verify alpha diversity row count
        with (results / "diversity" / "alpha_diversity.csv").open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 6)

        # Verify Kruskal-Wallis results
        kw = pd.read_csv(results / "diversity" / "kruskal_wallis_results.csv")
        self.assertGreater(len(kw), 0)
        self.assertIn("H_statistic", kw.columns)

        # Verify NMDS has correct shape
        nmds_df = pd.read_csv(results / "diversity" / "nmds_coordinates.csv")
        self.assertEqual(len(nmds_df), 6)
        self.assertIn("NMDS1", nmds_df.columns)

        # Verify DA results
        da = pd.read_csv(results / "differential_abundance" / "da_results.csv")
        if len(da) > 0:
            self.assertIn("log2_fold_change", da.columns)
            self.assertIn("p_adjusted", da.columns)

        # Verify core microbiome
        core = pd.read_csv(results / "taxonomy" / "core_microbiome.csv")
        if len(core) > 0:
            self.assertIn("taxon", core.columns)
            self.assertIn("prevalence", core.columns)

        # Verify indicator species
        indval = pd.read_csv(results / "taxonomy" / "indicator_species.csv")
        if len(indval) > 0:
            self.assertIn("indval", indval.columns)

        # R bridge outputs
        self.assertTrue((results / "qc" / "filtered_feature_ids.csv").exists())
        phylo_status = pd.read_csv(results / "diversity" / "phylogenetic_metrics_status.csv")
        self.assertIn("r_script_status", phylo_status.columns)

    def test_export_filtered_feature_ids(self):
        from microbiome_analysis import export_filtered_feature_ids

        out = self.base / "qc" / "filtered_feature_ids.csv"
        export_filtered_feature_ids(["ASV1", "ASV2", "ASV3"], out)
        self.assertTrue(out.exists())
        df = pd.read_csv(out)
        self.assertIn("OTUID", df.columns)
        self.assertEqual(list(df["OTUID"]), ["ASV1", "ASV2", "ASV3"])

    def test_run_phyloseq_analysis_missing_rscript(self):
        from unittest.mock import patch

        from microbiome_analysis import run_phyloseq_analysis

        with patch("microbiome_analysis.shutil.which", return_value=None):
            result = run_phyloseq_analysis(
                Path("fake_config.json"), Path("fake_features.csv"), self.base,
            )
        self.assertEqual(result["status"], "skipped")
        self.assertEqual(result["reason"], "rscript_not_found")

    def test_integrate_r_outputs_with_faith_pd(self):
        from microbiome_analysis import integrate_r_outputs

        div_dir = self.base / "diversity"
        div_dir.mkdir(parents=True, exist_ok=True)
        faith = pd.DataFrame({
            "sample_id": ["S1", "S2", "S3"],
            "faith_pd": [5.2, 8.1, 3.4],
            "phylogenetic_richness": [10, 15, 7],
        })
        faith.to_csv(div_dir / "faith_pd.csv", index=False)

        alpha = pd.DataFrame({
            "run_accession": ["S1", "S2", "S3"],
            "shannon": [2.1, 3.0, 1.5],
        })
        result = integrate_r_outputs(self.base, alpha)
        self.assertIn("faith_pd", result.columns)
        self.assertAlmostEqual(result.loc[0, "faith_pd"], 5.2)
        self.assertAlmostEqual(result.loc[1, "faith_pd"], 8.1)

    def test_integrate_r_outputs_no_faith_file(self):
        from microbiome_analysis import integrate_r_outputs

        alpha = pd.DataFrame({
            "run_accession": ["S1", "S2"],
            "shannon": [2.1, 3.0],
        })
        result = integrate_r_outputs(self.base, alpha)
        self.assertNotIn("faith_pd", result.columns)
        self.assertEqual(len(result), 2)

    def test_extract_timepoint_denmark_samples(self):
        from microbiome_analysis import extract_timepoint

        metadata = pd.DataFrame({
            "run_accession": ["D1", "D2", "D3", "D4"],
            "study_alias": ["PRJNA649549"] * 4,
            "sample_description": [
                "seed:T0.Seed:NA:something",
                "seed:11 DAS:stuff:here",
                "seed:17DAS:other:data",
                "soil:3DBS:bulk:sample",
            ],
        })
        result = extract_timepoint(metadata)
        self.assertEqual(list(result["Timepoint"]), ["T0", "11_DAS", "17_DAS", "3_DBS"])

    def test_extract_timepoint_no_denmark(self):
        from microbiome_analysis import extract_timepoint

        metadata = pd.DataFrame({
            "run_accession": ["S1", "S2"],
            "study_alias": ["PRJNA926007", "PRJNA926007"],
            "sample_description": ["some desc", "other desc"],
        })
        result = extract_timepoint(metadata)
        self.assertEqual(list(result["Timepoint"]), ["", ""])


if __name__ == "__main__":
    unittest.main()
