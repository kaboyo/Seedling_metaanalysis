#!/usr/bin/env python3
"""
Independent recomputation of every headline statistic in Manuscript.docx.

Deliberately NOT a re-run of phyloseq_analysis.qmd: the statistics are
reimplemented from the raw inputs (asv-counts.csv, taxa.csv, metadata.csv)
in numpy/scipy so that agreement is evidence rather than tautology.

Inputs come from analysis_config.json (`inputs.counts`, `inputs.taxonomy`,
`inputs.metadata`), resolved relative to the config file like the rest of the
repo:
  asv-counts.csv   2457 rows x 129071 cols
  taxa.csv         129070 x 8
  metadata.csv     2457 x 11 (latin-1!)

Usage:
  python3 validation/verify_manuscript.py --config analysis_config.json

Run order: STEP 1 is slow (~5 min, 638 MB text parse) and checkpoints to
<outputs.directory>/validation/asv_sparse.npz; later steps reload the
checkpoint, so a second run is fast.

Note: step6_alpha's depth models need statsmodels. If it is not installed that
one table is skipped and everything else still runs, matching how the Python
pipeline skips the R phylogenetic step when Rscript is absent.
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from scipy.spatial.distance import pdist, squareform
import argparse, csv, json, os
from pathlib import Path

try:
    import statsmodels.formula.api as smf
except ModuleNotFoundError:                         # optional, see step6_alpha
    smf = None

# Populated by configure() from analysis_config.json. Module-level so the step
# functions below keep referring to plain COUNTS / TAXA / META as they always did.
COUNTS = TAXA = META = None
WORK = None


def _resolve(base_dir: Path, configured_path: str) -> Path:
    """Same base-dir rule as microbiome_analysis.resolve_path and the Quarto notebook."""
    path = Path(configured_path)
    return path if path.is_absolute() else (base_dir / path).resolve()


def configure(config_path: Path) -> dict:
    """Point the recomputation at the inputs and output dir named by the config."""
    global COUNTS, TAXA, META, WORK
    config_path = Path(config_path).resolve()
    with config_path.open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    base_dir = config_path.parent

    inputs = config["inputs"]
    COUNTS = str(_resolve(base_dir, inputs["counts"]))
    TAXA = str(_resolve(base_dir, inputs["taxonomy"]))
    META = str(_resolve(base_dir, inputs["metadata"]))

    WORK = _resolve(base_dir, config["outputs"]["directory"]) / "validation"
    WORK.mkdir(parents=True, exist_ok=True)
    return config


def w(name: str) -> str:
    """Path to a checkpoint or result file inside the validation work directory."""
    if WORK is None:
        raise RuntimeError("configure() must run before any work path is used")
    return str(WORK / name)


# ─────────────────────────────────────────────────────────────────────────────
# STEP 1. Parse the count matrix into CSR.
#
# Why not pd.read_csv whole: 2456 x 129070 int32 dense is ~1.2 GB, and pandas
# needs ~2-3x that transiently on a 16 GB machine. Chunked read -> per-chunk
# csr_matrix -> vstack keeps peak RSS at ~800 MB. The matrix is 99.93% zero,
# so the sparse form is 0.6 MB on disk.
# ─────────────────────────────────────────────────────────────────────────────
def step1_parse():
    hdr = pd.read_csv(COUNTS, nrows=0)
    feature_ids = np.array(hdr.columns[1:], dtype=object)

    blocks, sample_ids = [], []
    for chunk in pd.read_csv(COUNTS, chunksize=150, dtype={0: str},
                             engine='c', low_memory=False):
        sample_ids.extend(chunk.iloc[:, 0].astype(str).tolist())
        arr = chunk.iloc[:, 1:].to_numpy(dtype=np.int32, copy=False)
        blocks.append(sp.csr_matrix(arr))
        del arr, chunk

    X = sp.vstack(blocks, format='csr')
    sample_ids = np.array(sample_ids, dtype=object)
    sp.save_npz(w('asv_sparse.npz'), X)
    np.save(w('feature_ids.npy'), feature_ids)
    np.save(w('sample_ids.npy'), sample_ids)
    return X, feature_ids, sample_ids
    # -> (2456, 129070), nnz=222811, total reads=34,993,355


def verify_parse(X, sample_ids, rows=(0, 700, 2455)):
    """Re-read selected rows with the stdlib csv reader and demand exact
    agreement on sample ID, row sum and nonzero count. Guards against a
    silent dtype/alignment error in the chunked path."""
    want = {}
    with open(COUNTS, newline='') as fh:
        rd = csv.reader(fh); next(rd)
        for i, row in enumerate(rd):
            if i in rows:
                want[i] = (row[0],
                           sum(int(v) for v in row[1:]),
                           sum(1 for v in row[1:] if v != '0'))
            if i > max(rows):
                break
    for i, (sid, tot, nz) in want.items():
        assert sid == sample_ids[i], (sid, sample_ids[i])
        assert tot == int(X[i].sum()), (i, tot, int(X[i].sum()))
        assert nz == X[i].nnz, (i, nz, X[i].nnz)


# ─────────────────────────────────────────────────────────────────────────────
# STEP 2. Load metadata / taxonomy and align to matrix order.
#
# metadata.csv is NOT utf-8: byte 0x96 inside the author string "Acuna"
# raises UnicodeDecodeError. latin-1 round-trips it.
# ─────────────────────────────────────────────────────────────────────────────
FAMILY_MAP = {  # transcribed from family_map tribble, phyloseq_analysis.qmd:1438
    'Brassica oleracea': 'Brassicaceae', 'Brassica napus': 'Brassicaceae',
    'Brassica rapa': 'Brassicaceae', 'Arabidopsis thaliana': 'Brassicaceae',
    'Diplotaxis tenuifolia': 'Brassicaceae', 'Raphanus sativus': 'Brassicaceae',
    'Daucus carota': 'Apiaceae', 'Petroselinum crispum': 'Apiaceae',
    'Glycine max': 'Fabaceae', 'Glycine clandestina': 'Fabaceae',
    'Phaseolus vulgaris': 'Fabaceae', 'Medicago sativa': 'Fabaceae',
    'Medicago truncatula': 'Fabaceae', 'Triticum aestivum': 'Poaceae',
    'Lactuca sativa': 'Asteraceae', 'Solanum lycopersicum': 'Solanaceae',
}


def step2_align(X, feature_ids, sample_ids):
    md = pd.read_csv(META, low_memory=False, encoding='latin-1')
    md = md.drop(columns=[c for c in md.columns if c.startswith('Unnamed')])
    tax = pd.read_csv(TAXA, low_memory=False, encoding='latin-1')

    # SRR18176465 (Glycine max seedling, PRJNA811248) is in metadata but has
    # no row in the count matrix. True starting N is 2456, not the 2457 the
    # manuscript reports as "samples retrieved".
    missing = sorted(set(md.run_accession) - set(sample_ids))
    assert missing == ['SRR18176465'], missing
    assert not (set(sample_ids) - set(md.run_accession))

    md = md.set_index('run_accession').loc[list(sample_ids)].reset_index()
    tax = tax.set_index('OTUID').loc[list(feature_ids)].reset_index()
    assert (md.run_accession.values == sample_ids).all()
    assert (tax.OTUID.values == feature_ids).all()

    md['Family'] = md.Scientific_Name.map(FAMILY_MAP).fillna('Other')
    md['depth_raw'] = np.asarray(X.sum(axis=1)).ravel()
    md['nobs_raw'] = np.diff(X.tocsr().indptr)
    return md, tax


# ─────────────────────────────────────────────────────────────────────────────
# STEP 3. Cross-study ASV comparability.
#
# The manuscript pools raw counts across 8 studies and claims comparisons are
# "made on common ASVs". That is testable: how many ASVs are actually common?
# ─────────────────────────────────────────────────────────────────────────────
def step3_comparability(X, md, studies):
    pres = {s: np.asarray((X[np.where(md.study_alias.values == s)[0]] > 0)
                          .sum(axis=0)).ravel() > 0
            for s in studies}
    P = pd.DataFrame(pres)
    nstud = P.sum(axis=1).values
    obs = nstud[nstud > 0]
    # -> 119135 features observed; 98.79% in exactly ONE study; max 4 of 8

    rows = []
    for s in studies:
        f_s = P[s].values
        rows.append(dict(
            study=s,
            n_samples=int((md.study_alias == s).sum()),
            n_features=int(f_s.sum()),
            n_private=int((f_s & (nstud == 1)).sum()),
            pct_private=100 * (f_s & (nstud == 1)).sum() / f_s.sum(),
            reads=int(X[np.where(md.study_alias.values == s)[0]].sum())))
    per_study = pd.DataFrame(rows)

    n = len(studies)
    J = np.zeros((n, n))
    for i, a in enumerate(studies):
        for j, b in enumerate(studies):
            fa, fb = P[a].values, P[b].values
            uni = (fa | fb).sum()
            J[i, j] = (fa & fb).sum() / uni if uni else 0.0
    jac = pd.DataFrame(J, index=studies, columns=studies)
    # -> median off-diagonal Jaccard 0.0000, max 0.0142

    return per_study, jac, P, nstud


def agg(rank_vals, Xm):
    """Sum columns of Xm within levels of rank_vals. Blank/NA -> 'Unassigned',
    matching aggregate_by_rank() at phyloseq_analysis.qmd:1365."""
    rv = pd.Series(rank_vals).fillna('Unassigned').replace('', 'Unassigned').values
    cats, inv = np.unique(rv, return_inverse=True)
    M = sp.csr_matrix((np.ones(len(inv)), (np.arange(len(inv)), inv)),
                      shape=(len(inv), len(cats)))
    return cats, (Xm @ M).tocsr()


# ─────────────────────────────────────────────────────────────────────────────
# STEP 4. Reproduce the filtering cascade from analysis_config.json.
#
# Order matters and is NOT the order the manuscript describes. The pipeline
# (qmd:1050-1067) removes soil BEFORE the library-size floor, so the 1886
# figure is already plant-only -- there is no later "1886 -> 1883 after
# removing 48 soil" step for the Results text to describe.
# ─────────────────────────────────────────────────────────────────────────────
def step4_cascade(X, md, tax):
    acc = []
    fmask = np.ones(X.shape[1], bool)
    smask = np.ones(X.shape[0], bool)

    def rec(step, note=''):
        acc.append(dict(step=step, n_samples=int(smask.sum()),
                        n_features=int(fmask.sum()),
                        total_reads=int(X[np.ix_(smask, fmask)].sum()), note=note))

    rec('0. raw matrix as parsed',
        'metadata has 2457 rows; SRR18176465 absent from counts')

    K, O = tax.Kingdom.fillna('').values, tax.Order.fillna('').values
    FAM, PH = tax.Family.fillna('').values, tax.Phyla.fillna('').values

    for newf, label, note in [
        (K == 'Bacteria', '1. keep Kingdom==Bacteria',
         f'removed {(K != "Bacteria").sum()} incl. {(K == "Unassigned").sum()} Unassigned'),
        (O != 'Chloroplast', '2. drop Order==Chloroplast', ''),
        (FAM != 'Mitochondria', '3. drop Family==Mitochondria', ''),
        (~np.isin(PH, ['Cyanobacteria', 'Cyanobacteriota']),
         '4. drop Cyanobacteria phyla', ''),
    ]:
        fmask &= newf
        rec(label, note)

    soil = (md.Developmental_Stage == 'soil').values
    smask &= ~soil
    rec('5. remove soil samples', f'removed {soil.sum()} soil')

    depth_f = np.asarray(X[:, fmask].sum(axis=1)).ravel()
    smask &= (depth_f >= 500)                       # min_library_size
    rec('6. library size >= 500', 'applied to Bacteria-filtered depths')

    def feat_step(pred_on_sub, label, note=''):
        nonlocal fmask
        idx = np.where(fmask)[0]
        full = np.zeros(X.shape[1], bool)
        full[idx[pred_on_sub]] = True
        fmask = full
        rec(label, note)

    S = X[np.ix_(smask, fmask)]
    feat_step(np.asarray(S.sum(axis=0)).ravel() > 0, '7. drop zero-sum features')
    S = X[np.ix_(smask, fmask)]
    feat_step(np.asarray((S > 0).sum(axis=0)).ravel() >= 3,
              '8. prevalence >= 3 samples')       # min_feature_prevalence
    S = X[np.ix_(smask, fmask)]
    feat_step(np.asarray(S.sum(axis=0)).ravel() >= 10,
              '9. feature total count >= 10')     # min_feature_total_count
    S = X[np.ix_(smask, fmask)]
    ts = np.asarray(S.sum(axis=0)).ravel()
    keep = np.zeros(len(ts), bool)
    keep[np.argsort(-ts)[:2000]] = True           # max_retained_features
    feat_step(keep, '10. cap at 2000 most abundant features')

    casc = pd.DataFrame(acc)
    casc['pct_reads_retained'] = (100 * casc.total_reads / casc.total_reads.iloc[0]).round(2)
    return casc, smask, fmask
    # -> 1886 samples x 2000 features, 52.14% of reads. 1886 reproduces exactly.


# ─────────────────────────────────────────────────────────────────────────────
# STEP 5. Composition (manuscript Tables 2 and 3).
# ─────────────────────────────────────────────────────────────────────────────
def step5_composition(Xf, taxf, n):
    pc, Ph = agg(taxf.Phyla.values, Xf)
    rs = np.asarray(Ph.sum(axis=1)).ravel()
    Prel = Ph.multiply(sp.csr_matrix(1 / np.where(rs == 0, 1, rs)).T).tocsr()
    Pd = pd.DataFrame(Prel.toarray(), columns=pc)
    phy = (pd.DataFrame({'phylum': pc,
                         'mean_rel': Pd.mean().values,
                         'median_rel': Pd.median().values,
                         'sd_rel': Pd.std(ddof=1).values,
                         'prevalence_n': (Pd > 0).sum().values})
           .assign(prevalence_pct=lambda d: 100 * d.prevalence_n / n)
           .sort_values('mean_rel', ascending=False).reset_index(drop=True))
    # -> Proteobacteria 0.7187/0.9257/1770; Firmicutes 0.1279; Actinobacteriota
    #    0.1255 (prevalence 940, manuscript says 937); Bacteroidota 0.0214

    gc_, Gt = agg(taxf.Genus.values, Xf)
    gtot = np.asarray(Gt.sum(axis=0)).ravel()
    gen = (pd.DataFrame({'genus': gc_, 'total_abundance': gtot,
                         'mean_abundance': gtot / n,
                         'prevalence_n': np.asarray((Gt > 0).sum(axis=0)).ravel()})
           .sort_values('total_abundance', ascending=False).reset_index(drop=True))
    # -> Table 3 reproduces exactly, all ten rows
    return phy, gen, Gt


# ─────────────────────────────────────────────────────────────────────────────
# STEP 6. Alpha diversity (manuscript Table S1 / Results para 3).
#
# The 3 samples whose reads all fall outside the 2000-feature cap have
# undefined Shannon. Imputing 0 for them reproduces the manuscript's
# H = 38.2 at n = 1886; dropping them gives H = 36.4 at n = 1883.
# ─────────────────────────────────────────────────────────────────────────────
def eps2(H, N, k):
    """Epsilon-squared effect size for Kruskal-Wallis."""
    return (H - k + 1) / (N - k)


def step6_alpha(Xd, mdf):
    obs = (Xd > 0).sum(axis=1).astype(float)
    dep = Xd.sum(axis=1)
    p_ = np.divide(Xd, np.where(dep[:, None] == 0, 1, dep[:, None]))
    with np.errstate(divide='ignore', invalid='ignore'):
        lp = np.where(p_ > 0, np.log(p_), 0.0)
    shan = -(p_ * lp).sum(axis=1)
    simp = 1 - (p_ ** 2).sum(axis=1)
    shan[dep == 0] = np.nan
    simp[dep == 0] = np.nan
    even = np.where(obs > 1, shan / np.log(np.where(obs > 1, obs, np.e)), np.nan)

    A = mdf.copy()
    A['Observed'], A['Shannon'], A['Simpson'] = obs, shan, simp
    A['Evenness'], A['depth2000'] = even, dep

    rows = []
    for var in ['Developmental_Stage', 'Scientific_Name', 'Cultivar',
                'study_alias', 'Author', 'Family']:
        for metric in ['Observed', 'Shannon']:
            d = A[A[var].notna() & (A[var].astype(str).str.strip() != '')
                  & A[metric].notna()]
            vc = d[var].value_counts()
            keep = vc[vc >= 3].index                # min_group_n
            for lab, dd in [('all_groups', d), ('min_n3', d[d[var].isin(keep)])]:
                g = [x[metric].values for _, x in dd.groupby(var) if len(x) > 0]
                if len(g) < 2:
                    continue
                H, p = stats.kruskal(*g)
                rows.append(dict(variable=var, metric=metric, subset=lab,
                                 H=round(H, 3), p=p, n_samples=len(dd),
                                 n_groups=len(g),
                                 epsilon_sq=round(eps2(H, len(dd), len(g)), 4)))
    alpha = pd.DataFrame(rows)
    alpha['q_BH'] = stats.false_discovery_control(
        np.clip(alpha.p.values, 1e-300, 1.0), method='bh')

    # Reconciliation of the manuscript's H = 38.2 / n = 1886
    A2 = A.assign(Shannon0=np.where(np.isfinite(shan), shan, 0.0))
    H0, p0 = stats.kruskal(*[x.Shannon0.values
                             for _, x in A2.groupby('Developmental_Stage')])
    # -> H0 = 38.234, p = 6.28e-10 : exactly the reported value

    # Depth control. min_library_size is only 500 and Observed correlates with
    # depth (rho = 0.298), so the stage effect needs a covariate check.
    A3 = A[A.Shannon.notna()].copy()
    A3['logdepth'] = np.log10(A3.depth2000.clip(lower=1))
    A3['stage'] = (A3.Developmental_Stage == 'Seedlings').astype(int)
    depth_rows = []
    for metric in ['Shannon', 'Observed']:
        for nm, formula in [('stage only', f'{metric} ~ stage'),
                            ('+ logdepth', f'{metric} ~ stage + logdepth'),
                            ('+ logdepth + study',
                             f'{metric} ~ stage + logdepth + C(study_alias)')]:
            if smf is None:
                continue
            m = smf.ols(formula, data=A3).fit()
            depth_rows.append(dict(metric=metric, model=nm,
                                   stage_beta=m.params['stage'],
                                   se=m.bse['stage'], p=m.pvalues['stage'],
                                   adj_r2=m.rsquared_adj))
    return A, alpha, (H0, p0), pd.DataFrame(depth_rows)


# ─────────────────────────────────────────────────────────────────────────────
# STEP 7. Beta diversity: PERMANOVA and PERMDISP from scratch.
#
# Reimplemented rather than imported so the basis of every number is explicit.
# Validated against vegan by matching all five manuscript Table 4 rows to
# 4 significant figures.
# ─────────────────────────────────────────────────────────────────────────────
def bray(M):
    """Bray-Curtis on row-normalised (relative-abundance) profiles."""
    R = M / np.where(M.sum(1, keepdims=True) == 0, 1, M.sum(1, keepdims=True))
    return squareform(pdist(R, metric='braycurtis'))


def braw(M):
    """Bray-Curtis on raw counts, no normalisation. This is what the
    manuscript's Table 4 actually used."""
    return squareform(pdist(M, metric='braycurtis'))


def gower(D):
    """Gower-centred matrix G = -0.5 * D^2, double-centred."""
    A_ = -0.5 * D ** 2
    return A_ - A_.mean(0, keepdims=True) - A_.mean(1, keepdims=True) + A_.mean()


def hat(dm):
    Q, _ = np.linalg.qr(dm)
    return Q @ Q.T


def dmat(df, terms):
    """Design matrix with intercept and treatment-coded dummies, plus a map
    from term name to its column span (needed for marginal SS)."""
    cols = [np.ones((len(df), 1))]
    idx = {}
    for t in terms:
        d = pd.get_dummies(df[t].astype(str), drop_first=True).values.astype(float)
        start = sum(c.shape[1] for c in cols)
        idx[t] = (start, start + d.shape[1])
        cols.append(d)
    return np.hstack(cols), idx


def permanova(D, df, terms, marginal=False, nperm=999, seed=42):
    """Distance-based PERMANOVA (McArdle-Anderson). marginal=True gives
    type-III SS in the presence of the other terms, i.e. adonis2(by='margin')
    -- the model the manuscript needs in order to separate study identity
    from biology, and does not report."""
    rng = np.random.default_rng(seed)
    G = gower(D)
    n = len(D)
    Xd_, idx = dmat(df, terms)
    Hf = hat(Xd_)
    SSt, SSf = np.trace(G), np.sum(G * Hf)
    dff = np.linalg.matrix_rank(Xd_) - 1
    dfr = n - 1 - dff

    out = []
    for t in terms:
        if marginal:
            keep = [i for i in range(Xd_.shape[1]) if not (idx[t][0] <= i < idx[t][1])]
            Hr = hat(Xd_[:, keep])
            ss = SSf - np.sum(G * Hr)
            dt = np.linalg.matrix_rank(Xd_) - np.linalg.matrix_rank(Xd_[:, keep])
            F = (ss / dt) / ((SSt - SSf) / dfr)
        else:
            Xu, _ = dmat(df, [t])
            Hu = hat(Xu)
            ss = np.sum(G * Hu)
            dt = np.linalg.matrix_rank(Xu) - 1
            F = (ss / dt) / ((SSt - ss) / (n - 1 - dt))

        cnt = 0
        for _ in range(nperm):
            p = rng.permutation(n)
            Gp = G[np.ix_(p, p)]
            if marginal:
                ssp = np.sum(Gp * Hf) - np.sum(Gp * Hr)
                Fp = (ssp / dt) / ((np.trace(Gp) - np.sum(Gp * Hf)) / dfr)
            else:
                ssp = np.sum(Gp * Hu)
                Fp = (ssp / dt) / ((np.trace(Gp) - ssp) / (n - 1 - dt))
            cnt += (Fp >= F - 1e-12)

        out.append(dict(term=t, df=dt, SS=ss, R2=ss / SSt, F=F,
                        p=(cnt + 1) / (nperm + 1),
                        model='marginal' if marginal else 'univariate', n=n))
    return pd.DataFrame(out)


def spatial_median(Y, tol=1e-10, itmax=500):
    """Weiszfeld iteration. vegan::betadisper defaults to type='median', and
    using the centroid instead changes stage F from 21.49 to 27.51 -- enough
    to look like a failed replication when it is only a different estimator."""
    m = Y.mean(0)
    for _ in range(itmax):
        d = np.sqrt(((Y - m) ** 2).sum(1))
        d[d < 1e-12] = 1e-12
        w = 1 / d
        mn = (Y * w[:, None]).sum(0) / w.sum()
        if np.abs(mn - m).max() < tol:
            return mn
        m = mn
    return m


def permdisp(D, groups, nperm=999, seed=42, kind='median'):
    """PERMDISP / betadisper. Bray-Curtis is non-Euclidean, so the Gower
    matrix has negative eigenvalues; vegan corrects distance-to-centre with
    sqrt(pos - neg), reproduced here."""
    rng = np.random.default_rng(seed)
    G = gower(D)
    w, v = np.linalg.eigh(G)
    pos, neg = w > 1e-9, w < -1e-9
    Yp = v[:, pos] * np.sqrt(w[pos])
    Yn = v[:, neg] * np.sqrt(-w[neg]) if neg.any() else None
    g = pd.Series(groups).astype(str).values
    lv = np.unique(g)
    ctr = spatial_median if kind == 'median' else (lambda Y: Y.mean(0))

    def stat(gg):
        d = np.empty(len(Yp))
        for l in lv:
            m = gg == l
            if m.sum() < 2:
                d[m] = 0
                continue
            dp = ((Yp[m] - ctr(Yp[m])) ** 2).sum(1)
            if Yn is not None:
                dn = ((Yn[m] - ctr(Yn[m])) ** 2).sum(1)
                d[m] = np.sqrt(np.clip(dp - dn, 0, None))
            else:
                d[m] = np.sqrt(dp)
        gr = [d[gg == l] for l in lv if (gg == l).sum() > 1]
        return stats.f_oneway(*gr).statistic, d

    F, d = stat(g)
    cnt = sum(stat(rng.permutation(g))[0] >= F - 1e-12 for _ in range(nperm))
    return F, (cnt + 1) / (nperm + 1), d


def step7_identify_table4_basis(Xd, Gden, mdf, ka, kg):
    """Table 4 is captioned 'genus-level Bray-Curtis' with n = 1883, but no
    genus-level variant reproduces its F. Enumerate the candidates."""
    variants = {
        'genus top100, rel':   (bray(Gden[kg]),  mdf.Developmental_Stage.values[kg]),
        'genus top100, raw':   (braw(Gden[kg]),  mdf.Developmental_Stage.values[kg]),
        'ASV top100, rel':     (bray(Xd[ka]),    mdf.Developmental_Stage.values[ka]),
        'ASV top100, raw':     (braw(Xd[ka]),    mdf.Developmental_Stage.values[ka]),
    }
    out = []
    for nm, (D, g) in variants.items():
        G = gower(D)
        Xu, _ = dmat(pd.DataFrame({'g': g}), ['g'])
        Hu = hat(Xu)
        ss = np.sum(G * Hu)
        dt = np.linalg.matrix_rank(Xu) - 1
        F = (ss / dt) / ((np.trace(G) - ss) / (len(D) - 1 - dt))
        out.append(dict(variant=nm, F=round(F, 3), R2=round(ss / np.trace(G), 4),
                        n=len(D)))
    return pd.DataFrame(out)
    # -> 'ASV top100, raw' gives F = 19.260, R2 = 0.0123, n = 1553.
    #    Manuscript reports F = 19.26, R2 = 0.012, n = 1883.
    #    All five Table 4 rows and all five PERMDISP F values match this
    #    basis. The caption and the n column are both wrong: the statistics
    #    are ASV-level, raw-count, and computed on 1553 samples.


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Independently recompute the manuscript's headline statistics.")
    parser.add_argument("--config", default="analysis_config.json",
                        help="Path to the analysis config JSON file.")
    args = parser.parse_args()
    configure(Path(args.config))
    if smf is None:
        print("statsmodels not installed - alpha depth models will be skipped")
    run()


def run() -> None:
    # All three checkpoint members must be present: the matrix alone is not enough
    # to resume, and a partial checkpoint used to fail on the np.load below.
    checkpoint = [w('asv_sparse.npz'), w('feature_ids.npy'), w('sample_ids.npy')]
    if all(os.path.exists(p) for p in checkpoint):
        X = sp.load_npz(w('asv_sparse.npz'))
        feature_ids = np.load(w('feature_ids.npy'), allow_pickle=True)
        sample_ids = np.load(w('sample_ids.npy'), allow_pickle=True)
    else:
        missing = [os.path.basename(p) for p in checkpoint if not os.path.exists(p)]
        if missing != [os.path.basename(p) for p in checkpoint]:
            print(f"incomplete checkpoint (missing {', '.join(missing)}) - reparsing")
        X, feature_ids, sample_ids = step1_parse()
    verify_parse(X, sample_ids)

    md, tax = step2_align(X, feature_ids, sample_ids)
    studies = sorted(md.study_alias.unique())

    per_study, jac, P, nstud = step3_comparability(X, md, studies)
    casc, smask, fmask = step4_cascade(X, md, tax)

    mdf = md[smask].reset_index(drop=True)
    Xf = X[np.ix_(smask, fmask)].tocsr()
    taxf = tax[fmask].reset_index(drop=True)
    Xd = Xf.toarray().astype(float)

    phy, gen, Gt = step5_composition(Xf, taxf, Xf.shape[0])
    A, alpha, (H0, p0), depth_models = step6_alpha(Xd, mdf)

    # top-100 subsets, as prepare_beta_phyloseq_rank(top_n=beta_top_taxa=100)
    asv100 = np.argsort(-np.asarray(Xf.sum(axis=0)).ravel())[:100]
    _, Gm = agg(taxf.Genus.values, Xf)
    Gden = Gm.toarray()
    gtop = np.argsort(-Gden.sum(0))[:100]
    Xa, Xg = Xd[:, asv100], Gden[:, gtop]
    ka, kg = Xa.sum(1) > 0, Xg.sum(1) > 0

    basis = step7_identify_table4_basis(Xa, Xg, mdf, ka, kg)

    terms = ['Developmental_Stage', 'Scientific_Name', 'Cultivar',
             'study_alias', 'Author']
    Draw_a = braw(Xa[ka])
    ma = mdf[ka].reset_index(drop=True)
    uni = permanova(Draw_a, ma, terms, marginal=False)
    marg = permanova(Draw_a, ma, terms, marginal=True)

    for df_, name in [(per_study, 'asv_sharing_by_study'), (jac, 'asv_jaccard_by_study'),
                      (casc, 'filtering_accounting'), (phy, 'phylum_composition_recomputed'),
                      (gen, 'genus_composition_recomputed'), (alpha, 'alpha_stats_recomputed'),
                      (depth_models, 'alpha_depth_models'), (basis, 'table4_basis_identification'),
                      (uni, 'permanova_univariate'), (marg, 'permanova_marginal')]:
        df_.to_csv(w(f'{name}.csv'), index=(name == 'asv_jaccard_by_study'))


if __name__ == '__main__':
    main()
