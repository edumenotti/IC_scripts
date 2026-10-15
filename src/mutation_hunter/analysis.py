"""PheWAS (one variant or gene vs. all phenotypes) and GWAS (one phenotype vs. all genes)."""

import numpy as np
import pandas as pd

from . import data, stats, variants

COLUMNS = ["phenotype", "gene", "test", "n", "n_carriers", "beta", "se", "p", "p_hypergeom", "q"]
ORF_PCS = 10  # accessory-genome structure that SNP kinship misses (calibration grid)


def _structure() -> tuple[pd.DataFrame, pd.DataFrame]:
    """SNP kinship (random effect) and top ORF presence/absence coordinates (fixed covariates)."""
    return stats.kinship(data.load_distance()), stats.pcs(stats.kinship(data.load_orf_distance()), ORF_PCS)


def _null(y: pd.Series, K: pd.DataFrame, C: pd.DataFrame):
    """Analysed strains, raw and rank-normalised phenotype, design X and the LMM whitening matrix."""
    strains = y.dropna().index.intersection(K.index).intersection(C.index)
    raw = y[strains].to_numpy(float)
    yv = stats.inverse_normal(raw)
    X = np.column_stack([np.ones(len(strains)), C.loc[strains].to_numpy()])
    return strains, raw, yv, X, stats.lmm_null(K.loc[strains, strains].to_numpy(), yv, X)


def _test(G: pd.DataFrame, null, min_carriers: int) -> pd.DataFrame:
    """Test each column of G (1 = carrier) against the phenotype in null."""
    strains, raw, yv, X, W = null
    g = G.reindex(strains)  # strains lacking the gene -> NaN: missing, not reference
    present = g.notna().to_numpy()
    n, n_car = present.sum(axis=0), g.sum().to_numpy()
    gi = g.fillna(g.mean()).to_numpy(float)  # mean imputation: missing strains add no signal
    beta, se, p = stats.assoc(W @ gi, W @ yv, W @ X)
    ok = (n_car >= min_carriers) & (n_car < n)
    p_hg = [stats.hypergeom_p(gi[present[:, i], i] > 0, raw[present[:, i]]) if ok[i] else np.nan
            for i in range(g.shape[1])]
    return pd.DataFrame({
        "test": G.columns, "n": n, "n_carriers": n_car.astype(int),
        "beta": beta, "se": se, "p": np.where(ok, p, np.nan), "p_hypergeom": p_hg,
    })


def phewas(query: str, max_freq=0.02, min_carriers=3) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Test a variant ('YNL079C:A123V') or a gene burden ('YNL079C') against every phenotype.
    Returns (results, the rare variants of the gene)."""
    gene, var = variants.parse_query(query)
    seqs = data.load_gene(gene)
    rv = variants.rare_variants(gene, seqs, max_freq, min_carriers)
    if var:
        carriers = variants.carriers_of(gene, seqs, *var)
    else:
        carriers = sorted({s for c in rv["carriers"] for s in c})  # burden: any rare allele (CAST)
        if not carriers:
            raise ValueError(f"{gene} has no rare variant at max_freq={max_freq}, min_carriers={min_carriers}")
    G = variants.carrier_matrix(seqs, {query: carriers})
    Y, (K, C) = data.load_phenotypes(), _structure()
    res = pd.concat([_test(G, _null(Y[ph], K, C), min_carriers).assign(phenotype=ph) for ph in Y], ignore_index=True)
    res["gene"], res["q"] = gene, stats.bh(res["p"])
    return res[COLUMNS], rv


def gwas(phenotype: pd.Series, level="gene", max_freq=0.02, min_carriers=3) -> pd.DataFrame:
    """Test every gene burden (level='gene') or every rare variant (level='variant') against one phenotype."""
    null, parts = _null(phenotype, *_structure()), []
    for gene in data.genes():
        seqs = data.load_gene(gene)
        rv = variants.rare_variants(gene, seqs, max_freq, min_carriers)
        if rv.empty:
            continue
        G = variants.carrier_matrix(seqs, dict(zip(rv["variant"], rv["carriers"])))
        if level == "gene":
            G = G.max(axis=1).to_frame(gene)
        parts.append(_test(G, null, min_carriers).assign(gene=gene))
    res = pd.concat(parts, ignore_index=True).assign(phenotype=phenotype.name)
    res["q"] = stats.bh(res["p"])
    return res[COLUMNS]
