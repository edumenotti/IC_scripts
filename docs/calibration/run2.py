"""Throwaway: calibration checks that do not reuse the model's own kinship.

  orf_sim     null phenotypes simulated from ORF presence/absence relatedness, tested with the SNP-kinship LMM
  perm_kN     real phenotypes permuted within N relatedness clusters (keeps lineage means and real distributions)
  real        real phenotypes
Each with the raw phenotype and with a rank-based inverse normal transform (INT).
"""
import sys

import numpy as np
import pandas as pd
from scipy import optimize, stats as sps
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

from mutation_hunter import stats
from mutation_hunter.report import lambda_gc

z = np.load(sys.argv[1], allow_pickle=True)  # our own file
orf_path, out = sys.argv[2], sys.argv[3]
strains, gene, freq = z["strains"], z["gene"], z["freq"]
from mutation_hunter import data  # noqa: E402

Y = data.load_phenotypes().loc[strains]
Dsnp = data.load_distance().loc[strains, strains]
K = stats.kinship(Dsnp).to_numpy()
s, U = np.linalg.eigh(K)
s = np.clip(s, 0, None)
Korf = stats.kinship(pd.read_csv(orf_path, sep="\t", index_col=0).loc[strains, strains]).to_numpy()
so, Uo = np.linalg.eigh(Korf)
so = np.clip(so, 0, None)
n = len(strains)
ONE = np.ones((n, 1))
rng = np.random.default_rng(7)


def prepare(M):
    present = ~np.isnan(M)
    car, npres = np.nansum(M, 0), present.sum(0)
    M = M[:, (car >= 3) & (car < npres)]
    mean = np.nanmean(M, 0)
    idx = np.where(np.isnan(M))
    M[idx] = mean[idx[1]]
    return M


m = freq <= 0.02
Gv = z["G"][:, m].astype(float)
mats = {"variant": prepare(Gv.copy()), "gene": prepare(pd.DataFrame(Gv.T).groupby(gene[m]).max().T.to_numpy())}
del Gv
rot = {k: U.T @ M for k, M in mats.items()}


def lmm_p(UtM, y):
    uy, ux = U.T @ y, U.T @ ONE

    def nll(ld):
        v = s + np.exp(ld)
        b = (ux[:, 0] / v) @ uy / ((ux[:, 0] / v) @ ux[:, 0])
        return 0.5 * (n * np.log(((uy - ux[:, 0] * b) ** 2 / v).sum()) + np.log(v).sum())

    w = 1 / np.sqrt(s + np.exp(optimize.minimize_scalar(nll, bounds=(-10, 10), method="bounded").x))
    return stats.assoc(UtM * w[:, None], uy * w, ux * w[:, None])[2]


def inorm(y):
    return sps.norm.ppf((sps.rankdata(y) - 0.5) / len(y))


def orf_null(h2):
    u = Uo @ (np.sqrt(so) * rng.normal(size=n))
    return np.sqrt(h2) * u / u.std() + np.sqrt(1 - h2) * rng.normal(size=n)


def permuted(y, labels):
    y = y.copy()
    for c in np.unique(labels):
        i = np.flatnonzero(labels == c)
        y[i] = y[rng.permutation(i)]
    return y


link = hierarchy.linkage(squareform(Dsnp.to_numpy(), checks=False), "average")
sets = {"orf_sim": [orf_null(h2) for h2 in [0.3] * 10 + [0.7] * 10],
        "real": [Y[c].to_numpy(float) for c in Y]}
for k in (10, 30):
    lab = hierarchy.fcluster(link, k, "maxclust")
    sets[f"perm_k{k}"] = [permuted(Y[c].to_numpy(float), lab) for c in Y]

rows = []
for level, UtM in rot.items():
    for name, ys in sets.items():
        for tr, f in [("raw", lambda y: y), ("INT", inorm)]:
            ps = [lmm_p(UtM, f(y)) for y in ys]
            lam = [lambda_gc(p) for p in ps]
            ps_ols = [stats.assoc(mats[level], f(y), ONE)[2] for y in ys[:5]] if name == "orf_sim" and tr == "raw" else []
            row = {"level": level, "set": name, "transform": tr, "tests": UtM.shape[1],
                   "lambda_median": np.median(lam), "lambda_max": np.max(lam),
                   "tail_1e-3": np.mean([np.mean(p < 1e-3) for p in ps]) / 1e-3,
                   "ols_lambda_median": np.median([lambda_gc(p) for p in ps_ols]) if ps_ols else np.nan}
            rows.append(row)
            print(row, flush=True)
            pd.DataFrame(rows).to_csv(out, sep="\t", index=False, float_format="%.3g")
