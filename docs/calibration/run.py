"""Throwaway calibration grid: rarity threshold x level x model.

Metrics per configuration:
  real_lambda  median / max lambda_GC over the 35 real conditions (part may be true polygenic signal)
  null_lambda  median lambda_GC over 20 simulated phenotypes with kinship structure and NO causal variant
  null_tail    observed / expected rate of p < 1e-3 under the null (1 = calibrated tail)
  power_b1/b2  fraction of 100 planted effects (beta = 1 or 2 SD) detected at Bonferroni 0.05
"""
import sys

import numpy as np
import pandas as pd
from scipy import optimize

from mutation_hunter import data, stats
from mutation_hunter.report import lambda_gc

z = np.load(sys.argv[1], allow_pickle=True)  # our own file
out = sys.argv[2]
strains, gene, freq = z["strains"], z["gene"], z["freq"]
G_all = z["G"].astype(float)
Y = data.load_phenotypes().loc[strains]
K = stats.kinship(data.load_distance().loc[strains, strains]).to_numpy()
s, U = np.linalg.eigh(K)
s = np.clip(s, 0, None)
top = np.argsort(s)[::-1][:10]
PC = U[:, top] * np.sqrt(s[top])
n = len(strains)
ONE = np.ones((n, 1))
MODELS = {"ols": (False, 0), "pcs5": (False, 5), "pcs10": (False, 10), "lmm": (True, 0), "lmm+pcs5": (True, 5)}


def prepare(M):
    """Keep columns with >= 3 carriers and >= 1 non-carrier among analysed strains; mean-impute absences."""
    present = ~np.isnan(M)
    car, npres = np.nansum(M, 0), present.sum(0)
    keep = (car >= 3) & (car < npres)
    M = M[:, keep]
    mean = np.nanmean(M, 0)
    idx = np.where(np.isnan(M))
    M[idx] = mean[idx[1]]
    return M


def matrices(F):
    m = freq <= F
    Gv = G_all[:, m]
    Gb = pd.DataFrame(Gv.T).groupby(gene[m]).max().T.to_numpy()  # all-NaN (gene absent) stays NaN
    return {"variant": prepare(Gv.copy()), "gene": prepare(Gb)}


def delta(uy, uX):
    def nll(ld):
        v = s + np.exp(ld)
        Xw = uX / v[:, None]
        b = np.linalg.solve(uX.T @ Xw, Xw.T @ uy)
        r = uy - uX @ b
        return 0.5 * (n * np.log((r * r / v).sum()) + np.log(v).sum())
    return np.exp(optimize.minimize_scalar(nll, bounds=(-10, 10), method="bounded").x)


def pvals(M, UtM, y, model):
    lmm, k = MODELS[model]
    X = np.column_stack([ONE, PC[:, :k]])
    if not lmm:
        return stats.assoc(M, y, X)[2]
    uy, uX = U.T @ y, U.T @ X
    w = 1 / np.sqrt(s + delta(uy, uX))
    return stats.assoc(UtM * w[:, None], uy * w, uX * w[:, None])[2]


rng = np.random.default_rng(42)


def null_pheno(h2):
    u = U @ (np.sqrt(s) * rng.normal(size=n))
    return np.sqrt(h2) * u / u.std() + np.sqrt(1 - h2) * rng.normal(size=n)


nulls = [null_pheno(h2) for h2 in [0.3] * 10 + [0.7] * 10]
rows = []
for F in [0.01, 0.02, 0.05, 0.10]:
    for level, M in matrices(F).items():
        UtM = U.T @ M
        plant = rng.choice(M.shape[1], 100)
        bases = [null_pheno(0.5) for _ in range(100)]
        for model in MODELS:
            real = [lambda_gc(pvals(M, UtM, Y[c].to_numpy(float), model)) for c in Y]
            nl = [pvals(M, UtM, y, model) for y in nulls]
            row = {"max_freq": F, "level": level, "model": model, "tests": M.shape[1],
                   "median_carriers": float(np.median(M.sum(0))),
                   "real_lambda_median": np.median(real), "real_lambda_max": np.max(real),
                   "null_lambda": np.median([lambda_gc(p) for p in nl]),
                   "null_tail": np.mean([np.mean(p < 1e-3) for p in nl]) / 1e-3}
            for b in (1, 2):
                hits = [pvals(M[:, [j]], UtM[:, [j]], y + b * M[:, j], model)[0] < 0.05 / M.shape[1]
                        for j, y in zip(plant, bases)]
                row[f"power_b{b}"] = np.mean(hits)
            rows.append(row)
            print(row, flush=True)
            pd.DataFrame(rows).to_csv(out, sep="\t", index=False, float_format="%.3g")
