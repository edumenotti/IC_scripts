"""Throwaway: which relatedness corrects both SNP and accessory-genome (ORF) structure? Gene burden, 2%, INT."""
import sys

import numpy as np
import pandas as pd
from scipy import optimize, stats as sps
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

from mutation_hunter import data, stats
from mutation_hunter.report import lambda_gc

z = np.load(sys.argv[1], allow_pickle=True)  # our own file
orf_path, out = sys.argv[2], sys.argv[3]
strains, gene, freq = z["strains"], z["gene"], z["freq"]
Y = data.load_phenotypes().loc[strains]
Dsnp = data.load_distance().loc[strains, strains]
Ksnp = stats.kinship(Dsnp).to_numpy()
Korf = stats.kinship(pd.read_csv(orf_path, sep="\t", index_col=0).loc[strains, strains]).to_numpy()
n = len(strains)
rng = np.random.default_rng(11)


def eig(K):
    s, U = np.linalg.eigh(K)
    return np.clip(s, 0, None), U


def top_pcs(K, k):
    s, U = eig(K)
    i = np.argsort(s)[::-1][:k]
    return U[:, i] * np.sqrt(s[i])


m = freq <= 0.02
Gv = z["G"][:, m].astype(float)
M = pd.DataFrame(Gv.T).groupby(gene[m]).max().T.to_numpy()
del Gv
car, npres = np.nansum(M, 0), (~np.isnan(M)).sum(0)
M = M[:, (car >= 3) & (car < npres)]
M[np.isnan(M)] = np.take(np.nanmean(M, 0), np.where(np.isnan(M))[1])

MODELS = {"snp": (Ksnp, 0), "snp+orfPC5": (Ksnp, 5), "snp+orfPC10": (Ksnp, 10), "snp+orfPC20": (Ksnp, 20)}
ORFPC = top_pcs(Korf, 20)
prep = {}
for name, (K, k) in MODELS.items():
    s, U = eig(K)
    X = np.column_stack([np.ones(n), ORFPC[:, :k]])
    prep[name] = (s, U, U.T @ X, U.T @ M)


def lmm_p(name, y):
    s, U, uX, UtM = prep[name]
    uy = U.T @ y

    def nll(ld):
        v = s + np.exp(ld)
        Xw = uX / v[:, None]
        r = uy - uX @ np.linalg.solve(uX.T @ Xw, Xw.T @ uy)
        return 0.5 * (n * np.log((r * r / v).sum()) + np.log(v).sum())

    w = 1 / np.sqrt(s + np.exp(optimize.minimize_scalar(nll, bounds=(-10, 10), method="bounded").x))
    return stats.assoc(UtM * w[:, None], uy * w, uX * w[:, None])[2]


def inorm(y):
    return sps.norm.ppf((sps.rankdata(y) - 0.5) / len(y))


def sim(K, h2):
    s, U = eig(K)
    u = U @ (np.sqrt(s) * rng.normal(size=n))
    return np.sqrt(h2) * u / u.std() + np.sqrt(1 - h2) * rng.normal(size=n)


def both(h2):  # structure from SNPs and accessory genome together
    return np.sqrt(h2) * (sim(Ksnp, 1) + sim(Korf, 1)) / np.sqrt(2) + np.sqrt(1 - h2) * rng.normal(size=n)


lab = hierarchy.fcluster(hierarchy.linkage(squareform(Dsnp.to_numpy(), checks=False), "average"), 30, "maxclust")


def permuted(y):
    y = y.copy()
    for c in np.unique(lab):
        i = np.flatnonzero(lab == c)
        y[i] = y[rng.permutation(i)]
    return y


hs = [0.3] * 10 + [0.7] * 10
sets = {"snp_sim": [sim(Ksnp, h) for h in hs], "orf_sim": [sim(Korf, h) for h in hs],
        "both_sim": [both(h) for h in hs], "perm_k30": [permuted(Y[c].to_numpy(float)) for c in Y],
        "real": [Y[c].to_numpy(float) for c in Y]}
plant = rng.choice(M.shape[1], 100)
rows = []
for model in MODELS:
    for name, ys in sets.items():
        ps = [lmm_p(model, inorm(y)) for y in ys]
        lam = [lambda_gc(p) for p in ps]
        row = {"model": model, "set": name, "lambda_median": np.median(lam), "lambda_max": np.max(lam),
               "tail_1e-3": np.mean([np.mean(p < 1e-3) for p in ps]) / 1e-3}
        if name == "both_sim":  # power: plant 1 SD in a random gene on a structured background
            hits = [lmm_p(model, inorm(ys[i % 20] + M[:, j]))[j] < 0.05 / M.shape[1]
                    for i, j in enumerate(plant)]
            row["power_b1"] = np.mean(hits)
        rows.append(row)
        print(row, flush=True)
pd.DataFrame(rows).to_csv(out, sep="\t", index=False, float_format="%.3g")
