"""Throwaway: carrier matrix for every rare variant (freq <= 0.10, >= 3 carriers) on the 970 analysed strains."""
import sys

import numpy as np
import pandas as pd

from mutation_hunter import data, variants

out = sys.argv[1]
Y, D = data.load_phenotypes(), data.load_distance()
strains = Y.index.intersection(D.index)
cols, genes, freqs = [], [], []
for gene in data.genes():
    seqs = data.load_gene(gene)
    rv = variants.rare_variants(gene, seqs, 0.10, 3)
    if rv.empty:
        continue
    G = variants.carrier_matrix(seqs, dict(zip(rv["variant"], rv["carriers"]))).reindex(strains)  # absent -> NaN
    cols.append(G.to_numpy(np.float32))
    genes += [gene] * len(rv)
    freqs += (rv["n_carriers"] / len(seqs)).tolist()
np.savez(out, G=np.hstack(cols), gene=np.array(genes), freq=np.array(freqs), strains=np.array(strains))
print(np.hstack(cols).shape)
