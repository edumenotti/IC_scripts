# Calibration of the association test (2026-10-15)

> **INVALID — do not use these numbers.** Every kinship here was built by squaring a distance that is already
> squared (Hamming-type), which made the matrices non-positive-semidefinite (2.7% negative eigenvalue mass for
> SNPs, 18.5% for ORFs). Found in review on 2026-10-15 and fixed in `stats.kinship`. The method choices below
> (model, ORF PCs, threshold) must be re-derived on the VCF-based panel. Kept as a record of the procedure.

Why the defaults are: linear mixed model (EMMAX) with SNP kinship, 10 accessory-genome (ORF) coordinates as
fixed covariates, rank-based inverse normal transform (INT) of the phenotype, gene-level burden, rare = ≤ 2%.

λGC = median observed χ² / expected; 1 = calibrated. *Tail* = observed / expected rate of p < 10⁻³.
*Power* = share of 100 planted effects (1 SD) detected at Bonferroni 0.05.
Null phenotypes are simulated with no causal variant, so any λ above 1 is false signal.

## 1. Rarity threshold × level × model (`run.py` → `grid.tsv`)

Null phenotypes simulated from the SNP kinship (h² 0.3 and 0.7, 20 replicates); 35 real conditions.

| Model | null λ (all thresholds) | null tail |
|---|---|---|
| OLS, no correction | 4.1–10.3 | 145–320× |
| OLS + 5 PCs | 2.8–3.2 | 65–119× |
| OLS + 10 PCs | 2.1–2.5 | 28–49× |
| LMM | 0.91–1.00 | 0.5–1.8× |
| LMM + 5 SNP PCs | 0.96–1.21 | 0.7–2.1× |

LMM only:

| max_freq | level | tests | median carriers | real λ | null λ | null tail | power 1 SD | power 2 SD |
|---|---|---|---|---|---|---|---|---|
| 0.01 | variant | 164,131 | 4 | 1.52 | 0.91 | 0.52× | 0.00 | 0.20 |
| 0.01 | gene | 6,341 | 74 | 1.64 | 0.98 | 1.81× | 0.85 | 0.97 |
| 0.02 | variant | 203,279 | 5 | 1.54 | 0.98 | 0.61× | 0.03 | 0.40 |
| **0.02** | **gene** | **6,388** | **122** | **1.72** | **0.96** | **1.49×** | **0.93** | **0.98** |
| 0.05 | variant | 247,332 | 6 | 1.64 | 0.99 | 0.67× | 0.03 | 0.40 |
| 0.05 | gene | 6,432 | 224 | 1.90 | 0.96 | 0.98× | 0.92 | 0.98 |
| 0.10 | variant | 258,507 | 6 | 1.63 | 1.00 | 0.68× | 0.08 | 0.44 |
| 0.10 | gene | 6,439 | 272 | 1.93 | 0.96 | 0.89× | 0.95 | 1.00 |

Gene level has power; variant level genome-wide has almost none (4–6 carriers, ~200k tests). Real λ rises
with the threshold; power plateaus from 2%.

Caveat: these nulls come from the same kinship the LMM uses, so they favour it by construction.

## 2. Checks independent of the model's kinship (`run2.py` → `grid2.tsv`, 2%)

- **ORF simulation**: nulls from accessory-genome relatedness (ORF presence/absence distance). SNP-kinship
  LMM: λ 1.73, tail 61× (gene). The SNP kinship misses this structure.
- **Within-cluster permutation**: real phenotypes shuffled inside 10 or 30 SNP clusters (keeps lineage means
  and the real, skewed distributions). λ ≈ 1 at gene level; tail 1.6–1.9× raw, 1.3–1.4× with INT. At variant
  level the tail is 3–3.8× raw and 0.8–1.0× with INT → INT fixes the skew problem.

## 3. Which relatedness (`run3.py`, `run4.py` → `grid3.tsv`, `grid4.tsv`; gene, 2%, INT)

| Relatedness | λ: SNP / ORF / both / perm | tail: SNP / ORF / both / perm | real λ | power |
|---|---|---|---|---|
| SNP kinship | 1.03 / 1.44 / 1.38 / 0.96 | 1.8 / 21 / 15 / 1.7× | 1.68 | 0.97 |
| ORF kinship | 2.42 / 1.00 / 1.70 / 1.47 | 33 / 1.6 / 14 / 7.5× | 2.14 | — |
| mean of both | 1.29 / 0.96 / 1.15 / 1.17 | 3.9 / 1.4 / 2.9 / 2.4× | 1.72 | — |
| SNP + 5 ORF PCs | 1.05 / 1.19 / 1.06 / 1.02 | 1.8 / 4.3 / 7.5 / 1.5× | 1.52 | 0.97 |
| **SNP + 10 ORF PCs** | **1.03 / 1.15 / 1.11 / 1.05** | **1.4 / 2.2 / 1.8 / 1.4×** | **1.36** | **0.98** |
| SNP + 20 ORF PCs | 1.02 / 1.16 / 1.10 / 1.07 | 1.3 / 2.3 / 1.7 / 1.8× | 1.34 | 0.96 |

The remaining real λ (≈1.36) may be polygenic signal or structure not captured here; it is reported with
every run as `lambda_gc`.

## Reproduce

```bash
python docs/calibration/build.py G.npz            # carrier matrix, all rare variants ≤ 10%
python docs/calibration/run.py G.npz grid.tsv     # ~45 min
python docs/calibration/run2.py G.npz <1011DistanceMatrixBasedOnORFs.tab.gz> grid2.tsv
python docs/calibration/run3.py G.npz <ORF matrix> grid3.tsv
python docs/calibration/run4.py G.npz <ORF matrix> grid4.tsv
```
