"""Rare amino-acid variants and carrier matrices from the aligned 1011 proteome."""

import re

import numpy as np
import pandas as pd

QUERY = re.compile(r"^(?P<gene>[A-Za-z0-9-]+)(?::(?P<ref>[A-Z*])(?P<pos>\d+)(?P<alt>[A-Z*]))?$")
MISSING = set("X-")


def parse_query(query: str) -> tuple[str, tuple[str, int, str] | None]:
    """'YNL079C' -> ('YNL079C', None); 'YNL079C:A123V' -> ('YNL079C', ('A', 123, 'V'))."""
    m = QUERY.match(query)
    if not m:
        raise ValueError(f"bad query {query!r}; expected GENE or GENE:REFposALT, e.g. YNL079C:A123V")
    if m["pos"] is None:
        return m["gene"], None
    return m["gene"], (m["ref"], int(m["pos"]), m["alt"])


def _matrix(seqs: dict[str, str]) -> np.ndarray:
    if len({len(s) for s in seqs.values()}) != 1:
        raise ValueError("sequences are not aligned (different lengths)")
    return np.array([list(s) for s in seqs.values()])


def _major(col: np.ndarray) -> str:
    values, counts = np.unique(col, return_counts=True)
    return values[counts.argmax()]


def rare_variants(gene: str, seqs: dict[str, str], max_freq: float, min_carriers: int) -> pd.DataFrame:
    """One row per rare allele. ref is the panel's major allele; pos is 1-based.
    Frequency is over strains that have the gene."""
    strains = np.array(list(seqs))
    A = _matrix(seqs)
    rows = []
    for j in np.flatnonzero((A != A[0]).any(axis=0)):
        col = A[:, j]
        ref = _major(col)
        values, counts = np.unique(col, return_counts=True)
        for alt, n in zip(values, counts):
            if alt != ref and alt not in MISSING and min_carriers <= n <= max_freq * len(col):
                rows.append((gene, f"{ref}{j + 1}{alt}", j + 1, ref, alt, int(n), strains[col == alt].tolist()))
    return pd.DataFrame(rows, columns=["gene", "variant", "pos", "ref", "alt", "n_carriers", "carriers"])


def carriers_of(gene: str, seqs: dict[str, str], ref: str, pos: int, alt: str) -> list[str]:
    """Strains carrying alt at pos. ref must be the panel's major allele there."""
    A = _matrix(seqs)
    if not 1 <= pos <= A.shape[1]:
        raise ValueError(f"{gene} has {A.shape[1]} positions; {pos} is out of range")
    col = A[:, pos - 1]
    if (major := _major(col)) != ref:
        raise ValueError(f"{gene} position {pos}: the major allele is {major}, not {ref}")
    found = np.array(list(seqs))[col == alt].tolist()
    if not found:
        raise ValueError(f"no strain carries {gene}:{ref}{pos}{alt}")
    return found


def carrier_matrix(seqs: dict[str, str], carriers: dict[str, list[str]]) -> pd.DataFrame:
    """Strains-with-the-gene x tests, 1.0 = carrier. Strains lacking the gene are not rows (missing, not reference)."""
    return pd.DataFrame({name: pd.Series(1.0, index=c) for name, c in carriers.items()}, index=list(seqs)).fillna(0.0)
