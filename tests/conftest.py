"""Fixtures: a 40-strain fake proteome panel (baseline analysis) and a tiny built panel store (no tools)."""

import gzip
import io
import tarfile

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from mutation_hunter import build, data

STRAINS = [f"A{i:02d}" for i in range(38)] + ["SACE_YBA", "SACE_YCL"]
CARRIERS = ["A00", "A01", "A02", "SACE_YBA"]  # carry YAL001C:A3V, 4/40 = 10%
ABSENT = ["A37"]  # lacks YAL002W


def _fasta(gene, seqs):
    return "".join(f">{s}_{gene}_NAME_1 chromosome1:1-30 + 1\n{q}\n" for s, q in seqs.items()).encode()


@pytest.fixture
def panel(tmp_path, monkeypatch):
    rng = np.random.default_rng(0)
    src = tmp_path / "upstream"
    src.mkdir()

    hi = [s in CARRIERS for s in STRAINS]
    rows = [f"{s}\t{2.0 + rng.normal(0, .1) if h else 1.0 + rng.normal(0, .1):.4f}\t{rng.normal(1, .1):.4f}"
            for s, h in zip(STRAINS, hi)]
    (src / data.PHENOTYPES).write_bytes(gzip.compress(("\tCOND1\tCOND2\n" + "\n".join(rows) + "\n").encode()))

    d = rng.uniform(0.1, 0.2, (40, 40))
    d = (d + d.T) / 2
    np.fill_diagonal(d, 0)
    lines = ["STD\t" + "\t".join(STRAINS)] + [s + "\t" + "\t".join(f"{x:.4f}" for x in r) for s, r in zip(STRAINS, d)]
    (src / data.DISTANCE).write_bytes(gzip.compress(("\n".join(lines) + "\n").encode()))

    o = rng.integers(20, 120, (40, 40))
    o = (o + o.T) // 2
    np.fill_diagonal(o, 0)
    lines = ["STD\t" + "\t".join(STRAINS)] + [s + "\t" + "\t".join(map(str, r)) for s, r in zip(STRAINS, o)]
    (src / data.ORF_DISTANCE).write_bytes(gzip.compress(("\n".join(lines) + "\n").encode()))

    g1 = {s: "MKVLE*" if s in CARRIERS else "MKALE*" for s in STRAINS}
    g2 = {s: "MSVP*" if s in ("A10", "A11", "A12") else "MSTP*" for s in STRAINS if s not in ABSENT}
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as t:
        for gene, seqs in [("YAL001C", g1), ("YAL002W", g2)]:
            raw = _fasta(gene, seqs)
            info = tarfile.TarInfo(f"Proteome_1011/{gene}.fasta")
            info.size = len(raw)
            t.addfile(info, io.BytesIO(raw))
    (src / data.PROTEOME).write_bytes(buf.getvalue())

    monkeypatch.setattr(data, "FILES", {n: data.md5(src / n) for n in data.FILES})
    monkeypatch.setenv("MH_BASE_URL", src.as_uri() + "/")
    monkeypatch.setenv("MH_CACHE", str(tmp_path / "cache"))
    return src


# --- tiny built panel ------------------------------------------------------------------------------------

STORE_STRAINS = ["S1", "S2", "S3", "S4", "S5", "SACE_YAA"]
GFF = """##gff-version 3
chromosome1\tt\tgene\t101\t400\t.\t+\t.\tID=gene:YAL901W;biotype=protein_coding
chromosome1\tt\tmRNA\t101\t400\t.\t+\t.\tID=transcript:YAL901W_mRNA;Parent=gene:YAL901W;biotype=protein_coding
chromosome1\tt\tCDS\t101\t400\t.\t+\t0\tParent=transcript:YAL901W_mRNA
chromosome2\tt\tgene\t1001\t1300\t.\t-\t.\tID=gene:YBL902C;biotype=protein_coding
chromosome2\tt\tmRNA\t1001\t1300\t.\t-\t.\tID=transcript:YBL902C_mRNA;Parent=gene:YBL902C;biotype=protein_coding
chromosome2\tt\tCDS\t1001\t1300\t.\t-\t0\tParent=transcript:YBL902C_mRNA
chromosome2\tt\tgene\t1500\t1800\t.\t+\t.\tID=gene:YBL903W;biotype=protein_coding
"""


def write_bed(prefix, dosage: np.ndarray, chroms: list[int], samples: list[str]) -> None:
    """Minimal PLINK 1 SNP-major writer (dosage = A1 copies, -1 missing), for tests only."""
    code = {2: 0, -1: 1, 1: 2, 0: 3}
    n = dosage.shape[1]
    out = bytearray([0x6C, 0x1B, 0x01])
    for row in dosage:
        b = [code[int(x)] for x in row] + [0] * (-n % 4)
        out += bytes(b[i] | b[i + 1] << 2 | b[i + 2] << 4 | b[i + 3] << 6 for i in range(0, len(b), 4))
    open(f"{prefix}.bed", "wb").write(bytes(out))
    pd.DataFrame({"c": chroms, "id": [f"rs{i}" for i in range(len(chroms))], "cm": 0, "pos": range(1, len(chroms) + 1),
                  "a1": "A", "a2": "G"}).to_csv(f"{prefix}.bim", sep="\t", header=False, index=False)
    open(f"{prefix}.fam", "w").write("".join(f"{s} {s} 0 0 0 -9\n" for s in samples))
