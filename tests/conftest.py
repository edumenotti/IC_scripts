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


@pytest.fixture
def store(tmp_path):
    """Run build.assemble on hand-made per-chromosome parts; returns the store folder."""
    n = len(STORE_STRAINS)
    parts, inputs, ref = tmp_path / "parts", tmp_path / "inputs", tmp_path / "reference"
    for d in (parts, inputs, ref):
        d.mkdir()
    var = {"chromosome1": [("YAL901W", 128, "missense", [0, 1, 2, 0, -1, 2]),
                           ("YAL901W", 159, "lof", [2, 0, 0, 0, 0, 0])],
           "chromosome2": [("YBL902C", 1273, "missense", [0, 0, 0, 1, 1, 0])]}
    rng = np.random.default_rng(1)
    for chrom, rows in var.items():
        pd.DataFrame([{"chrom": chrom, "pos": pos, "ref": "G", "alt": "A", "old_rec": "", "gene": g,
                       "transcript": g + "_mRNA", "consequence": c, "class": c, "aa_change": "", "aa_pos": 0,
                       "protein_len": 99, **build.site_stats(np.array(code, np.int8))}
                      for g, pos, c, code in rows], columns=build.VARIANT_COLUMNS).to_csv(
            parts / f"{chrom}.variants.tsv.gz", sep="\t", index=False)
        sparse.save_npz(parts / f"{chrom}.geno.npz", sparse.csr_matrix(np.array([r[3] for r in rows], np.int8)))
        Z = rng.normal(size=(50, n))
        np.savez(parts / f"{chrom}.kin.npz", zz=Z.T @ Z, m=50)
    # Placement parts: chromosome1 site 7 is multiallelic (two rows) and gets dropped.
    np.savez(parts / "chromosome1.place.npz", pos=[5, 7, 7], ref=["C", "A", "A"], alt=["T", "G", "C"],
             code=np.array([[2, 1, 0, -1, 0, 2], [0] * 6, [1] * 6], np.int8))
    np.savez(parts / "chromosome2.place.npz", pos=[9], ref=["T"], alt=["G"], code=np.array([[0, 0, 1, 1, 2, 0]], np.int8))
    (ref / "ref.gff3").write_text(GFF)
    (ref / "ref.fa").write_text(">chromosome1\nACGT\n")
    (ref / "aliases.tsv").write_text("alias\tchrom\nI\tchromosome1\n")
    strains = pd.DataFrame({"isolate": STORE_STRAINS, "ploidy": [2, 2, 2, 2, 1, 2]}, index=pd.Index(STORE_STRAINS, name="strain"))
    strains.to_csv(tmp_path / "strains.tsv", sep="\t")

    pheno = pd.DataFrame({"COND1": rng.normal(size=n)}, index=STORE_STRAINS)
    (inputs / "phenoMatrix_35ConditionsNormalizedByYPD.tab.gz").write_bytes(gzip.compress(pheno.to_csv(sep="\t").encode()))
    cols = ["X1.EC1118_1F14_0012g", "X2.YAL901W", "X3.YAL005C_NumOfGenes_3"]
    pa = pd.DataFrame([[1, 1, 1]] * n, index=pd.Index(STORE_STRAINS, name="STD_name"), columns=cols)
    pa.iloc[4, 1] = 0
    for name, m in (("genesMatrix_PresenceAbsence.tab.gz", pa), ("genesMatrix_CopyNumber.tab.gz", pa * 1.5)):
        (inputs / name).write_bytes(gzip.compress(m.to_csv(sep="\t").encode()))
    (inputs / "allORFs_pangenome.fasta.gz").write_bytes(
        gzip.compress(b">1-EC1118_1F14_0012g\nATG\n>2-YAL901W\nATG\n>3-YAL005C_NumOfGenes_3\nATG\n"))
    gw = tmp_path / "gw"
    gw.mkdir()
    # Common SNPs of the public GWAS matrix: chromosome2:11 is absent from the VCF, chromosome 17 rows are CNVs.
    (gw / "1011GWAS_matrix.bim").write_text("1\trs0\t0\t5\t1\t2\n1\trs1\t0\t7\t1\t2\n2\trs2\t0\t9\t1\t2\n"
                                           "2\trs3\t0\t11\t1\t2\n17\tcnv\t0\t1\t1\t2\n")
    with tarfile.open(inputs / "1011GWASMatrix.tar.gz", "w:gz") as t:
        for f in gw.iterdir():
            t.add(f, arcname=f.name)

    out = tmp_path / "panel"
    build.assemble(parts, inputs, ref, tmp_path / "strains.tsv", out, {"fake": "md5"}, {"panel_version": "0.0.1-test"})
    return out
