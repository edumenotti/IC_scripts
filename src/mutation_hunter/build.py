"""Build the panel store from the public 1011 genomes data. Maintainers only.

Steps are called by nextflow/build.nf; bcftools does normalisation and consequence annotation upstream.
Needs the optional extra `build` (cyvcf2, xlrd)."""

import argparse
import gzip
import hashlib
import json
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import connected_components

from . import __version__

ROMAN = "I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI".split()
CHROMS = [f"chromosome{i}" for i in range(1, 17)]
RENAME = dict(zip(ROMAN, CHROMS))
LOF = {"stop_gained", "frameshift", "splice_acceptor", "splice_donor", "start_lost"}
INFRAME = {"inframe_insertion", "inframe_deletion"}
SYSTEMATIC = re.compile(r"^Y[A-P][LR]\d{3}[WC](-[A-Z])?$")
VARIANT_COLUMNS = ["chrom", "pos", "ref", "alt", "old_rec", "gene", "transcript", "consequence", "class",
                   "aa_change", "aa_pos", "protein_len", "ac", "an", "n_het", "n_hom", "call_rate", "ref_is_minor"]
# Upstream inputs and their MD5 (md5.txt of 1002genomes on 2025-04-01, Ensembl release 115, Nature supplement).
INPUTS = {
    "1011Matrix.gvcf.gz": "42478e3e9dff4bd46993d82d8eab40d4",
    "Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz": "5ae49d2d22cd977abdd6b054b94c544d",
    "Saccharomyces_cerevisiae.R64-1-1.115.gff3.gz": "17ed999c1853e03562f89bbcaefa80b6",
    "phenoMatrix_35ConditionsNormalizedByYPD.tab.gz": "479862335552a050e46de163d1ee5f72",
    "genesMatrix_PresenceAbsence.tab.gz": "b3b96719de226eaa277589d8d01f9202",
    "genesMatrix_CopyNumber.tab.gz": "c36cdf023c554580fe5b67101742fc09",
    "allORFs_pangenome.fasta.gz": "4f1ebba402ae64e52f2af72fa6942f38",
    "1011GWASMatrix.tar.gz": "c68c2ee5d9f4b5c85b8436a4a2f31ba5",
    "41586_2018_30_MOESM3_ESM.xls": "23d0ff1ec5c13a6b2e20a8076eb588c1",
}


def file_hash(path, algo="md5") -> str:
    h = hashlib.new(algo)
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check_inputs(folder) -> dict:
    """Verify every upstream input in folder against its pinned MD5; return {name: md5}."""
    folder = Path(folder)
    for name, md5 in INPUTS.items():
        if not (folder / name).exists():
            raise FileNotFoundError(f"missing input {name} in {folder}")
        if (got := file_hash(folder / name)) != md5:
            raise ValueError(f"{name}: MD5 {got} != pinned {md5}; upstream changed or file is corrupt")
    return dict(INPUTS)


def reference(fasta_gz, gff_gz, out) -> None:
    """R64 without Mito, chromosomes renamed I..XVI -> chromosome1..16 (the panel's names), plus aliases.tsv."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    keep = False
    with gzip.open(fasta_gz, "rt") as f, open(out / "ref.fa", "w") as o:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if keep := name in RENAME:
                    o.write(f">{RENAME[name]}\n")
            elif keep:
                o.write(line)
    with gzip.open(gff_gz, "rt") as f, open(out / "ref.gff3", "w") as o:
        for line in f:
            if line.startswith("##sequence-region"):
                p = line.split()
                if p[1] in RENAME:
                    o.write(" ".join([p[0], RENAME[p[1]], *p[2:]]) + "\n")
            elif line.startswith("#"):
                o.write(line)
            elif (col := line.split("\t", 1))[0] in RENAME:
                o.write(RENAME[col[0]] + "\t" + col[1])
    rows = [(a, c) for r, c in RENAME.items() for a in (r, f"chr{r}", c, c.replace("chromosome", "chr"),
                                                         c.replace("chromosome", ""))]
    pd.DataFrame(rows, columns=["alias", "chrom"]).to_csv(out / "aliases.tsv", sep="\t", index=False)
    write_fai(out / "ref.fa")


def write_fai(fa) -> None:
    """samtools-style .fai, written once here so parallel bcftools jobs never race to create it."""
    rows, offset = [], 0
    with open(fa, "rb") as f:
        for line in f:
            if line.startswith(b">"):
                rows.append([line[1:].split()[0].decode(), 0, offset + len(line), 0, 0])
            else:
                row = rows[-1]
                if not row[3]:
                    row[3], row[4] = len(line.rstrip(b"\n")), len(line)
                row[1] += len(line.rstrip(b"\n"))
            offset += len(line)
    Path(f"{fa}.fai").write_text("".join("\t".join(map(str, r)) + "\n" for r in rows))


def protein_lengths(gff) -> dict[str, int]:
    """Transcript -> protein length in amino acids (CDS length / 3, stop codon excluded)."""
    cds = {}
    for line in open(gff):
        col = line.split("\t")
        if len(col) > 8 and col[2] == "CDS":
            t = re.search(r"Parent=transcript:([^;\n]+)", col[8])[1]
            cds[t] = cds.get(t, 0) + int(col[4]) - int(col[3]) + 1
    return {t: n // 3 - 1 for t, n in cds.items()}


def gene_spans(gff) -> pd.DataFrame:
    rows = []
    for line in open(gff):
        col = line.split("\t")
        if len(col) > 8 and col[2] == "gene" and "biotype=protein_coding" in col[8]:
            rows.append((re.search(r"ID=gene:([^;\n]+)", col[8])[1], col[0], int(col[3]), int(col[4])))
    return pd.DataFrame(rows, columns=["gene", "chrom", "start", "end"]).set_index("gene")


def parse_bcsq(bcsq: str) -> list[tuple[str, str, str]]:
    """BCSQ -> [(consequence, systematic gene, amino-acid change)] for protein-coding transcripts only."""
    out = []
    for item in bcsq.split(","):
        f = item.split("|")
        if len(f) >= 4 and f[2].endswith("_mRNA") and f[3] == "protein_coding":
            out.append((f[0], f[2][:-5], f[5] if len(f) > 5 else ""))
    return out


def classify(consequence: str, aa_pos: int, protein_len: int) -> str:
    """missense / lof / inframe / synonymous / other. Truncations in the last 5% of the protein are not LoF."""
    parts = set(consequence.split("&"))
    if parts & LOF:
        late = parts & {"stop_gained", "frameshift"} and not parts & {"start_lost", "splice_acceptor", "splice_donor"}
        return "other" if late and protein_len and aa_pos > 0.95 * protein_len else "lof"
    if parts & INFRAME:
        return "inframe"
    if "missense" in parts:
        return "missense"
    return "synonymous" if "synonymous" in parts else "other"


def encode(gt_types, gq, dp, haploid, min_gq=20, min_dp=5) -> np.ndarray:
    """cyvcf2 gt_types (0 ref, 1 het, 2 unknown, 3 hom alt) -> int8: -1 missing, 0 ref, 1 het, 2 hom alt.
    Low GQ/DP and heterozygous calls in haploid strains become missing."""
    code = np.select([gt_types == 1, gt_types == 3], [1, 2], 0).astype(np.int8)
    bad = (gt_types == 2) | (gq < min_gq) | (dp < min_dp) | ((gt_types == 1) & haploid)
    code[bad] = -1
    return code


def site_stats(code: np.ndarray) -> dict:
    called = code >= 0
    n_het, n_hom, an = int((code == 1).sum()), int((code == 2).sum()), int(2 * called.sum())
    ac = n_het + 2 * n_hom
    return {"ac": ac, "an": an, "n_het": n_het, "n_hom": n_hom, "call_rate": called.mean(),
            "ref_is_minor": bool(an and ac / an > 0.5)}


def standardised_carriers(code: np.ndarray, min_freq=0.05, min_call=0.95) -> np.ndarray | None:
    """Common-SNP carrier vector standardised for the GRM (missing -> 0, i.e. mean-imputed), or None."""
    called = code >= 0
    if called.mean() < min_call:
        return None
    f = (code[called] > 0).mean()
    if not min_freq <= f <= 1 - min_freq:
        return None
    return np.where(called, (code > 0) - f, 0.0) / np.sqrt(f * (1 - f))


def genotypes(bcf, strains_tsv, gff, out_prefix, min_gq=20, min_dp=5) -> None:
    """One chromosome: <out>.variants.tsv.gz, <out>.geno.npz (variant rows x strains) and <out>.kin.npz
    (sum of Z'Z over common SNPs and their count, for the GRM and its LOCO versions)."""
    from cyvcf2 import VCF

    strains = pd.read_csv(strains_tsv, sep="\t", index_col=0)
    vcf = VCF(str(bcf), strict_gt=True)
    if list(vcf.samples) != list(strains.index):
        raise ValueError("VCF samples and strains.tsv differ in content or order")
    haploid = (strains["ploidy"] == 1).to_numpy()
    plen = protein_lengths(gff)
    rows, mats, block = [], [], []
    zz, m = np.zeros((len(strains),) * 2), 0
    for rec in vcf:
        code = encode(rec.gt_types, rec.gt_quals, rec.format("DP")[:, 0], haploid, min_gq, min_dp)
        if len(rec.REF) == 1 and len(rec.ALT[0]) == 1 and (z := standardised_carriers(code)) is not None:
            block.append(z)
            if len(block) == 2000:
                Z = np.array(block)
                zz += Z.T @ Z
                m += len(block)
                block = []
        ann = parse_bcsq(rec.INFO.get("BCSQ") or "")
        if not ann:
            continue
        stats_ = site_stats(code)
        nz = sparse.csr_matrix(code)
        for consequence, gene, aa in ann:
            aa_pos = int(mt[1]) if (mt := re.match(r"(\d+)", aa)) else 0
            rows.append({"chrom": rec.CHROM, "pos": rec.POS, "ref": rec.REF, "alt": rec.ALT[0],
                         "old_rec": rec.INFO.get("OLD_REC") or "", "gene": gene, "transcript": gene + "_mRNA",
                         "consequence": consequence, "class": classify(consequence, aa_pos, plen.get(gene + "_mRNA", 0)),
                         "aa_change": aa, "aa_pos": aa_pos, "protein_len": plen.get(gene + "_mRNA", 0), **stats_})
            mats.append(nz)  # ponytail: overlapping genes duplicate the genotype row; fine at yeast scale
    if block:
        Z = np.array(block)
        zz += Z.T @ Z
        m += len(block)
    pd.DataFrame(rows, columns=VARIANT_COLUMNS).to_csv(f"{out_prefix}.variants.tsv.gz", sep="\t", index=False)
    geno = sparse.vstack(mats, format="csr") if mats else sparse.csr_matrix((0, len(strains)), dtype=np.int8)
    sparse.save_npz(f"{out_prefix}.geno.npz", geno.astype(np.int8))
    np.savez(f"{out_prefix}.kin.npz", zz=zz, m=m)


def strain_table(xls, samples: list[str], phenotypes_gz) -> pd.DataFrame:
    """Peter 2018 Table S1, in VCF sample order."""
    s1 = pd.read_excel(xls, "Table S1", header=3)
    s1 = s1[s1["Standardized name"].notna()].set_index("Standardized name")
    pheno = pd.read_csv(phenotypes_gz, sep="\t", index_col=0)
    t = pd.DataFrame({
        "isolate": s1["Isolate name"],
        "ploidy": pd.to_numeric(s1["Ploidy"], errors="coerce").astype("Int64"),  # 'Xn' -> unknown
        "zygosity": s1["Zygosity"],
        "aneuploidy": s1["Aneuploidies"],
        "clade": s1["Clades"].astype(str).str.strip(),
        "ecological_origin": s1["Ecological origins"],
        "mta_restricted": s1["MTA (no distribution)"].eq("yes"),
    }).reindex(samples)
    if t["isolate"].isna().any():
        raise ValueError(f"strains missing from Table S1: {t.index[t['isolate'].isna()].tolist()[:5]}")
    t.index.name = "strain"
    t["has_phenotype"] = t.index.isin(pheno.index)
    return t


def orf_table(pa_columns, pangenome_fasta_gz) -> pd.DataFrame:
    """Matrix column 'X7051.YOR083W' -> pangenome ORF 7051 -> real name and systematic gene (if any).
    The number links the R-mangled column to the FASTA header '7051-YOR083W'."""
    names = {}
    with gzip.open(pangenome_fasta_gz, "rt") as f:
        for line in f:
            if line.startswith(">"):
                num, name = line[1:].strip().split("-", 1)
                names[int(num)] = name
    num = pd.Series(pa_columns).str.extract(r"^X(\d+)\.")[0].astype(int)
    orf = num.map(names)
    if orf.isna().any():
        raise ValueError(f"{orf.isna().sum()} matrix columns have no pangenome ORF")
    base = orf.str.replace(r"_NumOfGenes_\d+$", "", regex=True)
    return pd.DataFrame({"column": list(pa_columns), "orf": orf, "family": orf != base,
                         "gene": base.where(base.str.match(SYSTEMATIC.pattern))})


def read_plink_bed(prefix) -> tuple[pd.DataFrame, np.ndarray]:
    """PLINK 1 SNP-major .bed -> (bim table, int8 sites x samples: alt-allele count 0/1/2, -1 missing)."""
    bim = pd.read_csv(f"{prefix}.bim", sep=r"\s+", header=None, names=["chrom", "id", "cm", "pos", "a1", "a2"])
    n = sum(1 for _ in open(f"{prefix}.fam"))
    raw = np.fromfile(f"{prefix}.bed", dtype=np.uint8)
    if raw[:3].tolist() != [0x6C, 0x1B, 0x01]:
        raise ValueError("not a SNP-major PLINK 1 .bed")
    packed = raw[3:].reshape(len(bim), -1)
    bits = np.stack([(packed >> s) & 3 for s in (0, 2, 4, 6)], axis=2).reshape(len(bim), -1)[:, :n]
    # PLINK codes: 0 hom A1, 1 missing, 2 het, 3 hom A2; count A1 copies.
    lut = np.array([2, -1, 1, 0], dtype=np.int8)
    return bim, lut[bits]
