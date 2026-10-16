"""bcftools + cyvcf2 on a hand-written VCF with known consequences. Skipped when the tools are missing."""

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from mutation_hunter import build

pytest.importorskip("cyvcf2")
pytestmark = pytest.mark.skipif(not shutil.which("bcftools"), reason="bcftools not installed")
ANNOTATE = Path(__file__).parents[1] / "nextflow" / "bin" / "annotate.sh"
COMP = str.maketrans("ACGT", "TGCA")
SAMPLES = ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "SH"]  # SH is haploid


def orf(changes: dict[int, str]) -> str:
    """100 codons: ATG, 98 x GCT (Ala) with overrides, TAA. Protein length 99 aa."""
    codons = ["ATG"] + ["GCT"] * 98 + ["TAA"]
    for k, c in changes.items():
        codons[k - 1] = c
    return "".join(codons)


def reference() -> str:
    seq = list(np.random.default_rng(5).choice(list("ACGT"), 2000))
    seq[100:400] = orf({20: "TCA", 97: "TCA"})  # YAL901W, + strand, CDS 101-400
    seq[1000:1300] = orf({}).translate(COMP)[::-1]  # YAL902C, - strand, CDS 1001-1300
    return "".join(seq)


def pos(codon: int, base: int) -> int:  # 1-based chromosome position in YAL901W
    return 101 + 3 * (codon - 1) + (base - 1)


GT = "0/0:99:30"
RECORDS = [  # pos, ref, alt, genotypes for S1..S7, SH ; expected (class, aa_change prefix) per alt
    (pos(10, 1), "G", "A", ["0/1:99:30", "1/1:99:30", GT, GT, GT, GT, GT, "1/1:99:30"]),  # missense 10A>10T
    (pos(20, 2), "C", "A", ["1/1:99:30", GT, GT, GT, GT, GT, GT, GT]),  # stop_gained, early -> lof
    (pos(30, 3), "T", "TG", [GT, "1/1:99:30", GT, GT, GT, GT, GT, GT]),  # frameshift -> lof
    (pos(40, 3), "T", "C", [GT, GT, "1/1:99:30", GT, GT, GT, GT, GT]),  # synonymous
    (pos(50, 1), "G", "A,C", [GT, GT, GT, "1/2:99:30", GT, GT, GT, GT]),  # multiallelic: two missense rows
    (pos(97, 2), "C", "A", [GT, GT, GT, GT, "1/1:99:30", GT, GT, GT]),  # stop in last 5% -> other
    (pos(60, 1), "G", "A", ["1/1:99:30", GT, GT, GT, GT, "0/1:10:30", "0/1:99:3", "0/1:99:30"]),  # low GQ, low DP, haploid het
    (1273, "C", "T", [GT, GT, GT, GT, GT, "1/1:99:30", GT, GT]),  # missense on the - strand gene
]


@pytest.fixture
def built(tmp_path):
    ref = reference()
    for p, r, *_ in RECORDS:
        assert ref[p - 1:p - 1 + len(r)] == r, p
    (tmp_path / "ref.fa").write_text(f">chromosome1\n{ref}\n")
    (tmp_path / "ref.gff3").write_text(
        "##gff-version 3\n" + "".join(
            f"chromosome1\tt\tgene\t{a}\t{b}\t.\t{s}\t.\tID=gene:{g};biotype=protein_coding\n"
            f"chromosome1\tt\tmRNA\t{a}\t{b}\t.\t{s}\t.\tID=transcript:{g}_mRNA;Parent=gene:{g};biotype=protein_coding\n"
            f"chromosome1\tt\texon\t{a}\t{b}\t.\t{s}\t.\tParent=transcript:{g}_mRNA\n"
            f"chromosome1\tt\tCDS\t{a}\t{b}\t.\t{s}\t0\tParent=transcript:{g}_mRNA\n"
            for g, a, b, s in [("YAL901W", 101, 400, "+"), ("YAL902C", 1001, 1300, "-")]))
    header = ("##fileformat=VCFv4.2\n##contig=<ID=chromosome1,length=2000>\n"
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
              '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="q">\n'
              '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="d">\n'
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(SAMPLES) + "\n")
    body = "".join(f"chromosome1\t{p}\t.\t{r}\t{a}\t100\t.\t.\tGT:GQ:DP\t" + "\t".join(g) + "\n"
                   for p, r, a, g in sorted(RECORDS))
    (tmp_path / "in.vcf").write_text(header + body)
    subprocess.run(["bcftools", "view", "-Oz", "-o", tmp_path / "in.vcf.gz", tmp_path / "in.vcf"], check=True)
    subprocess.run(["bcftools", "index", "-t", tmp_path / "in.vcf.gz"], check=True)
    subprocess.run([ANNOTATE, tmp_path / "in.vcf.gz", "chromosome1", tmp_path / "ref.fa", tmp_path / "ref.gff3",
                    tmp_path / "out.bcf"], check=True)
    pd.DataFrame({"ploidy": [2] * 7 + [1]}, index=pd.Index(SAMPLES, name="strain")).to_csv(tmp_path / "strains.tsv", sep="\t")
    # Placement sites: one SNV, one multiallelic SNV, one indel (not a SNV) and one position absent from the VCF.
    pd.DataFrame({"chrom": "chromosome1", "pos": [pos(10, 1), pos(50, 1), pos(30, 3), 1999]}).to_csv(
        tmp_path / "sites.tsv", sep="\t", index=False)
    build.genotypes(tmp_path / "out.bcf", tmp_path / "strains.tsv", tmp_path / "ref.gff3", tmp_path / "chromosome1",
                    sites=tmp_path / "sites.tsv")
    var = pd.read_csv(tmp_path / "chromosome1.variants.tsv.gz", sep="\t", keep_default_na=False)
    return var, sparse.load_npz(tmp_path / "chromosome1.geno.npz").toarray(), tmp_path


def test_consequence_classes(built):
    var, _, _ = built
    got = list(zip(var["pos"], var["gene"], var["class"]))
    assert got == [
        (pos(10, 1), "YAL901W", "missense"), (pos(20, 2), "YAL901W", "lof"), (pos(30, 3), "YAL901W", "lof"),
        (pos(40, 3), "YAL901W", "synonymous"), (pos(50, 1), "YAL901W", "missense"),
        (pos(50, 1), "YAL901W", "missense"), (pos(60, 1), "YAL901W", "missense"), (pos(97, 2), "YAL901W", "other"),
        (1273, "YAL902C", "missense")]
    assert var["aa_change"][0].startswith("10A>10T") and var["protein_len"][0] == 99


def test_multiallelic_split_keeps_origin_and_genotypes(built):
    var, geno, _ = built
    rows = var.index[var["pos"] == pos(50, 1)]
    assert set(var.loc[rows, "alt"]) == {"A", "C"} and var.loc[rows, "old_rec"].str.contains("A,C").all()
    assert (geno[rows, 3] > 0).all()  # S4 (1/2) carries both


def test_genotype_codes_and_site_stats(built):
    var, geno, _ = built
    assert geno[0].tolist() == [1, 2, 0, 0, 0, 0, 0, 2]  # haploid hom alt stays
    assert (var.loc[0, "ac"], var.loc[0, "an"], var.loc[0, "n_het"]) == (5, 16, 1)


def test_quality_rules_set_missing(built):
    var, geno, _ = built
    row = var.index[var["pos"] == pos(60, 1)][0]
    assert geno[row].tolist() == [2, 0, 0, 0, 0, -1, -1, -1]  # GQ 10, DP 3, het in haploid
    assert var.loc[row, "call_rate"] == 5 / 8


def test_placement_uses_panel_codes_at_common_snvs(built):
    _, geno, tmp_path = built
    pl = np.load(tmp_path / "chromosome1.place.npz")
    assert pl["pos"].tolist() == [pos(10, 1), pos(50, 1), pos(50, 1)]  # indel and absent site skipped
    assert pl["ref"].tolist() == ["G", "G", "G"] and sorted(pl["alt"].tolist()[1:]) == ["A", "C"]
    assert pl["code"][0].tolist() == geno[0].tolist() == [1, 2, 0, 0, 0, 0, 0, 2]  # alt-allele codes, same as the panel
