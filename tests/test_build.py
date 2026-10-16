import gzip
import json

import numpy as np
import pandas as pd
import pytest

from mutation_hunter import build


def test_parse_bcsq_keeps_protein_coding_and_systematic_name():
    s = ("missense|TFC3|YAL001C_mRNA|protein_coding|-|10A>10T|128G>A,"
         "non_coding|X|RNA1_ncRNA|ncRNA,synonymous|YAL002W|YAL002W_mRNA|protein_coding|+|40A|220T>C")
    assert build.parse_bcsq(s) == [("missense", "YAL001C", "10A>10T"), ("synonymous", "YAL002W", "40A")]
    assert build.parse_bcsq("") == []


@pytest.mark.parametrize("consequence,aa_pos,expected", [
    ("missense", 10, "missense"),
    ("stop_gained", 20, "lof"),
    ("stop_gained", 97, "other"),  # last 5% of a 99-aa protein
    ("frameshift&stop_retained", 30, "lof"),
    ("splice_donor", 0, "lof"),
    ("start_lost", 1, "lof"),
    ("stop_gained&inframe_insertion", 98, "other"),
    ("inframe_deletion", 12, "inframe"),
    ("synonymous", 40, "synonymous"),
    ("stop_lost", 99, "other"),
])
def test_classify(consequence, aa_pos, expected):
    assert build.classify(consequence, aa_pos, 99) == expected


def test_reference_renames_and_drops_mito(tmp_path):
    fa = tmp_path / "r.fa.gz"
    fa.write_bytes(gzip.compress(b">I dna:chromosome\nACGT\n>Mito dna:chromosome\nTTTT\n>XVI x\nGG\n"))
    gff = tmp_path / "r.gff3.gz"
    gff.write_bytes(gzip.compress(b"##gff-version 3\n##sequence-region   I 1 4\n##sequence-region   Mito 1 4\n"
                                  b"I\tsgd\tgene\t1\t4\t.\t+\t.\tID=gene:A\nMito\tsgd\tgene\t1\t4\t.\t+\t.\tID=gene:B\n"))
    build.reference(fa, gff, tmp_path / "out")
    assert (tmp_path / "out/ref.fa").read_text() == ">chromosome1\nACGT\n>chromosome16\nGG\n"
    g = (tmp_path / "out/ref.gff3").read_text()
    assert "chromosome1\tsgd\tgene" in g and "Mito" not in g and "##sequence-region chromosome1 1 4" in g
    assert (tmp_path / "out/ref.fa.fai").read_text() == "chromosome1\t4\t13\t4\t5\nchromosome16\t2\t32\t2\t3\n"
    aliases = pd.read_csv(tmp_path / "out/aliases.tsv", sep="\t").set_index("alias")["chrom"]
    assert aliases[["I", "chrI", "chr1", "1", "chromosome1"]].eq("chromosome1").all() and aliases["XVI"] == "chromosome16"


def test_check_inputs_rejects_changed_or_missing_files(tmp_path, monkeypatch):
    (tmp_path / "a.txt").write_text("hello")
    monkeypatch.setattr(build, "INPUTS", {"a.txt": build.file_hash(tmp_path / "a.txt")})
    assert build.check_inputs(tmp_path) == build.INPUTS
    (tmp_path / "a.txt").write_text("changed upstream")
    with pytest.raises(ValueError, match="a.txt: MD5"):
        build.check_inputs(tmp_path)
    monkeypatch.setattr(build, "INPUTS", {"b.txt": "0" * 32})
    with pytest.raises(FileNotFoundError, match="b.txt"):
        build.check_inputs(tmp_path)
