import gzip
import json

import numpy as np
import pandas as pd
import pytest

import tarfile

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


def test_encode_masks_low_quality_and_haploid_hets():
    gt = np.array([0, 1, 3, 2, 1, 3, 0])
    gq = np.array([99, 99, 99, -1, 99, 10, 99])
    dp = np.array([30, 30, 30, 0, 30, 30, 3])
    haploid = np.array([0, 0, 0, 0, 1, 0, 0], bool)
    assert build.encode(gt, gq, dp, haploid).tolist() == [0, 1, 2, -1, -1, -1, -1]


def test_site_stats():
    s = build.site_stats(np.array([0, 1, 2, -1, 2], np.int8))
    assert (s["ac"], s["an"], s["n_het"], s["n_hom"], s["call_rate"], s["ref_is_minor"]) == (5, 8, 1, 2, 0.8, True)


def test_standardised_carriers_filters_rare_and_badly_called():
    common = np.array([0, 1, 2, 0, 0, 1, 0, 0, 0, 1] * 10, np.int8)
    z = build.standardised_carriers(common)
    assert z is not None and np.isclose(z.mean(), 0) and np.isclose(z.std(), 1)
    assert build.standardised_carriers(np.array([1] + [0] * 99, np.int8)) is None  # 1% carriers
    assert build.standardised_carriers(np.array([-1] * 10 + [0, 1] * 45, np.int8)) is None  # 90% called


def test_orf_table_uses_the_number_to_recover_real_names(tmp_path):
    fa = tmp_path / "orfs.fa.gz"
    fa.write_bytes(gzip.compress(b">1-EC1118_1F14_0012g\nA\n>7-YAL016C-A\nA\n>9-YAL005C_NumOfGenes_3\nA\n"))
    t = build.orf_table(["X1.EC1118_1F14_0012g", "X7.YAL016C.A", "X9.YAL005C_NumOfGenes_3"], fa)
    assert t["orf"].tolist() == ["EC1118_1F14_0012g", "YAL016C-A", "YAL005C_NumOfGenes_3"]
    assert t["gene"].tolist()[1:] == ["YAL016C-A", "YAL005C"] and pd.isna(t["gene"][0])
    assert t["family"].tolist() == [False, False, True]
    with pytest.raises(ValueError, match="no pangenome ORF"):
        build.orf_table(["X2.SOMETHING"], fa)


def test_gwas_sites_keeps_snps_and_drops_cnv_rows(tmp_path):
    (tmp_path / "m.bim").write_text("1\trs0\t0\t31605\t1\t2\n16\trs1\t0\t9\t1\t2\n17\tcnv\t0\t1\t1\t2\n")
    with tarfile.open(tmp_path / "g.tar.gz", "w:gz") as t:
        t.add(tmp_path / "m.bim", arcname="dir/m.bim")
    s = build.gwas_sites(tmp_path / "g.tar.gz")
    assert list(zip(s["chrom"], s["pos"])) == [("chromosome1", 31605), ("chromosome16", 9)]


def test_clone_groups():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(5, 400))
    X[1] = X[0] + rng.normal(0, 0.001, 400)  # near-identical pair
    K = X @ X.T / 400
    g = build.clone_groups(K)
    assert g[0] == g[1] and len(set(g)) == 4


def test_assemble_store_layout(store):
    m = json.loads((store / "manifest.json").read_text())
    assert m["format"] == 1 and m["panel_version"] == "0.0.1-test" and m["counts"]["variant_rows"] == 3
    var = pd.read_csv(store / "variants.tsv.gz", sep="\t")
    assert var["gene"].tolist() == ["YAL901W", "YAL901W", "YBL902C"] and var["row"].tolist() == [0, 1, 2]
    genes = pd.read_csv(store / "genes.tsv", sep="\t", index_col=0)
    assert genes.loc["YAL901W", ["first", "last"]].tolist() == [0, 1] and genes.loc["YBL903W", "first"] == -1
    assert genes.loc["YAL901W", "orf"] == "YAL901W"
    sites = pd.read_csv(store / "placement_sites.tsv", sep="\t")
    assert sites.to_dict("list") == {"chrom": ["chromosome1", "chromosome2"], "pos": [5, 9], "ref": ["C", "T"], "alt": ["T", "G"]}
    assert np.load(store / "placement.npy").tolist() == [[2, 0], [1, 0], [0, 1], [-1, 1], [0, 2], [2, 0]]  # strains x sites
    assert (m["counts"]["placement_sites"], m["counts"]["placement_multiallelic"], m["counts"]["placement_missing"]) == (2, 1, 1)
    assert sorted(p.name for p in (store / "kinship_loco").iterdir()) == ["chromosome1.npy", "chromosome2.npy"]
    listed = {line.split("  ", 1)[1] for line in (store / "SHA256SUMS").read_text().splitlines()}
    assert {"manifest.json", "genotypes/data.npy", "reference/ref.fa"} <= listed
    assert not store.with_name("panel.part").exists()


def test_assemble_rejects_strains_missing_from_matrices(store, tmp_path):
    inputs = tmp_path / "inputs"
    pa = pd.read_csv(inputs / "genesMatrix_PresenceAbsence.tab.gz", sep="\t", index_col=0).drop("S3")
    (inputs / "genesMatrix_PresenceAbsence.tab.gz").write_bytes(gzip.compress(pa.to_csv(sep="\t").encode()))
    with pytest.raises(ValueError, match="1 VCF strains missing"):
        build.assemble(tmp_path / "parts", inputs, tmp_path / "reference", tmp_path / "strains.tsv", tmp_path / "p2",
                       {}, {"panel_version": "x"})


def test_assemble_is_byte_reproducible(store, tmp_path):
    build.assemble(tmp_path / "parts", tmp_path / "inputs", tmp_path / "reference", tmp_path / "strains.tsv",
                   tmp_path / "again", {"fake": "md5"}, {"panel_version": "0.0.1-test"})
    differ = [line for line, other in zip((store / "SHA256SUMS").read_text().splitlines(),
                                          (tmp_path / "again" / "SHA256SUMS").read_text().splitlines()) if line != other]
    assert {line.split("  ")[1] for line in differ} <= {"manifest.json"}  # at most the creation time changes


def test_aliases_include_refseq_and_genbank(tmp_path):
    fa = tmp_path / "r.fa.gz"
    fa.write_bytes(gzip.compress(b">I x\nACGT\n"))
    gff = tmp_path / "r.gff3.gz"
    gff.write_bytes(gzip.compress(b"##gff-version 3\n"))
    build.reference(fa, gff, tmp_path / "out")
    aliases = pd.read_csv(tmp_path / "out/aliases.tsv", sep="\t").set_index("alias")["chrom"]
    assert aliases["NC_001133.9"] == "chromosome1" and aliases["BK006935.2"] == "chromosome1"
    assert aliases["BK006934.2"] == "chromosome8" and aliases["NC_001148.4"] == "chromosome16"


def test_manifest_records_clone_threshold_and_input_sources(store):
    m = json.loads((store / "manifest.json").read_text())
    assert m["params"]["clone_frac"] == 0.01
    assert set(m["inputs_source"]) == set(build.INPUTS) and all(m["inputs_source"].values())
