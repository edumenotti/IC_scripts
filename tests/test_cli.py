import json

import pandas as pd

from mutation_hunter import cli


def test_fetch_lists_conditions(panel, capsys):
    assert cli.main(["fetch"]) == 0
    assert capsys.readouterr().out.split() == ["COND1", "COND2"]


def test_phewas_writes_three_files(panel, tmp_path):
    out = tmp_path / "r"
    assert cli.main(["phewas", "YAL001C:A3V", "--max-freq", "0.1", "-o", str(out)]) == 0
    assert len(pd.read_csv(f"{out}.tsv", sep="\t")) == 2
    prov = json.loads(open(f"{out}.json").read())
    assert prov["params"]["query"] == "YAL001C:A3V" and "1011proteome.tar.gz" in prov["data_md5"]
    assert "<img" in open(f"{out}.html").read()


def test_gwas_then_fdr(panel, tmp_path):
    for c in ["COND1", "COND2"]:
        assert cli.main(["gwas", c, "--max-freq", "0.1", "-o", str(tmp_path / c)]) == 0
    assert cli.main(["fdr", str(tmp_path / "COND1.tsv"), str(tmp_path / "COND2.tsv"), "-o", str(tmp_path / "all")]) == 0
    assert len(pd.read_csv(tmp_path / "all.tsv", sep="\t")) == 4


def test_user_phenotype_file(panel, tmp_path):
    f = tmp_path / "mine.tsv"
    f.write_text("strain\tgrowth\nA00\t3\nA01\t3\nA02\t3\nSACE_YBA\t3\n" + "".join(f"A{i:02d}\t1\n" for i in range(3, 30)))
    assert cli.main(["gwas", str(f), "--max-freq", "0.1", "-o", str(tmp_path / "mine")]) == 0


def test_bad_input_exits_1_with_message(panel, tmp_path, capsys):
    f = tmp_path / "other.tsv"
    f.write_text("strain\tgrowth\nZZZ\t1\n")
    assert cli.main(["gwas", str(f), "-o", str(tmp_path / "x")]) == 1
    assert cli.main(["phewas", "YAL001C:V3A", "-o", str(tmp_path / "x")]) == 1
    assert cli.main(["gwas", "NOPE", "-o", str(tmp_path / "x")]) == 1
    err = capsys.readouterr().err
    assert "no strain ID in common" in err and "major allele is A" in err and "unknown phenotype" in err
