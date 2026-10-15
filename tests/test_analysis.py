import pytest

from mutation_hunter import analysis, data


def test_phewas_variant_finds_planted_signal(panel):
    res, rare = analysis.phewas("YAL001C:A3V", max_freq=0.10)
    assert list(res.columns) == analysis.COLUMNS
    cond1 = res.set_index("phenotype").loc["COND1"]
    assert cond1.n_carriers == 4 and cond1.beta > 0.5 and cond1.q < 0.05
    assert res.set_index("phenotype").loc["COND2"].q > 0.05
    assert rare["variant"].tolist() == ["A3V"]


def test_phewas_gene_burden_excludes_absent_strains(panel):
    res, _ = analysis.phewas("YAL002W", max_freq=0.10)
    assert (res["n"] == 39).all()  # A37 lacks the gene: missing, not a non-carrier


def test_phewas_gene_without_rare_variants_errors(panel):
    with pytest.raises(ValueError, match="no rare variant"):
        analysis.phewas("YAL002W", max_freq=0.10, min_carriers=5)


def test_gwas_gene_level(panel):
    res = analysis.gwas(data.load_phenotypes()["COND1"], level="gene", max_freq=0.10)
    assert res["test"].tolist() == ["YAL001C", "YAL002W"]
    assert res.set_index("test").loc["YAL001C"].q < 0.05


def test_gwas_variant_level_respects_thresholds(panel):
    res = analysis.gwas(data.load_phenotypes()["COND1"], level="variant", max_freq=0.10, min_carriers=1)
    assert len(res) == 2 and res["p"].notna().all()
    res = analysis.gwas(data.load_phenotypes()["COND1"], level="variant", max_freq=0.5, min_carriers=4)
    assert res["test"].tolist() == ["A3V"]
