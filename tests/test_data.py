import pytest

from mutation_hunter import data


def test_fetch_verifies_and_caches(panel):
    path = data.fetch(data.PHENOTYPES)
    assert path.parent == data.cache_dir()
    (panel / data.PHENOTYPES).unlink()  # second call must not hit the network
    assert data.fetch(data.PHENOTYPES) == path


def test_fetch_rejects_md5_mismatch(panel, monkeypatch):
    monkeypatch.setitem(data.FILES, data.PHENOTYPES, "0" * 32)
    with pytest.raises(RuntimeError, match="MD5 mismatch"):
        data.fetch(data.PHENOTYPES)
    assert list(data.cache_dir().iterdir()) == []  # no corrupt or partial file left


def test_load_gene_keeps_underscored_ids_and_absent_strains(panel):
    seqs = data.load_gene("YAL002W")
    assert "SACE_YBA" in seqs and "A37" not in seqs
    assert data.genes() == ["YAL001C", "YAL002W"]
    with pytest.raises(KeyError, match="YXX999"):
        data.load_gene("YXX999")


def test_tables_load(panel):
    assert list(data.load_phenotypes().columns) == ["COND1", "COND2"]
    dist = data.load_distance()
    assert dist.shape == (40, 40) and (dist.index == dist.columns).all()
