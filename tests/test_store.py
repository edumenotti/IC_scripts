import tarfile

import numpy as np
import pytest

from mutation_hunter import store as st


def test_panel_reads_gene_blocks(store):
    p = st.Panel(store)
    var, g = p.gene("YAL901W")
    assert var["pos"].tolist() == [128, 159]
    assert g.tolist() == [[0, 1, 2, 0, -1, 2], [2, 0, 0, 0, 0, 0]]
    var, g = p.gene("YBL903W")  # gene without variants
    assert var.empty and g.shape == (0, 6)
    with pytest.raises(KeyError, match="systematic"):
        p.gene("TPM1")


def test_panel_tables(store):
    p = st.Panel(store)
    assert list(p.strains.index) == ["S1", "S2", "S3", "S4", "S5", "SACE_YAA"]
    assert p.kinship().shape == (6, 6) and p.kinship("chromosome2").shape == (6, 6)
    assert not np.allclose(p.kinship(), p.kinship("chromosome2"))
    assert p.presence_absence().loc["S5", "YAL901W"] == 0 and p.copy_number().loc["S1", "YAL901W"] == 1.5
    assert list(p.phenotypes().columns) == ["COND1"]
    p.verify()


def test_panel_detects_tampering(store):
    p = st.Panel(store)
    (store / "strains.tsv").write_text("strain\tisolate\n")
    with pytest.raises(RuntimeError, match="strains.tsv"):
        p.verify()


def test_fetch_store_verifies_and_caches(store, tmp_path, monkeypatch):
    monkeypatch.setenv("MH_CACHE", str(tmp_path / "cache"))
    tar = tmp_path / "panel-0.0.1.tar.gz"
    with tarfile.open(tar, "w:gz") as t:
        t.add(store, arcname="panel-0.0.1")
    with pytest.raises(RuntimeError, match="SHA256 mismatch"):
        st.fetch_store(tar.as_uri(), "0" * 64, "0.0.1")
    assert list((tmp_path / "cache").iterdir()) == []  # no partial store left behind
    path = st.fetch_store(tar.as_uri(), st.sha256(tar), "0.0.1")
    assert (path / ".complete").exists() and st.Panel(path).gene("YAL901W")[1].shape == (2, 6)
    tar.unlink()  # cached: no second download
    assert st.fetch_store(tar.as_uri(), "unused", "0.0.1") == path


def test_default_path(monkeypatch, store):
    monkeypatch.delenv("MH_PANEL", raising=False)
    with pytest.raises(RuntimeError, match="no released panel"):
        st.default_path()
    monkeypatch.setenv("MH_PANEL", str(store))
    assert st.Panel().gene("YAL901W")[0]["pos"].tolist() == [128, 159]


def test_panel_rejects_unknown_format(store):
    m = (store / "manifest.json").read_text().replace('"format": 1', '"format": 2')
    (store / "manifest.json").write_text(m)
    with pytest.raises(ValueError, match="format 2"):
        st.Panel(store)


def test_panel_rejects_partial_copy(store):
    (store / "kinship_loco" / "chromosome2.npy").unlink()
    with pytest.raises(RuntimeError, match="incomplete.*kinship_loco/chromosome2.npy"):
        st.Panel(store)


def test_panel_rejects_truncated_file(store):
    f = store / "presence_absence.npy"
    f.write_bytes(f.read_bytes()[:-8])
    with pytest.raises(RuntimeError, match="incomplete.*presence_absence.npy"):
        st.Panel(store)
