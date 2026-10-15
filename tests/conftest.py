"""A 40-strain fake panel in the exact upstream formats, served from file:// (no internet)."""

import gzip
import io
import tarfile

import numpy as np
import pytest

from mutation_hunter import data

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
