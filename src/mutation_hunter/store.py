"""Read a panel store built by mutation_hunter.build, and download/cache a released one."""

import hashlib
import json
import os
import shutil
import tarfile
import tempfile
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

FORMAT = 1
# The released store. Unset until the panel is published (needs the Schacherer/Liti group's permission).
RELEASE = {"version": None, "url": None, "sha256": None}


def cache_root() -> Path:
    default = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "mutation_hunter"
    return Path(os.environ.get("MH_CACHE", default))


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch_store(url: str, digest: str, version: str) -> Path:
    """Download a store tarball (one top-level folder) once into cache/<version>, verifying its SHA256.
    Extraction happens in a private temp dir and is renamed into place, so concurrent callers never see
    a partial store; '.complete' marks a finished one."""
    final = cache_root() / version
    if (final / ".complete").exists():
        return final
    cache_root().mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(dir=cache_root(), prefix=f".{version}-"))
    try:
        tar = work / "store.tar.gz"
        with urllib.request.urlopen(url, timeout=60) as r, open(tar, "wb") as f:
            shutil.copyfileobj(r, f)
        if sha256(tar) != digest:
            raise RuntimeError(f"panel {version}: SHA256 mismatch; download corrupt or release changed")
        with tarfile.open(tar) as t:
            t.extractall(work / "x", filter="data")
        (inner,) = (work / "x").iterdir()
        (inner / ".complete").touch()
        try:
            inner.rename(final)
        except OSError:  # another process finished first
            if not (final / ".complete").exists():
                raise
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return final


def default_path() -> Path:
    """MH_PANEL (a local store) or the released store."""
    if os.environ.get("MH_PANEL"):
        return Path(os.environ["MH_PANEL"])
    if RELEASE["url"] is None:
        raise RuntimeError("no released panel yet: build one (nextflow/build.nf) and set MH_PANEL to its folder")
    return fetch_store(RELEASE["url"], RELEASE["sha256"], RELEASE["version"])


class Panel:
    """A panel store on disk. Genotypes are memory-mapped and read one gene block at a time."""

    def __init__(self, path=None):
        self.path = Path(path) if path else default_path()
        self.manifest = json.loads((self.path / "manifest.json").read_text())
        if self.manifest["format"] != FORMAT:
            raise ValueError(f"panel format {self.manifest['format']} is not supported (expected {FORMAT})")
        for name, size in self.manifest["files"].items():
            f = self.path / name
            if not f.exists() or f.stat().st_size != size:
                raise RuntimeError(f"incomplete panel store: {name} is missing or truncated; copy it again")
        self.strains = pd.read_csv(self.path / "strains.tsv", sep="\t", index_col=0)
        self.genes = pd.read_csv(self.path / "genes.tsv", sep="\t", index_col=0)
        self.variants = pd.read_csv(self.path / "variants.tsv.gz", sep="\t", keep_default_na=False,
                                    dtype={"old_rec": str, "aa_change": str})
        g = self.path / "genotypes"
        self._geno = sparse.csr_matrix(
            (np.load(g / "data.npy", mmap_mode="r"), np.load(g / "indices.npy", mmap_mode="r"), np.load(g / "indptr.npy")),
            shape=tuple(json.loads((g / "shape.json").read_text())))

    def gene(self, gene: str) -> tuple[pd.DataFrame, np.ndarray]:
        """Variant rows of one gene and their genotypes: int8 variants x strains (-1 missing, 1 het, 2 hom alt)."""
        if gene not in self.genes.index:
            raise KeyError(f"gene {gene!r} is not a protein-coding R64 gene (use systematic names, e.g. YNL079C)")
        first, last = (int(x) for x in self.genes.loc[gene, ["first", "last"]])
        if first < 0:
            return self.variants.iloc[0:0], np.zeros((0, len(self.strains)), np.int8)
        return self.variants.iloc[first:last + 1], self._geno[first:last + 1].toarray()

    def kinship(self, loco: str | None = None) -> np.ndarray:
        """Genomic relationship matrix; loco='chromosome4' leaves that chromosome out."""
        return np.load(self.path / ("kinship.npy" if loco is None else f"kinship_loco/{loco}.npy"))

    def phenotypes(self) -> pd.DataFrame:
        return pd.read_csv(self.path / "phenotypes.tsv.gz", sep="\t", index_col=0)

    def orfs(self) -> pd.DataFrame:
        return pd.read_csv(self.path / "orfs.tsv", sep="\t")

    def presence_absence(self) -> pd.DataFrame:
        return pd.DataFrame(np.load(self.path / "presence_absence.npy"), index=self.strains.index, columns=self.orfs()["orf"])

    def copy_number(self) -> pd.DataFrame:
        return pd.DataFrame(np.load(self.path / "copy_number.npy"), index=self.strains.index, columns=self.orfs()["orf"])

    def verify(self) -> None:
        """Check every file against SHA256SUMS (slow; run after copying a store by hand)."""
        for line in (self.path / "SHA256SUMS").read_text().splitlines():
            digest, name = line.split("  ", 1)
            if sha256(self.path / name) != digest:
                raise RuntimeError(f"{name}: SHA256 mismatch")
