"""Download, verify and load the public 1011 Yeast Genomes files (Peter et al. 2018)."""

import hashlib
import os
import shutil
import tarfile
import urllib.request
from pathlib import Path

import pandas as pd

BASE_URL = "http://1002genomes.u-strasbg.fr/files/"
# Copied from BASE_URL + "md5.txt" on 2026-10-15. A mismatch means upstream changed.
FILES = {
    "phenoMatrix_35ConditionsNormalizedByYPD.tab.gz": "479862335552a050e46de163d1ee5f72",
    "1011DistanceMatrixBasedOnSNPs.tab.gz": "817f7456393dee33f6b41d3d59ddc985",
    "1011DistanceMatrixBasedOnORFs.tab.gz": "7fd9048da9b828041a4333c0e991e6d8",
    "1011proteome.tar.gz": "1bb03292cc91a87fb724d70a33ea3c67",
}
PHENOTYPES, DISTANCE, ORF_DISTANCE, PROTEOME = FILES


def cache_dir() -> Path:
    default = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "mutation_hunter"
    d = Path(os.environ.get("MH_CACHE", default))
    d.mkdir(parents=True, exist_ok=True)
    return d


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(name: str) -> Path:
    """Return the cached path of a public file, downloading and verifying it first if needed."""
    path = cache_dir() / name
    if path.exists():
        return path
    tmp = path.with_name(name + ".part")
    url = os.environ.get("MH_BASE_URL", BASE_URL) + name
    with urllib.request.urlopen(url) as r, open(tmp, "wb") as f:
        shutil.copyfileobj(r, f)
    if md5(tmp) != FILES[name]:
        tmp.unlink()
        raise RuntimeError(f"{name}: MD5 mismatch; the download is corrupt or upstream data changed")
    tmp.rename(path)
    return path


def proteome_dir() -> Path:
    """Directory with one aligned FASTA per gene, extracted once into the cache."""
    d = cache_dir() / "Proteome_1011"
    if not d.exists():
        tmp = cache_dir() / "proteome.part"
        shutil.rmtree(tmp, ignore_errors=True)
        with tarfile.open(fetch(PROTEOME)) as t:
            t.extractall(tmp, filter="data")
        (tmp / "Proteome_1011").rename(d)
        shutil.rmtree(tmp)
    return d


def genes() -> list[str]:
    return sorted(p.stem for p in proteome_dir().glob("*.fasta"))


def load_gene(gene: str) -> dict[str, str]:
    """Strain -> protein sequence for one systematic gene name. Strains lacking the gene are absent."""
    path = proteome_dir() / f"{gene}.fasta"
    if not path.exists():
        raise KeyError(f"gene {gene!r} is not in the 1011 proteome (use systematic names, e.g. YNL079C)")
    seqs, strain = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            # ">SACE_YBA_YNL079C_TPM1_1 ..." -> "SACE_YBA"; IDs may contain "_".
            strain = line[1:].split(f"_{gene}_")[0]
            seqs[strain] = ""
        else:
            seqs[strain] += line.strip()
    return seqs


def load_phenotypes() -> pd.DataFrame:
    """Strains x 35 conditions: growth ratio relative to YPD."""
    return pd.read_csv(fetch(PHENOTYPES), sep="\t", index_col=0)


def load_distance() -> pd.DataFrame:
    """Strains x strains: % non-identical SNP bases."""
    return pd.read_csv(fetch(DISTANCE), sep="\t", index_col=0)


def load_orf_distance() -> pd.DataFrame:
    """Strains x strains: number of ORFs present in only one of the two (accessory genome)."""
    return pd.read_csv(fetch(ORF_DISTANCE), sep="\t", index_col=0)
