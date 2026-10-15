# Mutation Hunter 2.0 — Plano de implementação

> **Para agentes:** SUB-SKILL OBRIGATÓRIA: use superpowers:subagent-driven-development (recomendado) ou superpowers:executing-plans para executar este plano tarefa por tarefa. Os passos usam checkbox (`- [ ]`).

**Padrões decididos pelo Eduardo em 15/10, após a calibração (`docs/calibration/README.md`):** modelo misto (EMMAX) com parentesco de SNPs + 10 PCs do genoma acessório (ORFs) como covariáveis, fenótipo com transformação normal por postos (INT), GWAS por gene, raridade ≤ 2%. Os testes fixam `max_freq=0.10` de propósito, para não depender do padrão.

**Objetivo:** transformar o protótipo da IC num pacote Python instalável que faz PheWAS (variante/gene → 35 fenótipos) e GWAS (fenótipo → genes/variantes) sobre o painel dos 1011 genomas de *S. cerevisiae*.

**Arquitetura:** seis módulos pequenos em `src/mutation_hunter/` (dados → variantes → estatística → análise → relatório → CLI). O Nextflow é uma casca que chama a CLI. Os testes usam um painel falso de 40 cepas nos formatos exatos do upstream, servido por `file://`.

**Stack:** Python ≥ 3.11, numpy, pandas, scipy, matplotlib; pytest; uv para o ambiente de desenvolvimento; Nextflow DSL2.

**Spec:** `docs/superpowers/specs/2026-10-15-mutation-hunter-2-design.md`

> [!note] Código já verificado
> Todo o código deste plano foi rodado num rascunho em 15/10/2026: 27 testes passando e execução com os dados reais (`phewas YNL079C`; `gwas YPDCUSO410MM` em gene e variante). Copie os blocos como estão.

## Restrições globais

- Python `>=3.11`. Dependências de execução só: `numpy>=1.26`, `pandas>=2.1`, `scipy>=1.11`, `matplotlib>=3.8`. Testes: `pytest>=8`.
- Nenhum teste acessa a internet. Tudo que baixa usa `MH_BASE_URL` e `MH_CACHE`, que o fixture `panel` redireciona.
- Nenhum arquivo de dados entra no git (cache fica em `~/.cache/mutation_hunter`).
- **Commits sem coautoria do Claude:** nada de `Co-Authored-By`, nem menção a Claude/IA na mensagem. O autor é o `git config` do repositório (Eduardo).
- Não fazer `git push` nem criar tag remota sem pedido explícito do Eduardo.
- Estilo ponytail: o mínimo que funciona, stdlib antes de dependência, sem abstração especulativa. Identificadores e docstrings em inglês.
- Nome provisório `mutation-hunter` / `mutation_hunter`; não espalhar o nome além de `pyproject.toml`, `__init__.py`, `cli.py`, `README.md`, `CITATION.cff`.
- Comandos rodam na raiz do repositório (`~/mutation_hunter`).

## Foco de revisão

1. **IDs com sublinhado** (`SACE_YBA`, 93 cepas): devem sobreviver inteiros na leitura do FASTA e cruzar com fenótipo e distância. → `test_load_gene_keeps_underscored_ids_and_absent_strains` (Tarefa 1).
2. **Cepa sem o gene**: tem de contar como ausente, nunca como referência. → `test_phewas_gene_burden_excludes_absent_strains` (Tarefa 4).
3. **Download corrompido ou dado upstream alterado**: erro claro, sem arquivo parcial no cache. → `test_fetch_rejects_md5_mismatch` (Tarefa 1).
4. **Entrada errada do usuário** (consulta malformada, alelo de referência errado, posição fora da proteína, fenótipo sem cepas em comum, condição inexistente): sair com código 1 e mensagem que diga o que corrigir. → `test_carriers_of_validates` (Tarefa 2) e `test_bad_input_exits_1_with_message` (Tarefa 5).
5. **Execução paralela no Nextflow**: vários processos não podem extrair o proteoma ao mesmo tempo nem disputar threads de BLAS. → `FETCH` roda antes, sozinho; `GWAS` fixa as threads em `task.cpus`; `test_pipeline_runs_and_merges` (Tarefa 6).

## Estrutura de arquivos

| Arquivo | Responsabilidade |
|---|---|
| `pyproject.toml` | metadados, dependências, comando `mutation-hunter` |
| `src/mutation_hunter/__init__.py` | versão |
| `src/mutation_hunter/data.py` | baixar, conferir MD5, cache, ler os 4 arquivos públicos |
| `src/mutation_hunter/variants.py` | consulta `GENE[:REFposALT]`, variantes raras, portadoras |
| `src/mutation_hunter/stats.py` | parentesco, PCs, INT, modelo misto EMMAX, OLS vetorizado, hipergeométrico, BH |
| `src/mutation_hunter/analysis.py` | `phewas()` e `gwas()` |
| `src/mutation_hunter/report.py` | TSV, JSON de proveniência com λGC, HTML de arquivo único |
| `src/mutation_hunter/cli.py` | `fetch`, `phewas`, `gwas`, `fdr` |
| `nextflow/main.nf` | FETCH → GWAS por fenótipo → FDR |
| `tests/conftest.py` | fixture `panel`: painel falso de 40 cepas |
| `tests/test_*.py` | um arquivo por módulo + Nextflow |
| `docs/calibration/` | evidência e scripts da escolha do modelo (não é código do pacote) |
| `README.md`, `LICENSE`, `CITATION.cff`, `.gitignore`, `.github/workflows/tests.yml` | repositório aberto e reprodutível |

---

### Tarefa 1: Preservar a versão da IC, criar o pacote e a camada de dados

**Arquivos:**
- Remover: `Python scripts/`, `hyper_script.R`, `MHpipeline.nf`, `README.md` (antigo)
- Criar: `pyproject.toml`, `.gitignore`, `README.md` (provisório), `src/mutation_hunter/__init__.py`, `src/mutation_hunter/data.py`, `tests/conftest.py`, `tests/test_data.py`

**Interfaces:**
- Produz: `data.BASE_URL: str`, `data.FILES: dict[str, str]` (nome → MD5), `data.PHENOTYPES`, `data.DISTANCE`, `data.ORF_DISTANCE`, `data.PROTEOME` (nomes de arquivo), `data.cache_dir() -> Path`, `data.md5(path) -> str`, `data.fetch(name) -> Path`, `data.proteome_dir() -> Path`, `data.genes() -> list[str]`, `data.load_gene(gene) -> dict[str, str]` (cepa → sequência; levanta `KeyError`), `data.load_phenotypes() -> DataFrame` (cepas × condições), `data.load_distance() -> DataFrame` (cepas × cepas, SNPs), `data.load_orf_distance() -> DataFrame` (cepas × cepas, ORFs presentes em só uma das duas).
- Produz o fixture `panel` (em `conftest.py`): cepas `A00`…`A37`, `SACE_YBA`, `SACE_YCL`; condições `COND1` (sinal plantado), `COND2` (ruído); genes `YAL001C` (variante `A3V` em `A00, A01, A02, SACE_YBA`) e `YAL002W` (variante `T3V` em `A10, A11, A12`; `A37` não tem o gene).

- [ ] **Passo 1: Marcar a versão da IC e remover os scripts antigos**

```bash
git tag -a v0-ic -m "Original undergraduate (IC) scripts, FAPESP 2022/13111-4"
git rm -r "Python scripts" hyper_script.R MHpipeline.nf README.md
```

A tag fica só local. Enviar ao GitHub (`git push origin v0-ic`) apenas quando o Eduardo pedir.

- [ ] **Passo 2: Criar `pyproject.toml`, `.gitignore`, README provisório e `__init__.py`**

`pyproject.toml`:
```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "mutation-hunter"
version = "2.0.0.dev0"
description = "Rare amino-acid variants vs. phenotypes in the 1011 Saccharomyces cerevisiae genomes panel"
readme = "README.md"
license = "MIT"
authors = [{ name = "Eduardo Menoti Silva" }]
requires-python = ">=3.11"
dependencies = ["numpy>=1.26", "pandas>=2.1", "scipy>=1.11", "matplotlib>=3.8"]

[project.optional-dependencies]
test = ["pytest>=8"]

[project.scripts]
mutation-hunter = "mutation_hunter.cli:main"
```

`.gitignore`:
```
.venv/
__pycache__/
*.egg-info/
.pytest_cache/
work/
.nextflow*
results/
*.npz
```

`README.md` (provisório; a versão final vem na Tarefa 7):
```markdown
# mutation-hunter

Rare amino-acid variants vs. phenotypes in the 1011 *S. cerevisiae* genomes panel. Work in progress.
```

`src/mutation_hunter/__init__.py`:
```python
"""Rare amino-acid variants vs. phenotypes in the 1011 S. cerevisiae genomes panel."""

__version__ = "2.0.0.dev0"
```

- [ ] **Passo 3: Escrever o fixture e os testes que falham**

`tests/conftest.py`:
```python
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
```

`tests/test_data.py`:
```python
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
```

- [ ] **Passo 4: Rodar e ver falhar**

Run: `uv run --extra test pytest tests/test_data.py -v`
Esperado: erro de coleta, `ModuleNotFoundError: No module named 'mutation_hunter.data'` (ou `ImportError`).

- [ ] **Passo 5: Implementar `data.py`**

`src/mutation_hunter/data.py`:
```python
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
```

- [ ] **Passo 6: Rodar e ver passar**

Run: `uv run --extra test pytest tests/test_data.py -v`
Esperado: 4 passed.

- [ ] **Passo 7: Commit**

```bash
git add pyproject.toml uv.lock .gitignore README.md src tests
git commit -m "Restart as a Python package: data download, MD5 check and loaders

The original IC scripts are preserved at tag v0-ic."
```

---

### Tarefa 2: Variantes raras e portadoras

**Arquivos:**
- Criar: `src/mutation_hunter/variants.py`, `tests/test_variants.py`

**Interfaces:**
- Consome: nada dos outros módulos (recebe `seqs: dict[str, str]` como `data.load_gene` devolve).
- Produz:
  - `variants.parse_query(query: str) -> tuple[str, tuple[str, int, str] | None]` — levanta `ValueError("bad query ...")`.
  - `variants.rare_variants(gene, seqs, max_freq: float, min_carriers: int) -> DataFrame` com colunas `gene, variant, pos, ref, alt, n_carriers, carriers` (`carriers` é `list[str]`; `variant` é `"A3V"`; `pos` é 1-based; `ref` é o alelo majoritário).
  - `variants.carriers_of(gene, seqs, ref, pos, alt) -> list[str]` — levanta `ValueError` com `"out of range"`, `"major allele is X, not Y"` ou `"no strain carries"`.
  - `variants.carrier_matrix(seqs, carriers: dict[str, list[str]]) -> DataFrame` — linhas = cepas que têm o gene, colunas = testes, 1.0 = portadora.

- [ ] **Passo 1: Escrever os testes que falham**

`tests/test_variants.py`:
```python
import pytest

from mutation_hunter import variants

SEQS = {"S1": "MKA*", "S2": "MKA*", "S3": "MKA*", "S4": "MKV*", "S5": "MXA*"}


def test_parse_query():
    assert variants.parse_query("YNL079C") == ("YNL079C", None)
    assert variants.parse_query("YAR019W-A:A123*") == ("YAR019W-A", ("A", 123, "*"))
    with pytest.raises(ValueError, match="bad query"):
        variants.parse_query("YNL079C:123V")


def test_rare_variants_uses_major_allele_and_ignores_X():
    rv = variants.rare_variants("G", SEQS, max_freq=0.2, min_carriers=1)
    assert rv["variant"].tolist() == ["A3V"]
    assert rv.loc[0, "carriers"] == ["S4"]
    assert variants.rare_variants("G", SEQS, max_freq=0.1, min_carriers=1).empty  # 1/5 = 20% > 10%


def test_carriers_of_validates():
    assert variants.carriers_of("G", SEQS, "A", 3, "V") == ["S4"]
    with pytest.raises(ValueError, match="major allele is A, not V"):
        variants.carriers_of("G", SEQS, "V", 3, "A")
    with pytest.raises(ValueError, match="out of range"):
        variants.carriers_of("G", SEQS, "A", 99, "V")
    with pytest.raises(ValueError, match="no strain carries"):
        variants.carriers_of("G", SEQS, "A", 3, "W")


def test_carrier_matrix_rows_are_strains_with_the_gene():
    G = variants.carrier_matrix(SEQS, {"A3V": ["S4"]})
    assert G["A3V"].to_dict() == {"S1": 0.0, "S2": 0.0, "S3": 0.0, "S4": 1.0, "S5": 0.0}


def test_unaligned_input_is_rejected():
    with pytest.raises(ValueError, match="not aligned"):
        variants.rare_variants("G", {"S1": "MK", "S2": "MKA"}, 0.5, 1)
```

- [ ] **Passo 2: Rodar e ver falhar**

Run: `uv run --extra test pytest tests/test_variants.py -v`
Esperado: `ImportError: cannot import name 'variants'`.

- [ ] **Passo 3: Implementar `variants.py`**

`src/mutation_hunter/variants.py`:
```python
"""Rare amino-acid variants and carrier matrices from the aligned 1011 proteome."""

import re

import numpy as np
import pandas as pd

QUERY = re.compile(r"^(?P<gene>[A-Za-z0-9-]+)(?::(?P<ref>[A-Z*])(?P<pos>\d+)(?P<alt>[A-Z*]))?$")
MISSING = set("X-")


def parse_query(query: str) -> tuple[str, tuple[str, int, str] | None]:
    """'YNL079C' -> ('YNL079C', None); 'YNL079C:A123V' -> ('YNL079C', ('A', 123, 'V'))."""
    m = QUERY.match(query)
    if not m:
        raise ValueError(f"bad query {query!r}; expected GENE or GENE:REFposALT, e.g. YNL079C:A123V")
    if m["pos"] is None:
        return m["gene"], None
    return m["gene"], (m["ref"], int(m["pos"]), m["alt"])


def _matrix(seqs: dict[str, str]) -> np.ndarray:
    if len({len(s) for s in seqs.values()}) != 1:
        raise ValueError("sequences are not aligned (different lengths)")
    return np.array([list(s) for s in seqs.values()])


def _major(col: np.ndarray) -> str:
    values, counts = np.unique(col, return_counts=True)
    return values[counts.argmax()]


def rare_variants(gene: str, seqs: dict[str, str], max_freq: float, min_carriers: int) -> pd.DataFrame:
    """One row per rare allele. ref is the panel's major allele; pos is 1-based.
    Frequency is over strains that have the gene."""
    strains = np.array(list(seqs))
    A = _matrix(seqs)
    rows = []
    for j in np.flatnonzero((A != A[0]).any(axis=0)):
        col = A[:, j]
        ref = _major(col)
        values, counts = np.unique(col, return_counts=True)
        for alt, n in zip(values, counts):
            if alt != ref and alt not in MISSING and min_carriers <= n <= max_freq * len(col):
                rows.append((gene, f"{ref}{j + 1}{alt}", j + 1, ref, alt, int(n), strains[col == alt].tolist()))
    return pd.DataFrame(rows, columns=["gene", "variant", "pos", "ref", "alt", "n_carriers", "carriers"])


def carriers_of(gene: str, seqs: dict[str, str], ref: str, pos: int, alt: str) -> list[str]:
    """Strains carrying alt at pos. ref must be the panel's major allele there."""
    A = _matrix(seqs)
    if not 1 <= pos <= A.shape[1]:
        raise ValueError(f"{gene} has {A.shape[1]} positions; {pos} is out of range")
    col = A[:, pos - 1]
    if (major := _major(col)) != ref:
        raise ValueError(f"{gene} position {pos}: the major allele is {major}, not {ref}")
    found = np.array(list(seqs))[col == alt].tolist()
    if not found:
        raise ValueError(f"no strain carries {gene}:{ref}{pos}{alt}")
    return found


def carrier_matrix(seqs: dict[str, str], carriers: dict[str, list[str]]) -> pd.DataFrame:
    """Strains-with-the-gene x tests, 1.0 = carrier. Strains lacking the gene are not rows (missing, not reference)."""
    return pd.DataFrame({name: pd.Series(1.0, index=c) for name, c in carriers.items()}, index=list(seqs)).fillna(0.0)
```

- [ ] **Passo 4: Rodar e ver passar**

Run: `uv run --extra test pytest tests/test_variants.py -v`
Esperado: 5 passed.

- [ ] **Passo 5: Commit**

```bash
git add src/mutation_hunter/variants.py tests/test_variants.py
git commit -m "Find rare amino-acid variants and their carriers"
```

---

### Tarefa 3: Estatística — modelo misto, teste original e FDR

**Arquivos:**
- Criar: `src/mutation_hunter/stats.py`, `tests/test_stats.py`

**Interfaces:**
- Produz:
  - `stats.kinship(distance: DataFrame) -> DataFrame` — parentesco, diagonal média 1.
  - `stats.pcs(K: DataFrame, k: int) -> DataFrame` — k coordenadas principais do parentesco (colunas `PC1`…).
  - `stats.inverse_normal(y: ndarray) -> ndarray` — transformação normal por postos.
  - `stats.lmm_null(K: ndarray, y: ndarray, X: ndarray) -> ndarray` — `X` inclui o intercepto; devolve `W` (n × n); `W @ M` rotaciona e branqueia.
  - `stats.assoc(G: ndarray, y: ndarray, X: ndarray) -> (beta, se, p)` — `X` **inclui** o intercepto; uma saída por coluna de `G`; NaN quando a coluna é constante.
  - `stats.hypergeom_p(carrier: ndarray[bool], y: ndarray, sd=1.5) -> float` — mesmo resultado do `phyper` do `hyper_script.R`.
  - `stats.bh(p) -> ndarray` — q de Benjamini-Hochberg; NaN fica NaN e não conta como teste.

Por que EMMAX: na grade de calibração (`docs/calibration/README.md`), só o modelo misto ficou calibrado em fenótipos nulos simulados (λ 0,91–1,00); OLS sem correção deu 4,1–10,3 e PCs, 2,1–3,2. Os 10 PCs de ORFs corrigem a estrutura do genoma acessório que o parentesco de SNPs não vê; o INT conserta a cauda quando há poucas portadoras. O modelo misto estima a razão de variâncias uma vez por fenótipo e reduz cada teste a um OLS sobre dados rotacionados, o que permite reaproveitar `assoc`.

- [ ] **Passo 1: Escrever os testes que falham**

`tests/test_stats.py`:
```python
import numpy as np
import pandas as pd
from scipy import stats as sps

from mutation_hunter import stats


def test_assoc_matches_full_ols():
    rng = np.random.default_rng(1)
    n = 60
    C = rng.normal(size=(n, 2))
    g = (rng.random(n) < 0.2).astype(float)
    y = 0.5 * g + C @ [0.3, -0.2] + rng.normal(size=n)
    beta, se, p = stats.assoc(g[:, None], y, np.column_stack([np.ones(n), C]))
    X = np.column_stack([np.ones(n), C, g])
    coef, rss, *_ = np.linalg.lstsq(X, y, rcond=None)
    se_full = np.sqrt(rss[0] / (n - 4) * np.linalg.inv(X.T @ X)[3, 3])
    assert np.isclose(beta[0], coef[3]) and np.isclose(se[0], se_full)
    assert np.isclose(p[0], 2 * sps.t.sf(abs(coef[3] / se_full), n - 4))


def test_assoc_constant_column_is_nan():
    beta, se, p = stats.assoc(np.zeros((10, 1)), np.arange(10.0), np.ones((10, 1)))
    assert np.isnan(p[0])


def test_hypergeom_matches_original_r_script():
    y = np.array([0, 0, 0, 0, 0, 0, 0, 0, 10, 10.0])  # 2 extreme strains
    carrier = np.array([0, 0, 0, 0, 0, 0, 0, 1, 1, 1], bool)
    # R: phyper(k-1, N, M-N, n, lower.tail=FALSE) with M=10, N=2, n=3, k=2
    assert np.isclose(stats.hypergeom_p(carrier, y), sps.hypergeom.sf(1, 10, 2, 3))


def test_bh_skips_nan():
    q = stats.bh([0.01, np.nan, 0.04])
    assert np.isnan(q[1]) and np.allclose(q[[0, 2]], [0.02, 0.04])


def test_lmm_with_identity_kinship_equals_ols():
    rng = np.random.default_rng(2)
    n = 50
    g = (rng.random(n) < 0.3).astype(float)[:, None]
    y = g[:, 0] + rng.normal(size=n)
    ones = np.ones((n, 1))
    W = stats.lmm_null(np.eye(n), y, ones)  # K = I: no structure, W only rescales
    assert np.isclose(stats.assoc(W @ g, W @ y, W @ ones)[2][0], stats.assoc(g, y, ones)[2][0])


def test_kinship_relatives_share_more():
    d = pd.DataFrame([[0, .1, 1, 1], [.1, 0, 1, 1], [1, 1, 0, .1], [1, 1, .1, 0]], dtype=float)
    K = stats.kinship(d).to_numpy()
    assert K[0, 1] > 0 > K[0, 2] and np.isclose(np.diag(K).mean(), 1)


def test_inverse_normal_is_rank_based_and_symmetric():
    z = stats.inverse_normal(np.array([0.1, 5.0, 0.2, 100.0]))
    assert list(np.argsort(z)) == [0, 2, 1, 3] and np.isclose(z.sum(), 0)
```

- [ ] **Passo 2: Rodar e ver falhar**

Run: `uv run --extra test pytest tests/test_stats.py -v`
Esperado: `ImportError: cannot import name 'stats' from 'mutation_hunter'`.

- [ ] **Passo 3: Implementar `stats.py`**

`src/mutation_hunter/stats.py`:
```python
"""Association tests: EMMAX-style linear mixed model, hypergeometric baseline, BH.

Model choice and defaults come from the calibration grid in docs/calibration/README.md."""

import numpy as np
import pandas as pd
from scipy import optimize, stats


def kinship(distance: pd.DataFrame) -> pd.DataFrame:
    """Genetic relatedness from a distance matrix: double-centred -D^2/2 (classical MDS), mean diagonal 1."""
    D2 = distance.to_numpy(float) ** 2
    J = np.eye(len(D2)) - 1 / len(D2)
    B = -0.5 * J @ D2 @ J
    return pd.DataFrame(B / np.diag(B).mean(), index=distance.index, columns=distance.columns)


def pcs(K: pd.DataFrame, k: int) -> pd.DataFrame:
    """Top-k principal coordinates of a kinship matrix."""
    s, U = np.linalg.eigh(K.to_numpy())
    top = np.argsort(s)[::-1][:k]
    return pd.DataFrame(U[:, top] * np.sqrt(np.clip(s[top], 0, None)), index=K.index,
                        columns=[f"PC{i + 1}" for i in range(k)])


def inverse_normal(y: np.ndarray) -> np.ndarray:
    """Rank-based inverse normal transform; keeps the t-test calibrated when carriers are few and y is skewed."""
    return stats.norm.ppf((stats.rankdata(y) - 0.5) / len(y))


def lmm_null(K: np.ndarray, y: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Fit y = X b + u + e, u ~ N(0, sg2*K), e ~ N(0, se2*I) once under the null (EMMAX, Kang 2010).
    Returns W such that W @ M whitens any strains-x-k matrix M; then OLS on W-rotated data is the LMM test."""
    s, U = np.linalg.eigh(K)
    s = np.clip(s, 0, None)
    yt, Xt = U.T @ y, U.T @ X

    def nll(log_delta):  # ML profile likelihood of delta = se2/sg2
        v = s + np.exp(log_delta)
        Xw = Xt / v[:, None]
        r = yt - Xt @ np.linalg.solve(Xt.T @ Xw, Xw.T @ yt)
        return 0.5 * (len(y) * np.log((r * r / v).sum()) + np.log(v).sum())

    delta = np.exp(optimize.minimize_scalar(nll, bounds=(-10, 10), method="bounded").x)
    return U.T / np.sqrt(s + delta)[:, None]


def assoc(G: np.ndarray, y: np.ndarray, X: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """OLS y ~ X + g for each column g of G, all at once (Frisch-Waugh-Lovell).
    X must include the intercept. Returns beta, se, p; NaN where g is constant after adjustment."""
    Q, _ = np.linalg.qr(X)
    yr = y - Q @ (Q.T @ y)
    Gr = G - Q @ (Q.T @ G)
    ss = (Gr**2).sum(axis=0)
    df = len(y) - X.shape[1] - 1
    with np.errstate(divide="ignore", invalid="ignore"):
        ss = np.where(ss > 1e-12, ss, np.nan)
        beta = Gr.T @ yr / ss
        se = np.sqrt(((yr**2).sum() - beta**2 * ss) / df / ss)
        p = 2 * stats.t.sf(np.abs(beta / se), df)
    return beta, se, p


def hypergeom_p(carrier: np.ndarray, y: np.ndarray, sd: float = 1.5) -> float:
    """Original Mutation Hunter test: are carriers enriched among strains with y > mean + sd*SD?"""
    extreme = y > y.mean() + sd * y.std(ddof=1)
    return float(stats.hypergeom.sf((carrier & extreme).sum() - 1, len(y), extreme.sum(), carrier.sum()))


def bh(p) -> np.ndarray:
    """Benjamini-Hochberg q-values; NaN p-values stay NaN and do not count as tests."""
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.any():
        q[ok] = stats.false_discovery_control(p[ok])
    return q
```

- [ ] **Passo 4: Rodar e ver passar**

Run: `uv run --extra test pytest tests/test_stats.py -v`
Esperado: 7 passed.

- [ ] **Passo 5: Commit**

```bash
git add src/mutation_hunter/stats.py tests/test_stats.py
git commit -m "Add mixed-model association test, hypergeometric baseline and BH"
```

---

### Tarefa 4: PheWAS e GWAS

**Arquivos:**
- Criar: `src/mutation_hunter/analysis.py`, `tests/test_analysis.py`

**Interfaces:**
- Consome: `data.load_gene`, `data.genes`, `data.load_phenotypes`, `data.load_distance`, `data.load_orf_distance` (Tarefa 1); `variants.parse_query`, `variants.rare_variants`, `variants.carriers_of`, `variants.carrier_matrix` (Tarefa 2); `stats.kinship`, `stats.pcs`, `stats.inverse_normal`, `stats.lmm_null`, `stats.assoc`, `stats.hypergeom_p`, `stats.bh` (Tarefa 3).
- Produz:
  - `analysis.ORF_PCS = 10`
  - `analysis.COLUMNS = ["phenotype", "gene", "test", "n", "n_carriers", "beta", "se", "p", "p_hypergeom", "q"]`
  - `analysis.phewas(query: str, max_freq=0.02, min_carriers=3) -> (results: DataFrame, rare: DataFrame)` — `results` com `COLUMNS`, uma linha por fenótipo; `rare` é a saída de `rare_variants` do gene.
  - `analysis.gwas(phenotype: Series, level="gene" | "variant", max_freq=0.02, min_carriers=3) -> DataFrame` com `COLUMNS`. `phenotype.name` vira a coluna `phenotype`.

- [ ] **Passo 1: Escrever os testes que falham**

`tests/test_analysis.py`:
```python
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
```

- [ ] **Passo 2: Rodar e ver falhar**

Run: `uv run --extra test pytest tests/test_analysis.py -v`
Esperado: `ImportError: cannot import name 'analysis'`.

- [ ] **Passo 3: Implementar `analysis.py`**

`src/mutation_hunter/analysis.py`:
```python
"""PheWAS (one variant or gene vs. all phenotypes) and GWAS (one phenotype vs. all genes)."""

import numpy as np
import pandas as pd

from . import data, stats, variants

COLUMNS = ["phenotype", "gene", "test", "n", "n_carriers", "beta", "se", "p", "p_hypergeom", "q"]
ORF_PCS = 10  # accessory-genome structure that SNP kinship misses (calibration grid)


def _structure() -> tuple[pd.DataFrame, pd.DataFrame]:
    """SNP kinship (random effect) and top ORF presence/absence coordinates (fixed covariates)."""
    return stats.kinship(data.load_distance()), stats.pcs(stats.kinship(data.load_orf_distance()), ORF_PCS)


def _null(y: pd.Series, K: pd.DataFrame, C: pd.DataFrame):
    """Analysed strains, raw and rank-normalised phenotype, design X and the LMM whitening matrix."""
    strains = y.dropna().index.intersection(K.index).intersection(C.index)
    raw = y[strains].to_numpy(float)
    yv = stats.inverse_normal(raw)
    X = np.column_stack([np.ones(len(strains)), C.loc[strains].to_numpy()])
    return strains, raw, yv, X, stats.lmm_null(K.loc[strains, strains].to_numpy(), yv, X)


def _test(G: pd.DataFrame, null, min_carriers: int) -> pd.DataFrame:
    """Test each column of G (1 = carrier) against the phenotype in null."""
    strains, raw, yv, X, W = null
    g = G.reindex(strains)  # strains lacking the gene -> NaN: missing, not reference
    present = g.notna().to_numpy()
    n, n_car = present.sum(axis=0), g.sum().to_numpy()
    gi = g.fillna(g.mean()).to_numpy(float)  # mean imputation: missing strains add no signal
    beta, se, p = stats.assoc(W @ gi, W @ yv, W @ X)
    ok = (n_car >= min_carriers) & (n_car < n)
    p_hg = [stats.hypergeom_p(gi[present[:, i], i] > 0, raw[present[:, i]]) if ok[i] else np.nan
            for i in range(g.shape[1])]
    return pd.DataFrame({
        "test": G.columns, "n": n, "n_carriers": n_car.astype(int),
        "beta": beta, "se": se, "p": np.where(ok, p, np.nan), "p_hypergeom": p_hg,
    })


def phewas(query: str, max_freq=0.02, min_carriers=3) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Test a variant ('YNL079C:A123V') or a gene burden ('YNL079C') against every phenotype.
    Returns (results, the rare variants of the gene)."""
    gene, var = variants.parse_query(query)
    seqs = data.load_gene(gene)
    rv = variants.rare_variants(gene, seqs, max_freq, min_carriers)
    if var:
        carriers = variants.carriers_of(gene, seqs, *var)
    else:
        carriers = sorted({s for c in rv["carriers"] for s in c})  # burden: any rare allele (CAST)
        if not carriers:
            raise ValueError(f"{gene} has no rare variant at max_freq={max_freq}, min_carriers={min_carriers}")
    G = variants.carrier_matrix(seqs, {query: carriers})
    Y, (K, C) = data.load_phenotypes(), _structure()
    res = pd.concat([_test(G, _null(Y[ph], K, C), min_carriers).assign(phenotype=ph) for ph in Y], ignore_index=True)
    res["gene"], res["q"] = gene, stats.bh(res["p"])
    return res[COLUMNS], rv


def gwas(phenotype: pd.Series, level="gene", max_freq=0.02, min_carriers=3) -> pd.DataFrame:
    """Test every gene burden (level='gene') or every rare variant (level='variant') against one phenotype."""
    null, parts = _null(phenotype, *_structure()), []
    for gene in data.genes():
        seqs = data.load_gene(gene)
        rv = variants.rare_variants(gene, seqs, max_freq, min_carriers)
        if rv.empty:
            continue
        G = variants.carrier_matrix(seqs, dict(zip(rv["variant"], rv["carriers"])))
        if level == "gene":
            G = G.max(axis=1).to_frame(gene)
        parts.append(_test(G, null, min_carriers).assign(gene=gene))
    res = pd.concat(parts, ignore_index=True).assign(phenotype=phenotype.name)
    res["q"] = stats.bh(res["p"])
    return res[COLUMNS]
```

- [ ] **Passo 4: Rodar e ver passar**

Run: `uv run --extra test pytest tests/test_analysis.py -v`
Esperado: 5 passed.

- [ ] **Passo 5: Commit**

```bash
git add src/mutation_hunter/analysis.py tests/test_analysis.py
git commit -m "Add PheWAS (variant or gene vs. all conditions) and GWAS (condition vs. all genes)"
```

---

### Tarefa 5: Relatório e linha de comando

**Arquivos:**
- Criar: `src/mutation_hunter/report.py`, `src/mutation_hunter/cli.py`, `tests/test_cli.py`
- Modificar: nenhum (o comando `mutation-hunter` já está declarado no `pyproject.toml`)

**Interfaces:**
- Consome: `analysis.phewas`, `analysis.gwas` (Tarefa 4); `stats.bh` (Tarefa 3); `data.FILES`, `data.BASE_URL`, `data.fetch`, `data.proteome_dir`, `data.load_phenotypes`, `data.load_distance` (Tarefa 1).
- Produz:
  - `report.lambda_gc(p) -> float`, `report.provenance(params: dict) -> dict`, `report.write(prefix, res, params, rare=None)` — grava `<prefix>.tsv`, `.json`, `.html`.
  - `cli.main(argv: list[str] | None = None) -> int` — 0 = sucesso, 1 = erro de entrada (mensagem em stderr com prefixo `mutation-hunter: error:`).
  - Subcomandos: `fetch` (imprime as condições, uma por linha), `phewas QUERY -o PREFIX`, `gwas PHENOTYPE [--level gene|variant] -o PREFIX`, `fdr TSV... -o PREFIX`; `phewas`/`gwas` aceitam `--max-freq` e `--min-carriers`.

O teste de `report` é feito pela CLI: é por ela que o usuário vê os arquivos.

- [ ] **Passo 1: Escrever os testes que falham**

`tests/test_cli.py`:
```python
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
```

- [ ] **Passo 2: Rodar e ver falhar**

Run: `uv run --extra test pytest tests/test_cli.py -v`
Esperado: `ImportError: cannot import name 'cli'`.

- [ ] **Passo 3: Implementar `report.py`**

`src/mutation_hunter/report.py`:
```python
"""Write results as TSV (the product), JSON (provenance) and a single-file HTML report."""

import base64
import html
import io
import json
import sys
from datetime import datetime, timezone

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy  # noqa: E402
from scipy import stats  # noqa: E402

from . import __version__, data  # noqa: E402


def lambda_gc(p) -> float:
    """Genomic inflation: median observed chi2 / expected (0.455). ~1 is calibrated; >>1 means confounding."""
    p = np.asarray(p, float)
    p = p[~np.isnan(p)]
    return float(np.median(stats.chi2.isf(p, 1)) / stats.chi2.ppf(0.5, 1)) if p.size else float("nan")


def provenance(params: dict) -> dict:
    return {
        "mutation_hunter": __version__,
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "params": params,
        "data_source": data.BASE_URL,
        "data_md5": data.FILES,
        "python": sys.version.split()[0],
        "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__,
    }


def _png(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return f'<img src="data:image/png;base64,{base64.b64encode(buf.getvalue()).decode()}">'


def _forest(res: pd.DataFrame) -> str:
    r = res.dropna(subset=["beta"]).sort_values("beta")
    fig, ax = plt.subplots(figsize=(6, 0.25 * len(r) + 1))
    ax.errorbar(r["beta"], r["phenotype"], xerr=1.96 * r["se"], fmt="o")
    ax.axvline(0, color="grey", lw=0.8)
    ax.set_xlabel("effect on growth ratio (beta ± 95% CI)")
    return _png(fig)


def _manhattan(res: pd.DataFrame) -> str:
    r = res.dropna(subset=["p"]).sort_values("gene")
    fig, ax = plt.subplots(figsize=(9, 3))
    ax.scatter(range(len(r)), -np.log10(r["p"]), s=4)
    ax.set_xlabel("tests, ordered by systematic gene name")
    ax.set_ylabel("-log10 p")
    return _png(fig)


def _lollipop(rv: pd.DataFrame) -> str:
    fig, ax = plt.subplots(figsize=(9, 2.5))
    ax.stem(rv["pos"], rv["n_carriers"])
    ax.set_xlabel("protein position")
    ax.set_ylabel("carriers")
    return _png(fig)


def write(prefix: str, res: pd.DataFrame, params: dict, rare: pd.DataFrame | None = None) -> None:
    """<prefix>.tsv, <prefix>.json and <prefix>.html. rare is given for PheWAS, None for GWAS."""
    res.to_csv(f"{prefix}.tsv", sep="\t", index=False)
    prov = provenance(params) | {"lambda_gc": lambda_gc(res["p"])}
    with open(f"{prefix}.json", "w") as f:
        json.dump(prov, f, indent=2)
    figs = [_forest(res)] if rare is not None else [_manhattan(res)]
    if rare is not None and not rare.empty:
        figs.append(_lollipop(rare))
    top = res.nsmallest(50, "p").to_html(index=False, float_format="%.3g")
    title = html.escape(f"mutation-hunter {params.get('command', '')}: {params.get('query', params.get('phenotype', ''))}")
    with open(f"{prefix}.html", "w") as f:
        f.write(f"<!doctype html><meta charset='utf-8'><title>{title}</title><h1>{title}</h1>"
                f"<p>Genomic inflation &lambda;<sub>GC</sub> = {prov['lambda_gc']:.2f} (~1 is calibrated)</p>"
                f"{''.join(figs)}<h2>Top tests</h2>{top}"
                f"<h2>Provenance</h2><pre>{html.escape(json.dumps(prov, indent=2))}</pre>")
```

- [ ] **Passo 4: Implementar `cli.py`**

`src/mutation_hunter/cli.py`:
```python
"""Command line: mutation-hunter fetch | phewas | gwas | fdr."""

import argparse
import sys
from pathlib import Path

import pandas as pd

from . import __version__, analysis, data, report, stats


def _phenotype(arg: str) -> pd.Series:
    """A condition name from the panel, or a TSV file: strain<TAB>value with a header line."""
    if Path(arg).is_file():
        y = pd.read_csv(arg, sep="\t", index_col=0).iloc[:, 0].rename(Path(arg).stem)
        if y.index.intersection(data.load_distance().index).empty:
            raise ValueError(f"{arg}: no strain ID in common with the 1011 panel (expected IDs like AAA, SACE_YBA)")
        return y
    Y = data.load_phenotypes()
    if arg not in Y:
        raise ValueError(f"unknown phenotype {arg!r}; run 'mutation-hunter fetch' to list the 35 conditions")
    return Y[arg]


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="mutation-hunter", description=__doc__)
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("fetch", help="download and verify the public data; print the 35 conditions")
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("-o", "--out", required=True, help="output prefix (writes .tsv, .json, .html)")
    common.add_argument("--max-freq", type=float, default=0.02, help="rare = allele frequency <= this (default 0.02)")
    common.add_argument("--min-carriers", type=int, default=3, help="minimum carriers to test (default 3)")
    ph = sub.add_parser("phewas", parents=[common], help="one variant or gene vs. all 35 phenotypes")
    ph.add_argument("query", help="GENE (burden of rare variants) or GENE:REFposALT, e.g. YNL079C:A123V")
    gw = sub.add_parser("gwas", parents=[common], help="one phenotype vs. all genes")
    gw.add_argument("phenotype", help="one of the 35 conditions, or a TSV file strain<TAB>value")
    gw.add_argument("--level", choices=["gene", "variant"], default="gene")
    fd = sub.add_parser("fdr", help="merge result TSVs and recompute q over all tests")
    fd.add_argument("tsv", nargs="+")
    fd.add_argument("-o", "--out", required=True, help="output prefix")
    a = p.parse_args(argv)

    try:
        if a.command == "fetch":
            for name in data.FILES:
                data.fetch(name)
            data.proteome_dir()
            print("\n".join(data.load_phenotypes().columns))
        elif a.command == "phewas":
            res, rare = analysis.phewas(a.query, a.max_freq, a.min_carriers)
            report.write(a.out, res, vars(a), rare)
        elif a.command == "gwas":
            res = analysis.gwas(_phenotype(a.phenotype), a.level, a.max_freq, a.min_carriers)
            report.write(a.out, res, vars(a))
        else:
            res = pd.concat([pd.read_csv(f, sep="\t") for f in a.tsv], ignore_index=True)
            res["q"] = stats.bh(res["p"])
            res.to_csv(f"{a.out}.tsv", sep="\t", index=False)
    except (ValueError, KeyError, RuntimeError) as e:
        print(f"mutation-hunter: error: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Passo 5: Rodar a suíte inteira**

Run: `uv run --extra test pytest -v`
Esperado: 26 passed.

- [ ] **Passo 6: Conferir o comando instalado**

Run: `uv run mutation-hunter --help`
Esperado: uso com os subcomandos `{fetch,phewas,gwas,fdr}`.

- [ ] **Passo 7: Commit**

```bash
git add src/mutation_hunter/report.py src/mutation_hunter/cli.py tests/test_cli.py
git commit -m "Add CLI and TSV/JSON/HTML reports with provenance and genomic inflation"
```

---

### Tarefa 6: Casca Nextflow

**Arquivos:**
- Criar: `nextflow/main.nf`, `tests/test_nextflow.py`

**Interfaces:**
- Consome: o comando `mutation-hunter` no `PATH` (`fetch`, `gwas`, `fdr`, Tarefa 5).
- Produz: `nextflow run nextflow/main.nf [--phenotypes all|A,B] [--level gene|variant] [--max_freq 0.02] [--min_carriers 3] [--outdir results]` → `results/<condição>.{tsv,json,html}` e `results/gwas_all.tsv` (q recalculado sobre todos os testes).

- [ ] **Passo 1: Escrever o teste que falha**

`tests/test_nextflow.py`:
```python
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from mutation_hunter import data

MAIN = Path(__file__).parents[1] / "nextflow" / "main.nf"


@pytest.mark.skipif(not shutil.which("nextflow"), reason="nextflow not installed")
def test_pipeline_runs_and_merges(panel, tmp_path):
    for name in data.FILES:  # warm the cache here: the subprocesses do not see the patched MD5s
        data.fetch(name)
    env = os.environ | {"PATH": f"{Path(sys.executable).parent}{os.pathsep}{os.environ['PATH']}"}
    subprocess.run(["nextflow", "run", str(MAIN), "--outdir", str(tmp_path / "res"), "--max_freq", "0.1"],
                   cwd=tmp_path, env=env, check=True, capture_output=True)
    merged = pd.read_csv(tmp_path / "res" / "gwas_all.tsv", sep="\t")
    assert set(merged["phenotype"]) == {"COND1", "COND2"} and merged["q"].notna().any()
```

- [ ] **Passo 2: Rodar e ver falhar**

Run: `uv run --extra test pytest tests/test_nextflow.py -v`
Esperado: FAIL — `nextflow` sai com erro porque `nextflow/main.nf` não existe (`CalledProcessError`). Se o Nextflow não estiver instalado, o teste é pulado: instale antes (`nextflow -version` deve responder).

- [ ] **Passo 3: Implementar `main.nf`**

`nextflow/main.nf`:
```groovy
// Thin wrapper: all logic lives in the mutation-hunter package. Runs one GWAS per phenotype in parallel,
// then recomputes q over every test. Requires `mutation-hunter` on PATH.
nextflow.enable.dsl = 2

params.phenotypes = 'all'      // 'all' or comma-separated condition names
params.level = 'gene'          // 'gene' or 'variant'
params.max_freq = 0.02
params.min_carriers = 3
params.outdir = 'results'

process FETCH {
    output:
    stdout

    script:
    'mutation-hunter fetch'
}

process GWAS {
    tag "$pheno"
    publishDir params.outdir, mode: 'copy'

    input:
    val pheno

    output:
    path "${pheno}.tsv", emit: tsv
    path "${pheno}.{json,html}"

    script:  // one BLAS thread per task: parallelism comes from running phenotypes side by side
    """
    export OMP_NUM_THREADS=$task.cpus OPENBLAS_NUM_THREADS=$task.cpus MKL_NUM_THREADS=$task.cpus
    mutation-hunter gwas $pheno --level $params.level --max-freq $params.max_freq \\
        --min-carriers $params.min_carriers -o $pheno
    """
}

process FDR {
    publishDir params.outdir, mode: 'copy'

    input:
    path tsvs

    output:
    path 'gwas_all.tsv'

    script:
    "mutation-hunter fdr $tsvs -o gwas_all"
}

workflow {
    wanted = params.phenotypes == 'all' ? null : params.phenotypes.tokenize(',')
    phenos = FETCH().splitText().map { it.trim() }.filter { it && (!wanted || it in wanted) }
    GWAS(phenos)
    FDR(GWAS.out.tsv.collect())
}
```

- [ ] **Passo 4: Rodar e ver passar**

Run: `uv run --extra test pytest tests/test_nextflow.py -v`
Esperado: 1 passed.

- [ ] **Passo 5: Commit**

```bash
git add nextflow/main.nf tests/test_nextflow.py
git commit -m "Add thin Nextflow wrapper: one GWAS per condition, global FDR"
```

---

### Tarefa 7: Repositório aberto — README, licença, citação e CI

**Arquivos:**
- Criar: `LICENSE`, `CITATION.cff`, `.github/workflows/tests.yml`
- Já no disco, só commitar: `docs/calibration/` (README com as tabelas, 5 scripts, 4 TSVs), `docs/superpowers/`
- Modificar: `README.md` (substituir o provisório)

**Interfaces:** nenhuma de código. A CI roda a mesma suíte (o teste de Nextflow é pulado lá, porque a imagem não tem Nextflow).

> [!warning] Decisões do Eduardo
> A licença MIT e o README em inglês são propostas da spec. Confirme antes do commit. Se a licença mudar, troque `LICENSE`, `license` no `pyproject.toml` e no `CITATION.cff`.

- [ ] **Passo 1: Escrever `LICENSE`**

```
MIT License

Copyright (c) 2026 Eduardo Menoti Silva

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

- [ ] **Passo 2: Escrever `CITATION.cff`**

```yaml
cff-version: 1.2.0
message: "If you use this software, please cite it, and cite Peter et al. 2018 for the data."
title: mutation-hunter
version: 2.0.0.dev0
license: MIT
repository-code: "https://github.com/edumenotti/IC_scripts"
authors:
  - family-names: Menoti Silva
    given-names: Eduardo
references:
  - type: article
    title: "Genome evolution across 1,011 Saccharomyces cerevisiae isolates"
    authors:
      - family-names: Peter
        given-names: Jackson
    journal: Nature
    year: 2018
    doi: 10.1038/s41586-018-0030-5
```

- [ ] **Passo 3: Escrever `.github/workflows/tests.yml`**

```yaml
name: tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python: ["3.11", "3.13"]
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python }}
      - run: pip install -e ".[test]"
      - run: pytest -v
```

- [ ] **Passo 4: Substituir `README.md`**

````markdown
# mutation-hunter

Rare amino-acid variants vs. phenotypes in the 1011 *Saccharomyces cerevisiae* genomes panel
([Peter et al. 2018](https://doi.org/10.1038/s41586-018-0030-5)).

Two questions, one engine:

- **PheWAS** — *I have a variant (or a gene); which of the 35 growth conditions does it affect?*
- **GWAS** — *I have a phenotype; which genes or rare variants explain it?*

## Install

```bash
pip install git+https://github.com/edumenotti/IC_scripts
mutation-hunter fetch        # downloads ~65 MB once, verifies MD5, prints the 35 conditions
```

## Use

```bash
mutation-hunter phewas YNL079C:A123V -o tpm1_a123v   # one variant vs. 35 conditions
mutation-hunter phewas YNL079C -o tpm1               # burden of all rare variants in the gene
mutation-hunter gwas YPDCUSO410MM -o copper          # one condition vs. every gene
mutation-hunter gwas my_phenotype.tsv --level variant -o mine   # your own phenotype: strain<TAB>value
```

Python API: `mutation_hunter.analysis.phewas(...)` and `.gwas(...)` return pandas DataFrames.

All 35 conditions in parallel, with one FDR over everything:

```bash
nextflow run nextflow/main.nf --outdir results   # needs mutation-hunter on PATH
```

## Output

Each run writes `<out>.tsv`, `<out>.json` (parameters, versions, data MD5, genomic inflation λGC) and a
self-contained `<out>.html`.

| column | meaning |
|---|---|
| `n`, `n_carriers` | strains tested / strains carrying the variant (or any rare variant of the gene) |
| `beta`, `se`, `p` | carrier effect on the rank-normalised phenotype (SD units), linear mixed model |
| `p_hypergeom` | the original, unadjusted hypergeometric test (extreme = mean + 1.5 SD), for comparison |
| `q` | Benjamini–Hochberg over all tests in the run |

## Method

- Reference allele = the panel's major allele; positions are 1-based in the aligned proteome.
- Rare = frequency ≤ `--max-freq` (0.02, chosen by calibration; see docs) among strains that have the gene, with ≥ `--min-carriers` (3).
- Gene test = burden (carrier of any rare variant).
- Population structure: EMMAX-style mixed model with SNP kinship, plus the top 10 coordinates of the
  accessory genome (ORF presence/absence) as covariates. Phenotypes are rank-normalised (INT).
  Chosen by calibration: see [docs/calibration](docs/calibration/README.md).
  Check `lambda_gc` in the JSON; values well above 1 mean residual confounding.
- Strains lacking a gene are missing, not reference.

Limitations of the public proteome: SNPs only (no indels); at heterozygous sites one allele was drawn at random.

## Development

```bash
uv run --extra test pytest -v
```

## Acknowledgements

Started as an undergraduate project funded by FAPESP (2022/13111-4). The original scripts are kept at tag `v0-ic`.

## License

MIT
````

- [ ] **Passo 5: Conferir que o README não promete o que o código não faz**

Run: `uv run mutation-hunter phewas --help && uv run mutation-hunter gwas --help && uv run --extra test pytest -q`
Esperado: as opções citadas no README (`--level`, `--max-freq`, `--min-carriers`, `-o`) aparecem; 27 passed (ou 26 passed, 1 skipped sem Nextflow).

- [ ] **Passo 6: Commit**

```bash
git add LICENSE CITATION.cff .github/workflows/tests.yml README.md docs/calibration docs/superpowers
git commit -m "Add README, MIT license, citation file, CI and calibration evidence"
```

---

### Tarefa 8: Verificação com os dados reais

Não é teste automático: confere se o pacote reproduz os números da spec com os dados verdadeiros. Nada desta tarefa entra no git.

**Arquivos:** nenhum no repositório. Saídas em `/tmp/mh-check/`.

- [ ] **Passo 1: Baixar os dados**

Run: `uv run mutation-hunter fetch | wc -l`
Esperado: `35` (cerca de 66 MB baixados uma vez; MD5 conferido).

- [ ] **Passo 2: PheWAS de um gene**

Run: `mkdir -p /tmp/mh-check && uv run mutation-hunter phewas YNL079C -o /tmp/mh-check/tpm1 && head -3 /tmp/mh-check/tpm1.tsv`
Esperado: 35 linhas de resultado; `n` = 970; nenhum q < 0,05; `lambda_gc` ≈ 0,88; arquivos `.json` e `.html` criados. Abrir o HTML no navegador e conferir o forest plot e o lollipop.

- [ ] **Passo 3: GWAS de cobre e inflação**

Run: `uv run mutation-hunter gwas YPDCUSO410MM -o /tmp/mh-check/cu && python3 -c "import json; print(json.load(open('/tmp/mh-check/cu.json'))['lambda_gc'])"`
Esperado: 6.388 testes com p (limiar padrão de 2%); `lambda_gc` ≈ 1,41; 2 genes com q < 0,05 (YOL003C e YPL028W). Leva ~5 min.

- [ ] **Passo 4: Reportar ao Eduardo**

Mandar os números (`n`, testes, `lambda_gc`, genes com q < 0,05) e o link do HTML. Divergência dos números acima indica mudança no código ou nos dados: investigar antes de seguir.

---

## Fora deste plano

Relatório de cepa nova, ESM/SIFT, CNVs e genes acessórios, frequência por clado, release pré-calculada no Zenodo, atlas, bioconda e contêiner. Cada um vira spec própria quando o Eduardo decidir.
