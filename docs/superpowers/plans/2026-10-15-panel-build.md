# Subprojeto 1 — Construção do painel: plano de implementação

> **Para agentes:** SUB-SKILL OBRIGATÓRIA: use superpowers:subagent-driven-development (recomendado) ou superpowers:executing-plans para executar este plano tarefa por tarefa. Os passos usam checkbox (`- [ ]`).

**Objetivo:** construir, a partir dos dados públicos do painel 1011, um pacote de dados versionado (~140 MB comprimido) com variantes anotadas, genótipos, parentesco, tabelas de cepas e de ORFs, e um leitor que o usa gene por gene.

**Arquitetura:** `bcftools` normaliza e anota cada cromossomo (`nextflow/bin/annotate.sh`); `mutation_hunter.build` (Python) extrai genótipos com `cyvcf2`, aplica os filtros de qualidade, monta as tabelas e grava o pacote; `nextflow/build.nf` orquestra 16 cromossomos em paralelo; `mutation_hunter.store` lê o pacote com `mmap` e baixa versões publicadas com cache atômico.

**Stack:** Python ≥ 3.11 (numpy, pandas, scipy; extra `build`: cyvcf2, xlrd), bcftools 1.20, Nextflow DSL2, pytest.

**Spec:** `docs/superpowers/specs/2026-10-15-panel-build-design.md`

> [!note] Código já verificado
> Todo o código deste plano rodou num rascunho em 15/10/2026: 60 testes passando (28 do baseline + 32 novos) e a construção real completa (16 cromossomos; 3 min 37 s à mão e 3 min 55 s via Nextflow, numa máquina de 16 núcleos, com resultados idênticos byte a byte). Copie os blocos como estão.

## Restrições globais

- Partir do branch `baseline/proteome` (commit `f5898cb`); trabalhar no branch novo `feat/panel-build`.
- Python `>=3.11`. Dependências de execução não mudam (numpy, pandas, scipy, matplotlib). Novo extra opcional `build = ["cyvcf2>=0.31", "xlrd>=2"]`, só para mantenedores.
- `bcftools >= 1.20` (testado com 1.20). Sem bcftools ou cyvcf2, os testes que dependem deles são pulados.
- Referência fixada: Ensembl R64-1-1 **release 115**. Todas as entradas têm MD5 fixado em `build.INPUTS`.
- Filtros de genótipo (decisão do Eduardo, 15/10): GQ < 20 ou DP < 5 → ausente; heterozigoto em cepa haploide → ausente.
- Codificação no pacote: −1 ausente, 0 referência, 1 het, 2 homozigoto alternativo. Não é dose aditiva.
- Nenhum dado entra no git: `work/`, `panel*/` e pastas de saída do Nextflow ficam fora.
- **Commits sem coautoria do Claude** (sem `Co-Authored-By`, sem menção a IA). Sem `git push`. **Nada vai para o Zenodo** antes da autorização do grupo Schacherer/Liti.
- Estilo ponytail: mínimo que funciona, sem abstração especulativa; identificadores e docstrings em inglês.
- Comandos rodam na raiz do repositório. `BCFTOOLS` = qualquer `bcftools` ≥ 1.20 no `PATH` (na máquina do Eduardo: `~/ic/vcf/bcftools-1.20/bcftools`; ou `conda env create -f environment.yml`).

## Foco de revisão

1. **Entrada de origem alterada ou corrompida, ou cepas que não batem entre arquivos** (o site é HTTP, sem licença declarada): a construção deve parar dizendo o arquivo e o motivo. → `test_check_inputs_rejects_changed_or_missing_files` (Tarefa 1) e `test_assemble_rejects_strains_missing_from_matrices` (Tarefa 4).
2. **Processos paralelos criando o índice da referência ao mesmo tempo** (falhou de verdade no primeiro teste com Nextflow: `EOF marker is absent`): o `.fai` é escrito uma vez pela etapa de referência. → `test_reference_renames_and_drops_mito` confere o `.fai` byte a byte (Tarefa 1).
3. **Cepa sem chamada, com GQ/DP baixo ou "het" em haploide** não pode virar referência nem portadora. → `test_encode_masks_low_quality_and_haploid_hets` e `test_quality_rules_set_missing` (Tarefa 2).
4. **Pacote parcial ou adulterado** (construção interrompida, cópia manual, download corrompido): nunca deve ser lido como válido. → `test_assemble_store_layout` (sem pasta `.part` restante), `test_panel_detects_tampering`, `test_fetch_store_verifies_and_caches` (Tarefas 4 e 5).
5. **Pacote de formato futuro lido por código antigo** → erro explícito. → `test_panel_rejects_unknown_format` (Tarefa 5).

## Estrutura de arquivos

| Arquivo | Responsabilidade |
|---|---|
| `src/mutation_hunter/build.py` | entradas, referência, consequências, genótipos, tabelas, montagem do pacote, CLI de etapas |
| `src/mutation_hunter/store.py` | leitura do pacote (`Panel`), download e cache (`fetch_store`) |
| `nextflow/bin/annotate.sh` | `bcftools view | norm | csq` de uma região — os mesmos flags servirão à cepa externa |
| `nextflow/build.nf` | orquestração: CHECK → REFERENCE, STRAINS → ANNOTATE ×16 → GENOTYPES ×16 → ASSEMBLE |
| `environment.yml` | ambiente dos mantenedores (python, bcftools, nextflow, extras) |
| `.github/workflows/tests.yml` | job `test` (pip) + job `tools` (micromamba com bcftools) |
| `tests/conftest.py` | + `write_bed` e fixture `store` (pacote sintético construído por `assemble`) |
| `tests/test_build.py`, `tests/test_build_tools.py`, `tests/test_store.py` | testes das peças acima |
| `docs/build-check/` | verificação manual contra o proteoma público |

---

### Tarefa 1: Entradas, referência e consequências

**Arquivos:**
- Modificar: `pyproject.toml` (extra `build`)
- Criar: `src/mutation_hunter/build.py`, `tests/test_build.py`

**Interfaces:**
- Produz: `build.INPUTS: dict[str, str]`, `build.file_hash(path, algo="md5") -> str`, `build.check_inputs(folder) -> dict`, `build.reference(fasta_gz, gff_gz, out) -> None` (grava `ref.fa`, `ref.fa.fai`, `ref.gff3`, `aliases.tsv`), `build.write_fai(fa)`, `build.protein_lengths(gff) -> dict[str, int]`, `build.gene_spans(gff) -> DataFrame` (índice `gene`; colunas `chrom, start, end`), `build.parse_bcsq(str) -> list[(consequence, gene, aa_change)]`, `build.classify(consequence, aa_pos, protein_len) -> "missense"|"lof"|"inframe"|"synonymous"|"other"`, e as constantes `CHROMS`, `RENAME`, `LOF`, `INFRAME`, `SYSTEMATIC`, `VARIANT_COLUMNS`.

- [ ] **Passo 1: Criar o branch e o extra**

```bash
git switch baseline/proteome && git switch -c feat/panel-build
```

Em `pyproject.toml`, logo abaixo de `test = ["pytest>=8"]`, acrescentar:
```toml
build = ["cyvcf2>=0.31", "xlrd>=2"]
```

- [ ] **Passo 2: Escrever os testes que falham** — `tests/test_build.py`:

```python
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
```

- [ ] **Passo 3: Rodar e ver falhar**

Run: `uv run --extra test --extra build pytest tests/test_build.py -v`
Esperado: erro de coleta `ImportError: cannot import name 'build' from 'mutation_hunter'`.

- [ ] **Passo 4: Implementar** — `src/mutation_hunter/build.py` (cabeçalho completo; as tarefas seguintes acrescentam funções no fim do arquivo):

```python
"""Build the panel store from the public 1011 genomes data. Maintainers only.

Steps are called by nextflow/build.nf; bcftools does normalisation and consequence annotation upstream.
Needs the optional extra `build` (cyvcf2, xlrd)."""

import argparse
import gzip
import hashlib
import json
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import connected_components

from . import __version__

ROMAN = "I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI".split()
CHROMS = [f"chromosome{i}" for i in range(1, 17)]
RENAME = dict(zip(ROMAN, CHROMS))
LOF = {"stop_gained", "frameshift", "splice_acceptor", "splice_donor", "start_lost"}
INFRAME = {"inframe_insertion", "inframe_deletion"}
SYSTEMATIC = re.compile(r"^Y[A-P][LR]\d{3}[WC](-[A-Z])?$")
VARIANT_COLUMNS = ["chrom", "pos", "ref", "alt", "old_rec", "gene", "transcript", "consequence", "class",
                   "aa_change", "aa_pos", "protein_len", "ac", "an", "n_het", "n_hom", "call_rate", "ref_is_minor"]
# Upstream inputs and their MD5 (md5.txt of 1002genomes on 2025-04-01, Ensembl release 115, Nature supplement).
INPUTS = {
    "1011Matrix.gvcf.gz": "42478e3e9dff4bd46993d82d8eab40d4",
    "Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz": "5ae49d2d22cd977abdd6b054b94c544d",
    "Saccharomyces_cerevisiae.R64-1-1.115.gff3.gz": "17ed999c1853e03562f89bbcaefa80b6",
    "phenoMatrix_35ConditionsNormalizedByYPD.tab.gz": "479862335552a050e46de163d1ee5f72",
    "genesMatrix_PresenceAbsence.tab.gz": "b3b96719de226eaa277589d8d01f9202",
    "genesMatrix_CopyNumber.tab.gz": "c36cdf023c554580fe5b67101742fc09",
    "allORFs_pangenome.fasta.gz": "4f1ebba402ae64e52f2af72fa6942f38",
    "1011GWASMatrix.tar.gz": "c68c2ee5d9f4b5c85b8436a4a2f31ba5",
    "41586_2018_30_MOESM3_ESM.xls": "23d0ff1ec5c13a6b2e20a8076eb588c1",
}


def file_hash(path, algo="md5") -> str:
    h = hashlib.new(algo)
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check_inputs(folder) -> dict:
    """Verify every upstream input in folder against its pinned MD5; return {name: md5}."""
    folder = Path(folder)
    for name, md5 in INPUTS.items():
        if not (folder / name).exists():
            raise FileNotFoundError(f"missing input {name} in {folder}")
        if (got := file_hash(folder / name)) != md5:
            raise ValueError(f"{name}: MD5 {got} != pinned {md5}; upstream changed or file is corrupt")
    return dict(INPUTS)


def reference(fasta_gz, gff_gz, out) -> None:
    """R64 without Mito, chromosomes renamed I..XVI -> chromosome1..16 (the panel's names), plus aliases.tsv."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    keep = False
    with gzip.open(fasta_gz, "rt") as f, open(out / "ref.fa", "w") as o:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if keep := name in RENAME:
                    o.write(f">{RENAME[name]}\n")
            elif keep:
                o.write(line)
    with gzip.open(gff_gz, "rt") as f, open(out / "ref.gff3", "w") as o:
        for line in f:
            if line.startswith("##sequence-region"):
                p = line.split()
                if p[1] in RENAME:
                    o.write(" ".join([p[0], RENAME[p[1]], *p[2:]]) + "\n")
            elif line.startswith("#"):
                o.write(line)
            elif (col := line.split("\t", 1))[0] in RENAME:
                o.write(RENAME[col[0]] + "\t" + col[1])
    rows = [(a, c) for r, c in RENAME.items() for a in (r, f"chr{r}", c, c.replace("chromosome", "chr"),
                                                         c.replace("chromosome", ""))]
    pd.DataFrame(rows, columns=["alias", "chrom"]).to_csv(out / "aliases.tsv", sep="\t", index=False)
    write_fai(out / "ref.fa")


def write_fai(fa) -> None:
    """samtools-style .fai, written once here so parallel bcftools jobs never race to create it."""
    rows, offset = [], 0
    with open(fa, "rb") as f:
        for line in f:
            if line.startswith(b">"):
                rows.append([line[1:].split()[0].decode(), 0, offset + len(line), 0, 0])
            else:
                row = rows[-1]
                if not row[3]:
                    row[3], row[4] = len(line.rstrip(b"\n")), len(line)
                row[1] += len(line.rstrip(b"\n"))
            offset += len(line)
    Path(f"{fa}.fai").write_text("".join("\t".join(map(str, r)) + "\n" for r in rows))


def protein_lengths(gff) -> dict[str, int]:
    """Transcript -> protein length in amino acids (CDS length / 3, stop codon excluded)."""
    cds = {}
    for line in open(gff):
        col = line.split("\t")
        if len(col) > 8 and col[2] == "CDS":
            t = re.search(r"Parent=transcript:([^;\n]+)", col[8])[1]
            cds[t] = cds.get(t, 0) + int(col[4]) - int(col[3]) + 1
    return {t: n // 3 - 1 for t, n in cds.items()}


def gene_spans(gff) -> pd.DataFrame:
    rows = []
    for line in open(gff):
        col = line.split("\t")
        if len(col) > 8 and col[2] == "gene" and "biotype=protein_coding" in col[8]:
            rows.append((re.search(r"ID=gene:([^;\n]+)", col[8])[1], col[0], int(col[3]), int(col[4])))
    return pd.DataFrame(rows, columns=["gene", "chrom", "start", "end"]).set_index("gene")


def parse_bcsq(bcsq: str) -> list[tuple[str, str, str]]:
    """BCSQ -> [(consequence, systematic gene, amino-acid change)] for protein-coding transcripts only."""
    out = []
    for item in bcsq.split(","):
        f = item.split("|")
        if len(f) >= 4 and f[2].endswith("_mRNA") and f[3] == "protein_coding":
            out.append((f[0], f[2][:-5], f[5] if len(f) > 5 else ""))
    return out


def classify(consequence: str, aa_pos: int, protein_len: int) -> str:
    """missense / lof / inframe / synonymous / other. Truncations in the last 5% of the protein are not LoF."""
    parts = set(consequence.split("&"))
    if parts & LOF:
        late = parts & {"stop_gained", "frameshift"} and not parts & {"start_lost", "splice_acceptor", "splice_donor"}
        return "other" if late and protein_len and aa_pos > 0.95 * protein_len else "lof"
    if parts & INFRAME:
        return "inframe"
    if "missense" in parts:
        return "missense"
    return "synonymous" if "synonymous" in parts else "other"
```

- [ ] **Passo 5: Rodar e ver passar**

Run: `uv run --extra test --extra build pytest tests/test_build.py -v`
Esperado: 13 passed (o `test_classify` tem 10 casos).

- [ ] **Passo 6: Commit**

```bash
git add pyproject.toml uv.lock src/mutation_hunter/build.py tests/test_build.py \
    docs/superpowers/specs/2026-10-15-panel-build-design.md docs/superpowers/plans/2026-10-15-panel-build.md
git commit -m "Panel build: pinned inputs, renamed R64 reference and consequence classes"
```

---

### Tarefa 2: Genótipos por cromossomo

**Arquivos:**
- Modificar: `src/mutation_hunter/build.py` (acrescentar no fim), `tests/test_build.py` (acrescentar no fim)
- Criar: `nextflow/bin/annotate.sh`, `tests/test_build_tools.py`

**Interfaces:**
- Consome: `parse_bcsq`, `classify`, `protein_lengths`, `VARIANT_COLUMNS` (Tarefa 1).
- Produz: `build.encode(gt_types, gq, dp, haploid, min_gq=20, min_dp=5) -> int8 array`, `build.site_stats(code) -> dict` (`ac, an, n_het, n_hom, call_rate, ref_is_minor`), `build.standardised_carriers(code, min_freq=0.05, min_call=0.95) -> array | None`, `build.genotypes(bcf, strains_tsv, gff, out_prefix, min_gq=20, min_dp=5)` que grava `<out_prefix>.variants.tsv.gz`, `<out_prefix>.geno.npz` (CSR linhas × cepas) e `<out_prefix>.kin.npz` (`zz`, `m`). `strains_tsv` precisa da coluna `ploidy` e das cepas na ordem do VCF.
- Produz `nextflow/bin/annotate.sh <vcf.gz> <região> <ref.fa> <ref.gff3> <out.bcf>`.

- [ ] **Passo 1: Escrever os testes que falham**

Acrescentar a `tests/test_build.py`:
```python


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
```

`tests/test_build_tools.py`:
```python
"""bcftools + cyvcf2 on a hand-written VCF with known consequences. Skipped when the tools are missing."""

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from mutation_hunter import build

pytest.importorskip("cyvcf2")
pytestmark = pytest.mark.skipif(not shutil.which("bcftools"), reason="bcftools not installed")
ANNOTATE = Path(__file__).parents[1] / "nextflow" / "bin" / "annotate.sh"
COMP = str.maketrans("ACGT", "TGCA")
SAMPLES = ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "SH"]  # SH is haploid


def orf(changes: dict[int, str]) -> str:
    """100 codons: ATG, 98 x GCT (Ala) with overrides, TAA. Protein length 99 aa."""
    codons = ["ATG"] + ["GCT"] * 98 + ["TAA"]
    for k, c in changes.items():
        codons[k - 1] = c
    return "".join(codons)


def reference() -> str:
    seq = list(np.random.default_rng(5).choice(list("ACGT"), 2000))
    seq[100:400] = orf({20: "TCA", 97: "TCA"})  # YAL901W, + strand, CDS 101-400
    seq[1000:1300] = orf({}).translate(COMP)[::-1]  # YAL902C, - strand, CDS 1001-1300
    return "".join(seq)


def pos(codon: int, base: int) -> int:  # 1-based chromosome position in YAL901W
    return 101 + 3 * (codon - 1) + (base - 1)


GT = "0/0:99:30"
RECORDS = [  # pos, ref, alt, genotypes for S1..S7, SH ; expected (class, aa_change prefix) per alt
    (pos(10, 1), "G", "A", ["0/1:99:30", "1/1:99:30", GT, GT, GT, GT, GT, "1/1:99:30"]),  # missense 10A>10T
    (pos(20, 2), "C", "A", ["1/1:99:30", GT, GT, GT, GT, GT, GT, GT]),  # stop_gained, early -> lof
    (pos(30, 3), "T", "TG", [GT, "1/1:99:30", GT, GT, GT, GT, GT, GT]),  # frameshift -> lof
    (pos(40, 3), "T", "C", [GT, GT, "1/1:99:30", GT, GT, GT, GT, GT]),  # synonymous
    (pos(50, 1), "G", "A,C", [GT, GT, GT, "1/2:99:30", GT, GT, GT, GT]),  # multiallelic: two missense rows
    (pos(97, 2), "C", "A", [GT, GT, GT, GT, "1/1:99:30", GT, GT, GT]),  # stop in last 5% -> other
    (pos(60, 1), "G", "A", ["1/1:99:30", GT, GT, GT, GT, "0/1:10:30", "0/1:99:3", "0/1:99:30"]),  # low GQ, low DP, haploid het
    (1273, "C", "T", [GT, GT, GT, GT, GT, "1/1:99:30", GT, GT]),  # missense on the - strand gene
]


@pytest.fixture
def built(tmp_path):
    ref = reference()
    for p, r, *_ in RECORDS:
        assert ref[p - 1:p - 1 + len(r)] == r, p
    (tmp_path / "ref.fa").write_text(f">chromosome1\n{ref}\n")
    (tmp_path / "ref.gff3").write_text(
        "##gff-version 3\n" + "".join(
            f"chromosome1\tt\tgene\t{a}\t{b}\t.\t{s}\t.\tID=gene:{g};biotype=protein_coding\n"
            f"chromosome1\tt\tmRNA\t{a}\t{b}\t.\t{s}\t.\tID=transcript:{g}_mRNA;Parent=gene:{g};biotype=protein_coding\n"
            f"chromosome1\tt\texon\t{a}\t{b}\t.\t{s}\t.\tParent=transcript:{g}_mRNA\n"
            f"chromosome1\tt\tCDS\t{a}\t{b}\t.\t{s}\t0\tParent=transcript:{g}_mRNA\n"
            for g, a, b, s in [("YAL901W", 101, 400, "+"), ("YAL902C", 1001, 1300, "-")]))
    header = ("##fileformat=VCFv4.2\n##contig=<ID=chromosome1,length=2000>\n"
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
              '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="q">\n'
              '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="d">\n'
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(SAMPLES) + "\n")
    body = "".join(f"chromosome1\t{p}\t.\t{r}\t{a}\t100\t.\t.\tGT:GQ:DP\t" + "\t".join(g) + "\n"
                   for p, r, a, g in sorted(RECORDS))
    (tmp_path / "in.vcf").write_text(header + body)
    subprocess.run(["bcftools", "view", "-Oz", "-o", tmp_path / "in.vcf.gz", tmp_path / "in.vcf"], check=True)
    subprocess.run(["bcftools", "index", "-t", tmp_path / "in.vcf.gz"], check=True)
    subprocess.run([ANNOTATE, tmp_path / "in.vcf.gz", "chromosome1", tmp_path / "ref.fa", tmp_path / "ref.gff3",
                    tmp_path / "out.bcf"], check=True)
    pd.DataFrame({"ploidy": [2] * 7 + [1]}, index=pd.Index(SAMPLES, name="strain")).to_csv(tmp_path / "strains.tsv", sep="\t")
    build.genotypes(tmp_path / "out.bcf", tmp_path / "strains.tsv", tmp_path / "ref.gff3", tmp_path / "chromosome1")
    var = pd.read_csv(tmp_path / "chromosome1.variants.tsv.gz", sep="\t", keep_default_na=False)
    return var, sparse.load_npz(tmp_path / "chromosome1.geno.npz").toarray(), tmp_path


def test_consequence_classes(built):
    var, _, _ = built
    got = list(zip(var["pos"], var["gene"], var["class"]))
    assert got == [
        (pos(10, 1), "YAL901W", "missense"), (pos(20, 2), "YAL901W", "lof"), (pos(30, 3), "YAL901W", "lof"),
        (pos(40, 3), "YAL901W", "synonymous"), (pos(50, 1), "YAL901W", "missense"),
        (pos(50, 1), "YAL901W", "missense"), (pos(60, 1), "YAL901W", "missense"), (pos(97, 2), "YAL901W", "other"),
        (1273, "YAL902C", "missense")]
    assert var["aa_change"][0].startswith("10A>10T") and var["protein_len"][0] == 99


def test_multiallelic_split_keeps_origin_and_genotypes(built):
    var, geno, _ = built
    rows = var.index[var["pos"] == pos(50, 1)]
    assert set(var.loc[rows, "alt"]) == {"A", "C"} and var.loc[rows, "old_rec"].str.contains("A,C").all()
    assert (geno[rows, 3] > 0).all()  # S4 (1/2) carries both


def test_genotype_codes_and_site_stats(built):
    var, geno, _ = built
    assert geno[0].tolist() == [1, 2, 0, 0, 0, 0, 0, 2]  # haploid hom alt stays
    assert (var.loc[0, "ac"], var.loc[0, "an"], var.loc[0, "n_het"]) == (5, 16, 1)


def test_quality_rules_set_missing(built):
    var, geno, _ = built
    row = var.index[var["pos"] == pos(60, 1)][0]
    assert geno[row].tolist() == [2, 0, 0, 0, 0, -1, -1, -1]  # GQ 10, DP 3, het in haploid
    assert var.loc[row, "call_rate"] == 5 / 8
```

- [ ] **Passo 2: Rodar e ver falhar**

Run: `uv run --extra test --extra build pytest tests/test_build.py tests/test_build_tools.py -v`
Esperado: `AttributeError: module 'mutation_hunter.build' has no attribute 'encode'` (e os testes de ferramentas falham por falta de `annotate.sh`, ou são pulados se não houver `bcftools` no `PATH`; neste caso, coloque o `bcftools` no `PATH` antes de seguir).

- [ ] **Passo 3: Implementar**

`nextflow/bin/annotate.sh` (depois `chmod +x nextflow/bin/annotate.sh`):
```bash
#!/usr/bin/env bash
# Normalise and annotate one region of the panel VCF. Same flags for the panel and for external strains.
# usage: annotate.sh <vcf.gz> <region> <ref.fa> <ref.gff3> <out.bcf>
set -euo pipefail
bcftools view -r "$2" "$1" -Ou \
  | bcftools norm -f "$3" -c x -m -any -a --atom-overlaps . --old-rec-tag OLD_REC -Ou 2>/dev/null \
  | bcftools csq -f "$3" -g "$4" -l -s - -v 0 -Ob -o "$5"
bcftools index -f "$5"
```

Acrescentar a `src/mutation_hunter/build.py`:
```python


def encode(gt_types, gq, dp, haploid, min_gq=20, min_dp=5) -> np.ndarray:
    """cyvcf2 gt_types (0 ref, 1 het, 2 unknown, 3 hom alt) -> int8: -1 missing, 0 ref, 1 het, 2 hom alt.
    Low GQ/DP and heterozygous calls in haploid strains become missing."""
    code = np.select([gt_types == 1, gt_types == 3], [1, 2], 0).astype(np.int8)
    bad = (gt_types == 2) | (gq < min_gq) | (dp < min_dp) | ((gt_types == 1) & haploid)
    code[bad] = -1
    return code


def site_stats(code: np.ndarray) -> dict:
    called = code >= 0
    n_het, n_hom, an = int((code == 1).sum()), int((code == 2).sum()), int(2 * called.sum())
    ac = n_het + 2 * n_hom
    return {"ac": ac, "an": an, "n_het": n_het, "n_hom": n_hom, "call_rate": called.mean(),
            "ref_is_minor": bool(an and ac / an > 0.5)}


def standardised_carriers(code: np.ndarray, min_freq=0.05, min_call=0.95) -> np.ndarray | None:
    """Common-SNP carrier vector standardised for the GRM (missing -> 0, i.e. mean-imputed), or None."""
    called = code >= 0
    if called.mean() < min_call:
        return None
    f = (code[called] > 0).mean()
    if not min_freq <= f <= 1 - min_freq:
        return None
    return np.where(called, (code > 0) - f, 0.0) / np.sqrt(f * (1 - f))


def genotypes(bcf, strains_tsv, gff, out_prefix, min_gq=20, min_dp=5) -> None:
    """One chromosome: <out>.variants.tsv.gz, <out>.geno.npz (variant rows x strains) and <out>.kin.npz
    (sum of Z'Z over common SNPs and their count, for the GRM and its LOCO versions)."""
    from cyvcf2 import VCF

    strains = pd.read_csv(strains_tsv, sep="\t", index_col=0)
    vcf = VCF(str(bcf), strict_gt=True)
    if list(vcf.samples) != list(strains.index):
        raise ValueError("VCF samples and strains.tsv differ in content or order")
    haploid = (strains["ploidy"] == 1).to_numpy()
    plen = protein_lengths(gff)
    rows, mats, block = [], [], []
    zz, m = np.zeros((len(strains),) * 2), 0
    for rec in vcf:
        code = encode(rec.gt_types, rec.gt_quals, rec.format("DP")[:, 0], haploid, min_gq, min_dp)
        if len(rec.REF) == 1 and len(rec.ALT[0]) == 1 and (z := standardised_carriers(code)) is not None:
            block.append(z)
            if len(block) == 2000:
                Z = np.array(block)
                zz += Z.T @ Z
                m += len(block)
                block = []
        ann = parse_bcsq(rec.INFO.get("BCSQ") or "")
        if not ann:
            continue
        stats_ = site_stats(code)
        nz = sparse.csr_matrix(code)
        for consequence, gene, aa in ann:
            aa_pos = int(mt[1]) if (mt := re.match(r"(\d+)", aa)) else 0
            rows.append({"chrom": rec.CHROM, "pos": rec.POS, "ref": rec.REF, "alt": rec.ALT[0],
                         "old_rec": rec.INFO.get("OLD_REC") or "", "gene": gene, "transcript": gene + "_mRNA",
                         "consequence": consequence, "class": classify(consequence, aa_pos, plen.get(gene + "_mRNA", 0)),
                         "aa_change": aa, "aa_pos": aa_pos, "protein_len": plen.get(gene + "_mRNA", 0), **stats_})
            mats.append(nz)  # ponytail: overlapping genes duplicate the genotype row; fine at yeast scale
    if block:
        Z = np.array(block)
        zz += Z.T @ Z
        m += len(block)
    pd.DataFrame(rows, columns=VARIANT_COLUMNS).to_csv(f"{out_prefix}.variants.tsv.gz", sep="\t", index=False)
    geno = sparse.vstack(mats, format="csr") if mats else sparse.csr_matrix((0, len(strains)), dtype=np.int8)
    sparse.save_npz(f"{out_prefix}.geno.npz", geno.astype(np.int8))
    np.savez(f"{out_prefix}.kin.npz", zz=zz, m=m)
```

- [ ] **Passo 4: Rodar e ver passar**

Run: `uv run --extra test --extra build pytest tests/test_build.py tests/test_build_tools.py -v`
Esperado: 16 + 4 passed.

- [ ] **Passo 5: Commit**

```bash
git add nextflow/bin/annotate.sh src/mutation_hunter/build.py tests/test_build.py tests/test_build_tools.py
git commit -m "Panel build: per-chromosome annotation, genotype quality rules and kinship sums"
```

---

### Tarefa 3: Tabelas de cepas, de ORFs e de posicionamento

**Arquivos:**
- Modificar: `src/mutation_hunter/build.py`, `tests/test_build.py`, `tests/conftest.py` (acrescentar no fim)

**Interfaces:**
- Produz: `build.strain_table(xls, samples, phenotypes_gz) -> DataFrame` (índice `strain`; colunas `isolate, ploidy, zygosity, aneuploidy, clade, ecological_origin, mta_restricted, has_phenotype`; erro se alguma cepa do VCF faltar na Tabela S1), `build.orf_table(pa_columns, pangenome_fasta_gz) -> DataFrame` (`column, orf, family, gene`), `build.read_plink_bed(prefix) -> (bim, int8 sites × amostras)`.
- Produz em `conftest.py`: `write_bed(prefix, dosage, chroms, samples)` (só para testes).

O número no nome da coluna (`X7051.YOR083W`) liga a coluna desfigurada pelo R ao cabeçalho `7051-YOR083W` do FASTA do pangenoma; é assim que `YAL016C.A` volta a ser `YAL016C-A`. Na base real: 7.796/7.796 colunas resolvidas; 5.267 dos 6.600 genes da Ensembl têm coluna (o resto são famílias `_NumOfGenes_n` ou genes sem coluna).

- [ ] **Passo 1: Escrever os testes que falham**

Em `tests/conftest.py`: trocar a docstring da primeira linha por
`"""Fixtures: a 40-strain fake proteome panel (baseline analysis) and a tiny built panel store (no tools)."""`;
acrescentar aos imports `import pandas as pd` e `from scipy import sparse`; trocar `from mutation_hunter import data`
por `from mutation_hunter import build, data`. Depois, acrescentar no fim:
```python


# --- tiny built panel ------------------------------------------------------------------------------------

STORE_STRAINS = ["S1", "S2", "S3", "S4", "S5", "SACE_YAA"]
GFF = """##gff-version 3
chromosome1\tt\tgene\t101\t400\t.\t+\t.\tID=gene:YAL901W;biotype=protein_coding
chromosome1\tt\tmRNA\t101\t400\t.\t+\t.\tID=transcript:YAL901W_mRNA;Parent=gene:YAL901W;biotype=protein_coding
chromosome1\tt\tCDS\t101\t400\t.\t+\t0\tParent=transcript:YAL901W_mRNA
chromosome2\tt\tgene\t1001\t1300\t.\t-\t.\tID=gene:YBL902C;biotype=protein_coding
chromosome2\tt\tmRNA\t1001\t1300\t.\t-\t.\tID=transcript:YBL902C_mRNA;Parent=gene:YBL902C;biotype=protein_coding
chromosome2\tt\tCDS\t1001\t1300\t.\t-\t0\tParent=transcript:YBL902C_mRNA
chromosome2\tt\tgene\t1500\t1800\t.\t+\t.\tID=gene:YBL903W;biotype=protein_coding
"""


def write_bed(prefix, dosage: np.ndarray, chroms: list[int], samples: list[str]) -> None:
    """Minimal PLINK 1 SNP-major writer (dosage = A1 copies, -1 missing), for tests only."""
    code = {2: 0, -1: 1, 1: 2, 0: 3}
    n = dosage.shape[1]
    out = bytearray([0x6C, 0x1B, 0x01])
    for row in dosage:
        b = [code[int(x)] for x in row] + [0] * (-n % 4)
        out += bytes(b[i] | b[i + 1] << 2 | b[i + 2] << 4 | b[i + 3] << 6 for i in range(0, len(b), 4))
    open(f"{prefix}.bed", "wb").write(bytes(out))
    pd.DataFrame({"c": chroms, "id": [f"rs{i}" for i in range(len(chroms))], "cm": 0, "pos": range(1, len(chroms) + 1),
                  "a1": "A", "a2": "G"}).to_csv(f"{prefix}.bim", sep="\t", header=False, index=False)
    open(f"{prefix}.fam", "w").write("".join(f"{s} {s} 0 0 0 -9\n" for s in samples))
```

Em `tests/test_build.py`, acrescentar ao bloco de imports `from conftest import write_bed` e, no fim:
```python


def test_orf_table_uses_the_number_to_recover_real_names(tmp_path):
    fa = tmp_path / "orfs.fa.gz"
    fa.write_bytes(gzip.compress(b">1-EC1118_1F14_0012g\nA\n>7-YAL016C-A\nA\n>9-YAL005C_NumOfGenes_3\nA\n"))
    t = build.orf_table(["X1.EC1118_1F14_0012g", "X7.YAL016C.A", "X9.YAL005C_NumOfGenes_3"], fa)
    assert t["orf"].tolist() == ["EC1118_1F14_0012g", "YAL016C-A", "YAL005C_NumOfGenes_3"]
    assert t["gene"].tolist()[1:] == ["YAL016C-A", "YAL005C"] and pd.isna(t["gene"][0])
    assert t["family"].tolist() == [False, False, True]
    with pytest.raises(ValueError, match="no pangenome ORF"):
        build.orf_table(["X2.SOMETHING"], fa)


def test_read_plink_bed_roundtrip(tmp_path):
    dos = np.array([[0, 1, 2, -1, 2], [2, 2, 0, 0, 1]])
    write_bed(tmp_path / "x", dos, [1, 17], list("abcde"))
    bim, got = build.read_plink_bed(tmp_path / "x")
    assert got.tolist() == dos.tolist() and bim["chrom"].tolist() == [1, 17]
```

- [ ] **Passo 2: Rodar e ver falhar**

Run: `uv run --extra test --extra build pytest tests/test_build.py -v`
Esperado: `AttributeError: ... 'orf_table'` e `... 'read_plink_bed'`.

- [ ] **Passo 3: Implementar** — acrescentar a `src/mutation_hunter/build.py`:

```python


def strain_table(xls, samples: list[str], phenotypes_gz) -> pd.DataFrame:
    """Peter 2018 Table S1, in VCF sample order."""
    s1 = pd.read_excel(xls, "Table S1", header=3)
    s1 = s1[s1["Standardized name"].notna()].set_index("Standardized name")
    pheno = pd.read_csv(phenotypes_gz, sep="\t", index_col=0)
    t = pd.DataFrame({
        "isolate": s1["Isolate name"],
        "ploidy": pd.to_numeric(s1["Ploidy"], errors="coerce").astype("Int64"),  # 'Xn' -> unknown
        "zygosity": s1["Zygosity"],
        "aneuploidy": s1["Aneuploidies"],
        "clade": s1["Clades"].astype(str).str.strip(),
        "ecological_origin": s1["Ecological origins"],
        "mta_restricted": s1["MTA (no distribution)"].eq("yes"),
    }).reindex(samples)
    if t["isolate"].isna().any():
        raise ValueError(f"strains missing from Table S1: {t.index[t['isolate'].isna()].tolist()[:5]}")
    t.index.name = "strain"
    t["has_phenotype"] = t.index.isin(pheno.index)
    return t


def orf_table(pa_columns, pangenome_fasta_gz) -> pd.DataFrame:
    """Matrix column 'X7051.YOR083W' -> pangenome ORF 7051 -> real name and systematic gene (if any).
    The number links the R-mangled column to the FASTA header '7051-YOR083W'."""
    names = {}
    with gzip.open(pangenome_fasta_gz, "rt") as f:
        for line in f:
            if line.startswith(">"):
                num, name = line[1:].strip().split("-", 1)
                names[int(num)] = name
    num = pd.Series(pa_columns).str.extract(r"^X(\d+)\.")[0].astype(int)
    orf = num.map(names)
    if orf.isna().any():
        raise ValueError(f"{orf.isna().sum()} matrix columns have no pangenome ORF")
    base = orf.str.replace(r"_NumOfGenes_\d+$", "", regex=True)
    return pd.DataFrame({"column": list(pa_columns), "orf": orf, "family": orf != base,
                         "gene": base.where(base.str.match(SYSTEMATIC.pattern))})


def read_plink_bed(prefix) -> tuple[pd.DataFrame, np.ndarray]:
    """PLINK 1 SNP-major .bed -> (bim table, int8 sites x samples: alt-allele count 0/1/2, -1 missing)."""
    bim = pd.read_csv(f"{prefix}.bim", sep=r"\s+", header=None, names=["chrom", "id", "cm", "pos", "a1", "a2"])
    n = sum(1 for _ in open(f"{prefix}.fam"))
    raw = np.fromfile(f"{prefix}.bed", dtype=np.uint8)
    if raw[:3].tolist() != [0x6C, 0x1B, 0x01]:
        raise ValueError("not a SNP-major PLINK 1 .bed")
    packed = raw[3:].reshape(len(bim), -1)
    bits = np.stack([(packed >> s) & 3 for s in (0, 2, 4, 6)], axis=2).reshape(len(bim), -1)[:, :n]
    # PLINK codes: 0 hom A1, 1 missing, 2 het, 3 hom A2; count A1 copies.
    lut = np.array([2, -1, 1, 0], dtype=np.int8)
    return bim, lut[bits]
```

`strain_table` lê o Excel do suplemento e não tem teste unitário (exigiria escrever `.xls`); ela é exercitada na construção real (Tarefa 6), que falha se alguma cepa do VCF faltar.

- [ ] **Passo 4: Rodar e ver passar**

Run: `uv run --extra test --extra build pytest tests/test_build.py -v`
Esperado: 18 passed.

- [ ] **Passo 5: Commit**

```bash
git add src/mutation_hunter/build.py tests/test_build.py tests/conftest.py
git commit -m "Panel build: strain metadata, pangenome ORF names and PLINK placement sites"
```

---

### Tarefa 4: Montagem do pacote

**Arquivos:**
- Modificar: `src/mutation_hunter/build.py`, `tests/test_build.py`, `tests/conftest.py`

**Interfaces:**
- Consome: tudo das Tarefas 1–3.
- Produz: `build.kinships(parts) -> (K, {chrom: K_loco})`, `build.clone_groups(K, frac=0.01) -> int array`, `build.assemble(parts_dir, inputs_dir, reference_dir, strains_tsv, out, input_md5: dict, params: dict)` (`params` precisa de `panel_version`). Escreve em `<out>.part` e renomeia no fim.
- Layout do pacote (formato 1): `manifest.json`, `SHA256SUMS`, `strains.tsv` (+ `clone_group`), `variants.tsv.gz` (+ coluna `row`, ordenado por gene), `genotypes/{data,indices,indptr}.npy` + `shape.json`, `genes.tsv` (`chrom, start, end, first, last, orf`; `first = -1` sem variantes), `kinship.npy`, `kinship_loco/<chrom>.npy`, `phenotypes.tsv.gz`, `orfs.tsv`, `presence_absence.npy`, `copy_number.npy`, `placement.npy`, `placement_sites.tsv`, `reference/`.
- Produz em `conftest.py`: fixture `store` (pacote sintético de 6 cepas, 3 linhas de variante em 2 cromossomos).

- [ ] **Passo 1: Escrever os testes que falham**

Acrescentar a `tests/conftest.py`:
```python


@pytest.fixture
def store(tmp_path):
    """Run build.assemble on hand-made per-chromosome parts; returns the store folder."""
    n = len(STORE_STRAINS)
    parts, inputs, ref = tmp_path / "parts", tmp_path / "inputs", tmp_path / "reference"
    for d in (parts, inputs, ref):
        d.mkdir()
    var = {"chromosome1": [("YAL901W", 128, "missense", [0, 1, 2, 0, -1, 2]),
                           ("YAL901W", 159, "lof", [2, 0, 0, 0, 0, 0])],
           "chromosome2": [("YBL902C", 1273, "missense", [0, 0, 0, 1, 1, 0])]}
    rng = np.random.default_rng(1)
    for chrom, rows in var.items():
        pd.DataFrame([{"chrom": chrom, "pos": pos, "ref": "G", "alt": "A", "old_rec": "", "gene": g,
                       "transcript": g + "_mRNA", "consequence": c, "class": c, "aa_change": "", "aa_pos": 0,
                       "protein_len": 99, **build.site_stats(np.array(code, np.int8))}
                      for g, pos, c, code in rows], columns=build.VARIANT_COLUMNS).to_csv(
            parts / f"{chrom}.variants.tsv.gz", sep="\t", index=False)
        sparse.save_npz(parts / f"{chrom}.geno.npz", sparse.csr_matrix(np.array([r[3] for r in rows], np.int8)))
        Z = rng.normal(size=(50, n))
        np.savez(parts / f"{chrom}.kin.npz", zz=Z.T @ Z, m=50)
    (ref / "ref.gff3").write_text(GFF)
    (ref / "ref.fa").write_text(">chromosome1\nACGT\n")
    (ref / "aliases.tsv").write_text("alias\tchrom\nI\tchromosome1\n")
    strains = pd.DataFrame({"isolate": STORE_STRAINS, "ploidy": [2, 2, 2, 2, 1, 2]}, index=pd.Index(STORE_STRAINS, name="strain"))
    strains.to_csv(tmp_path / "strains.tsv", sep="\t")

    pheno = pd.DataFrame({"COND1": rng.normal(size=n)}, index=STORE_STRAINS)
    (inputs / "phenoMatrix_35ConditionsNormalizedByYPD.tab.gz").write_bytes(gzip.compress(pheno.to_csv(sep="\t").encode()))
    cols = ["X1.EC1118_1F14_0012g", "X2.YAL901W", "X3.YAL005C_NumOfGenes_3"]
    pa = pd.DataFrame([[1, 1, 1]] * n, index=pd.Index(STORE_STRAINS, name="STD_name"), columns=cols)
    pa.iloc[4, 1] = 0
    for name, m in (("genesMatrix_PresenceAbsence.tab.gz", pa), ("genesMatrix_CopyNumber.tab.gz", pa * 1.5)):
        (inputs / name).write_bytes(gzip.compress(m.to_csv(sep="\t").encode()))
    (inputs / "allORFs_pangenome.fasta.gz").write_bytes(
        gzip.compress(b">1-EC1118_1F14_0012g\nATG\n>2-YAL901W\nATG\n>3-YAL005C_NumOfGenes_3\nATG\n"))
    gw = tmp_path / "gw"
    gw.mkdir()
    write_bed(gw / "1011GWAS_matrix", np.array([[0, 1, 2, -1, 0, 2], [2, 2, 0, 0, 1, 1], [1, 0, 0, 0, 0, 0]]),
              [1, 2, 17], STORE_STRAINS)
    with tarfile.open(inputs / "1011GWASMatrix.tar.gz", "w:gz") as t:
        for f in gw.iterdir():
            t.add(f, arcname=f.name)

    out = tmp_path / "panel"
    build.assemble(parts, inputs, ref, tmp_path / "strains.tsv", out, {"fake": "md5"}, {"panel_version": "0.0.1-test"})
    return out
```

Acrescentar a `tests/test_build.py`:
```python


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
    assert np.load(store / "placement.npy").shape == (6, 2)  # chromosome 17 (CNV) rows dropped
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
```

- [ ] **Passo 2: Rodar e ver falhar**

Run: `uv run --extra test --extra build pytest tests/test_build.py -v`
Esperado: `AttributeError: ... 'clone_groups'` e erro na fixture `store` (`... 'assemble'`).

- [ ] **Passo 3: Implementar** — acrescentar a `src/mutation_hunter/build.py`:

```python


def kinships(parts: list[Path]) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """GRM over all chromosomes and leave-one-chromosome-out versions (float32)."""
    zz = {p.name.split(".")[0]: np.load(p) for p in parts}
    total, m = sum(z["zz"] for z in zz.values()), sum(int(z["m"]) for z in zz.values())
    loco = {c: ((total - z["zz"]) / (m - int(z["m"]))).astype(np.float32) for c, z in zz.items()}
    return (total / m).astype(np.float32), loco


def clone_groups(K: np.ndarray, frac=0.01) -> np.ndarray:
    """Connected components of strains whose GRM distance is below frac x the median distance."""
    d = np.diag(K)[:, None] + np.diag(K)[None, :] - 2 * K
    thr = frac * np.median(d[np.triu_indices(len(K), 1)])
    return connected_components(sparse.csr_matrix(d < thr), directed=False)[1]


def assemble(parts_dir, inputs_dir, reference_dir, strains_tsv, out, input_md5: dict, params: dict) -> None:
    parts_dir, inputs, out = Path(parts_dir), Path(inputs_dir), Path(out)
    tmp = out.with_name(out.name + ".part")
    shutil.rmtree(tmp, ignore_errors=True)
    (tmp / "genotypes").mkdir(parents=True)
    chroms = [c for c in CHROMS if (parts_dir / f"{c}.variants.tsv.gz").exists()]
    var = pd.concat([pd.read_csv(parts_dir / f"{c}.variants.tsv.gz", sep="\t", keep_default_na=False,
                                 dtype={"old_rec": str, "aa_change": str}) for c in chroms], ignore_index=True)
    geno = sparse.vstack([sparse.load_npz(parts_dir / f"{c}.geno.npz") for c in chroms], format="csr")
    order = np.lexsort((var["alt"], var["pos"], var["gene"]))  # each gene becomes one contiguous block
    var, geno = var.iloc[order].reset_index(drop=True), geno[order]
    var.insert(0, "row", np.arange(len(var)))
    var.to_csv(tmp / "variants.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0})  # reproducible bytes
    for name in ("data", "indices", "indptr"):
        np.save(tmp / "genotypes" / f"{name}.npy", getattr(geno, name))
    (tmp / "genotypes" / "shape.json").write_text(json.dumps(list(geno.shape)))

    spans = gene_spans(Path(reference_dir) / "ref.gff3")
    blocks = var.groupby("gene")["row"].agg(first="min", last="max")
    genes = spans.join(blocks, how="left")
    genes["first"], genes["last"] = genes["first"].fillna(-1).astype(int), genes["last"].fillna(-1).astype(int)

    K, loco = kinships(sorted(parts_dir.glob("*.kin.npz")))
    np.save(tmp / "kinship.npy", K)
    (tmp / "kinship_loco").mkdir()
    for c, k in loco.items():
        np.save(tmp / "kinship_loco" / f"{c}.npy", k)
    strains = pd.read_csv(strains_tsv, sep="\t", index_col=0)
    strains["clone_group"] = clone_groups(K)
    strains.to_csv(tmp / "strains.tsv", sep="\t")

    shutil.copy(inputs / "phenoMatrix_35ConditionsNormalizedByYPD.tab.gz", tmp / "phenotypes.tsv.gz")
    pa = pd.read_csv(inputs / "genesMatrix_PresenceAbsence.tab.gz", sep="\t", index_col=0)
    if len(missing := strains.index.difference(pa.index)):
        raise ValueError(f"{len(missing)} VCF strains missing from the presence/absence matrix: {list(missing)[:5]}")
    pa = pa.reindex(strains.index)
    cn = pd.read_csv(inputs / "genesMatrix_CopyNumber.tab.gz", sep="\t", index_col=0).reindex(strains.index)
    orfs = orf_table(pa.columns, inputs / "allORFs_pangenome.fasta.gz")
    if list(cn.columns) != list(pa.columns):
        raise ValueError("presence/absence and copy-number matrices have different columns")
    orfs.to_csv(tmp / "orfs.tsv", sep="\t", index=False)
    np.save(tmp / "presence_absence.npy", pa.to_numpy(np.float32))
    np.save(tmp / "copy_number.npy", cn.to_numpy(np.float32))
    genes = genes.join(orfs.dropna(subset=["gene"]).drop_duplicates("gene").set_index("gene")["orf"])
    genes.to_csv(tmp / "genes.tsv", sep="\t")

    with __import__("tarfile").open(inputs / "1011GWASMatrix.tar.gz") as t:
        t.extractall(tmp / "_gwas", filter="data")
    prefix = next((tmp / "_gwas").rglob("*.bed")).with_suffix("")
    bim, dos = read_plink_bed(prefix)
    fam = [line.split()[1] for line in open(f"{prefix}.fam")]
    snp = (bim["chrom"] <= 16).to_numpy()  # chromosome 17 rows are CNVs
    sites = bim[snp].assign(chrom=lambda b: "chromosome" + b["chrom"].astype(str))
    sites.to_csv(tmp / "placement_sites.tsv", sep="\t", index=False)
    np.save(tmp / "placement.npy", pd.DataFrame(dos[snp].T, index=fam).reindex(strains.index).fillna(-1)
            .to_numpy(np.int8))
    shutil.rmtree(tmp / "_gwas")

    (tmp / "reference").mkdir()
    for name in ("ref.fa", "ref.gff3", "aliases.tsv"):
        shutil.copy(Path(reference_dir) / name, tmp / "reference" / name)
    shutil.copy(inputs / "allORFs_pangenome.fasta.gz", tmp / "reference")

    manifest = {
        "format": 1, "panel_version": params.pop("panel_version"), "built_with": f"mutation-hunter {__version__}",
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"), "inputs_md5": input_md5,
        "params": params, "counts": {"strains": len(strains), "variant_rows": len(var), "genes": int((genes["first"] >= 0).sum()),
                                     "clone_groups": int(strains["clone_group"].nunique()),
                                     "classes": var["class"].value_counts().to_dict()},
    }
    (tmp / "manifest.json").write_text(json.dumps(manifest, indent=2))
    files = sorted(p for p in tmp.rglob("*") if p.is_file())
    (tmp / "SHA256SUMS").write_text("".join(f"{file_hash(p, 'sha256')}  {p.relative_to(tmp)}\n" for p in files))
    shutil.rmtree(out, ignore_errors=True)
    tmp.rename(out)
```

- [ ] **Passo 4: Rodar e ver passar**

Run: `uv run --extra test --extra build pytest tests/test_build.py -v`
Esperado: 22 passed.

- [ ] **Passo 5: Commit**

```bash
git add src/mutation_hunter/build.py tests/test_build.py tests/conftest.py
git commit -m "Panel build: assemble the versioned store with manifest, kinship and clone groups"
```

---

### Tarefa 5: Leitor do pacote e download com cache

**Arquivos:**
- Criar: `src/mutation_hunter/store.py`, `tests/test_store.py`

**Interfaces:**
- Consome: o layout da Tarefa 4.
- Produz: `store.Panel(path=None)` com `.manifest`, `.strains`, `.genes`, `.variants`, `.gene(gene) -> (DataFrame, int8 variantes × cepas)`, `.kinship(loco=None)`, `.phenotypes()`, `.orfs()`, `.presence_absence()`, `.copy_number()`, `.verify()`; `store.fetch_store(url, digest, version) -> Path`; `store.default_path()` (`MH_PANEL` ou a versão publicada; `RELEASE` fica vazio até o Zenodo); `store.FORMAT = 1`.

- [ ] **Passo 1: Escrever os testes que falham** — `tests/test_store.py`:

```python
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
```

- [ ] **Passo 2: Rodar e ver falhar**

Run: `uv run --extra test --extra build pytest tests/test_store.py -v`
Esperado: `ImportError: cannot import name 'store'`.

- [ ] **Passo 3: Implementar** — `src/mutation_hunter/store.py`:

```python
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
```

- [ ] **Passo 4: Rodar e ver passar**

Run: `uv run --extra test --extra build pytest tests/test_store.py -v`
Esperado: 6 passed.

- [ ] **Passo 5: Commit**

```bash
git add src/mutation_hunter/store.py tests/test_store.py
git commit -m "Panel store reader with memory-mapped gene blocks and atomic cached download"
```

---

### Tarefa 6: CLI de etapas, Nextflow, ambiente e CI — e a construção real

**Arquivos:**
- Modificar: `src/mutation_hunter/build.py` (acrescentar `main`), `.github/workflows/tests.yml`, `.gitignore`
- Criar: `nextflow/build.nf`, `environment.yml`

**Interfaces:**
- Produz: `python -m mutation_hunter.build {check-inputs,reference,strains,genotypes,assemble} ...`; `nextflow run nextflow/build.nf --inputs <pasta> --panel_version <x.y.z> --outdir <saída>` → `<saída>/panel-<x.y.z>/`.

- [ ] **Passo 1: Implementar a CLI** — acrescentar a `src/mutation_hunter/build.py`:

```python


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="python -m mutation_hunter.build", description=__doc__)
    sub = p.add_subparsers(dest="step", required=True)
    s = sub.add_parser("check-inputs"); s.add_argument("inputs")
    s = sub.add_parser("reference"); s.add_argument("inputs"); s.add_argument("out")
    s = sub.add_parser("strains"); s.add_argument("inputs"); s.add_argument("samples"); s.add_argument("out")
    s = sub.add_parser("genotypes"); s.add_argument("bcf"); s.add_argument("strains"); s.add_argument("gff")
    s.add_argument("out_prefix"); s.add_argument("--min-gq", type=int, default=20); s.add_argument("--min-dp", type=int, default=5)
    s = sub.add_parser("assemble")
    for a in ("parts", "inputs", "reference", "strains", "out"):
        s.add_argument(a)
    s.add_argument("--panel-version", required=True); s.add_argument("--params", default="{}")
    a = p.parse_args(argv)
    if a.step == "check-inputs":
        print(json.dumps(check_inputs(a.inputs)))
    elif a.step == "reference":
        d = Path(a.inputs)
        reference(d / "Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz", d / "Saccharomyces_cerevisiae.R64-1-1.115.gff3.gz", a.out)
    elif a.step == "strains":
        d = Path(a.inputs)
        samples = Path(a.samples).read_text().split()
        strain_table(d / "41586_2018_30_MOESM3_ESM.xls", samples, d / "phenoMatrix_35ConditionsNormalizedByYPD.tab.gz").to_csv(a.out, sep="\t")
    elif a.step == "genotypes":
        genotypes(a.bcf, a.strains, a.gff, a.out_prefix, a.min_gq, a.min_dp)
    else:
        assemble(a.parts, a.inputs, a.reference, a.strains, a.out, check_inputs(a.inputs),
                 {"panel_version": a.panel_version, **json.loads(a.params)})
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Passo 2: Escrever `nextflow/build.nf`**

```groovy
// Build the panel store from the public 1011 genomes data (maintainers only).
// All logic lives in mutation_hunter.build and bin/annotate.sh. Needs, on PATH: bcftools >= 1.20 and
// mutation-hunter installed with the `build` extra (pip install 'mutation-hunter[build]').
nextflow.enable.dsl = 2

params.inputs = null          // folder with the pinned inputs (mutation_hunter.build.INPUTS); VCF indexed (.tbi)
params.outdir = 'panel'
params.panel_version = null
params.min_gq = 20
params.min_dp = 5

process CHECK {
    input:
    path inputs

    output:
    stdout

    script:
    "python -m mutation_hunter.build check-inputs $inputs"
}

process REFERENCE {
    input:
    path inputs
    val checked

    output:
    path 'reference'

    script:
    "python -m mutation_hunter.build reference $inputs reference"
}

process STRAINS {
    conda 'bioconda::bcftools=1.20'

    input:
    path inputs
    val checked

    output:
    path 'strains.tsv'

    script:
    """
    bcftools query -l $inputs/1011Matrix.gvcf.gz > samples.txt
    python -m mutation_hunter.build strains $inputs samples.txt strains.tsv
    """
}

process ANNOTATE {
    tag "$chrom"
    conda 'bioconda::bcftools=1.20'

    input:
    val chrom
    path inputs
    path reference

    output:
    tuple val(chrom), path("${chrom}.bcf"), path("${chrom}.bcf.csi")

    script:
    "annotate.sh $inputs/1011Matrix.gvcf.gz $chrom $reference/ref.fa $reference/ref.gff3 ${chrom}.bcf"
}

process GENOTYPES {
    tag "$chrom"
    memory '4 GB'

    input:
    tuple val(chrom), path(bcf), path(csi)
    path strains
    path reference

    output:
    path "${chrom}.*.{tsv.gz,npz}"

    script:
    """
    python -m mutation_hunter.build genotypes $bcf $strains $reference/ref.gff3 $chrom \\
        --min-gq $params.min_gq --min-dp $params.min_dp
    """
}

process ASSEMBLE {
    publishDir params.outdir, mode: 'copy'
    memory '16 GB'

    input:
    path parts
    path inputs
    path reference
    path strains
    val provenance

    output:
    path "panel-${params.panel_version}"

    script:
    """
    mkdir parts && mv ${parts} parts/
    python -m mutation_hunter.build assemble parts $inputs $reference $strains panel-${params.panel_version} \\
        --panel-version ${params.panel_version} \\
        --params '{"min_gq": $params.min_gq, "min_dp": $params.min_dp, "bcftools": "'"\$(bcftools --version | head -1)"'", ${provenance}}'
    """
}

workflow {
    if (!params.inputs || !params.panel_version) error "set --inputs <folder> and --panel_version <x.y.z>"
    inputs = file(params.inputs, checkIfExists: true)
    commit = "git -C $projectDir rev-parse HEAD".execute().text.trim()
    provenance = "\"nextflow\": \"${nextflow.version}\", \"pipeline_commit\": \"${commit}\""
    checked = CHECK(inputs)
    reference = REFERENCE(inputs, checked)
    strains = STRAINS(inputs, checked)
    chroms = Channel.of(1..16).map { "chromosome$it" }
    ANNOTATE(chroms, inputs, reference)
    GENOTYPES(ANNOTATE.out, strains, reference)
    ASSEMBLE(GENOTYPES.out.flatten().collect(), inputs, reference, strains, provenance)
}
```

- [ ] **Passo 3: Ambiente dos mantenedores** — `environment.yml`:

```yaml
# Maintainer environment: builds the panel and runs the tests that need bcftools.
#   conda env create -f environment.yml && conda activate mutation-hunter-build
name: mutation-hunter-build
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.11
  - bcftools=1.20
  - nextflow
  - pip
  - pip:
      - -e .[test,build]
```

- [ ] **Passo 4: CI** — substituir `.github/workflows/tests.yml` por:

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
  tools:
    runs-on: ubuntu-latest
    defaults:
      run:
        shell: bash -el {0}
    steps:
      - uses: actions/checkout@v4
      - uses: mamba-org/setup-micromamba@v2
        with:
          environment-file: environment.yml
      - run: pytest -v
```

E acrescentar ao `.gitignore`:
```
panel*/
*.bcf
*.csi
```

- [ ] **Passo 5: Suíte completa**

Run: `PATH=<pasta do bcftools>:$PATH uv run --extra test --extra build pytest -q`
Esperado: 60 passed.

- [ ] **Passo 6: Montar a pasta de entradas** (fora do git, ex.: `work/inputs/`)

Os 9 arquivos de `build.INPUTS`, com estes nomes exatos (links simbólicos servem):
- `1011Matrix.gvcf.gz` → `~/ic/vcf/1011Matrix.gvcf.gz`, e o índice ao lado: `bcftools index -t -o work/inputs/1011Matrix.gvcf.gz.tbi ~/ic/vcf/1011Matrix.gvcf.gz` (~90 s).
- Ensembl release 115: `https://ftp.ensembl.org/pub/release-115/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz` e `https://ftp.ensembl.org/pub/release-115/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.115.gff3.gz`.
- `http://1002genomes.u-strasbg.fr/files/` + `phenoMatrix_35ConditionsNormalizedByYPD.tab.gz`, `genesMatrix_PresenceAbsence.tab.gz`, `genesMatrix_CopyNumber.tab.gz`, `allORFs_pangenome.fasta.gz`, `1011GWASMatrix.tar.gz` (servidor lento: use `curl -C - --retry 10 -O`).
- `https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0030-5/MediaObjects/41586_2018_30_MOESM3_ESM.xls`.

Run: `uv run --extra build python -m mutation_hunter.build check-inputs work/inputs`
Esperado: um JSON com os 9 MD5, sem erro.

- [ ] **Passo 7: Construção real**

```bash
export PATH=$PWD/.venv/bin:<pasta do bcftools>:$PATH
nextflow run nextflow/build.nf --inputs work/inputs --panel_version 0.1.0 --outdir work/panel-out
```
Esperado (valores da construção de 15/10): ~4 min; o `manifest.json` registra bcftools, Nextflow e o commit; `manifest.json` com `strains` 1011, `variant_rows` 1.170.091, `genes` 6.516, `clone_groups` 811, classes `missense` 576.072, `synonymous` 540.642, `lof` 29.563, `inframe` 18.478, `other` 5.336; pacote ~556 MB em disco, ~138 MB em `.tgz`.

- [ ] **Passo 8: Conferir o pacote real**

```bash
uv run python -c "
import numpy as np
from mutation_hunter.store import Panel
p = Panel('work/panel-out/panel-0.1.0'); p.verify()
s = np.linalg.eigvalsh(p.kinship().astype(float))
print(round(float(np.diag(p.kinship()).mean()), 3), s.min() > -1e-6, p.gene('YNL079C')[1].shape)"
```
Esperado: `0.995 True (59, 1011)` — parentesco com diagonal média ~1 e sem autovalor negativo (o erro da calibração anterior não volta).

- [ ] **Passo 9: Commit**

```bash
git add src/mutation_hunter/build.py nextflow/build.nf environment.yml .github/workflows/tests.yml .gitignore
git commit -m "Panel build pipeline in Nextflow, maintainer environment and bcftools CI job"
```

---

### Tarefa 7: Verificação contra o proteoma público

**Arquivos:**
- Criar: `docs/build-check/compare_proteome.py`, `docs/build-check/README.md`

**Interfaces:**
- Consome: `store.Panel` (Tarefa 5) e o pacote real (Tarefa 6).

- [ ] **Passo 1: Escrever o script** — `docs/build-check/compare_proteome.py`:

```python
"""Manual check of a built panel against the public SNP-only proteome (1011proteome.tar.gz, extracted).

For every missense row on one chromosome, compare its carriers with the strains whose proteome sequence has the
alternative amino acid at that position. Disagreements come mostly from codons with two changes: bcftools csq -l
annotates each variant alone, while the proteome applies both (see README.md in this folder).

usage: python compare_proteome.py <panel folder> <Proteome_1011 folder> [chromosome1]
"""

import sys
from pathlib import Path

import numpy as np

from mutation_hunter.store import Panel


def load_gene(folder: Path, gene: str) -> dict[str, str]:
    seqs, strain = {}, None
    for line in (folder / f"{gene}.fasta").read_text().splitlines():
        if line.startswith(">"):
            strain = line[1:].split(f"_{gene}_")[0]
            seqs[strain] = ""
        else:
            seqs[strain] += line.strip()
    return seqs


def main(panel_dir, proteome_dir, chrom="chromosome1") -> None:
    p, prot = Panel(panel_dir), Path(proteome_dir)
    strains = np.array(p.strains.index)
    hom = het = n_rows = 0
    hom_ok = het_ok = prot_found = prot_total = 0
    for gene in p.genes.index[p.genes["chrom"] == chrom]:
        if not (prot / f"{gene}.fasta").exists():
            continue
        var, geno = p.gene(gene)
        seqs = load_gene(prot, gene)
        for i, (cls, aa) in enumerate(zip(var["class"], var["aa_change"])):
            ch = aa.split(">")
            if cls != "missense" or len(ch) != 2 or not ch[1][:-1].isdigit():
                continue
            pos, alt = int(ch[1][:-1]), ch[1][-1]
            in_prot = np.array([s in seqs for s in strains])
            prot_alt = np.array([s in seqs and len(seqs[s]) >= pos and seqs[s][pos - 1] == alt for s in strains])
            g = geno[i]
            n_rows += 1
            hom += ((g == 2) & in_prot).sum()
            hom_ok += ((g == 2) & prot_alt).sum()
            het += ((g == 1) & in_prot).sum()
            het_ok += ((g == 1) & prot_alt).sum()
            prot_total += (prot_alt & (g >= 0)).sum()
            prot_found += (prot_alt & (g > 0)).sum()
    print(f"{chrom}: {n_rows} single-codon missense rows")
    print(f"homozygous carriers that the proteome agrees with: {hom_ok}/{hom} = {hom_ok / hom:.3f}")
    print(f"heterozygous carriers that the proteome agrees with: {het_ok}/{het} = {het_ok / het:.3f}")
    print(f"proteome carriers (called in panel) that the panel also calls: {prot_found}/{prot_total} = {prot_found / prot_total:.3f}")


if __name__ == "__main__":
    main(*sys.argv[1:])
```

- [ ] **Passo 2: Baixar e rodar**

```bash
curl -C - --retry 10 -o work/1011proteome.tar.gz http://1002genomes.u-strasbg.fr/files/1011proteome.tar.gz
tar xzf work/1011proteome.tar.gz -C work/
uv run python docs/build-check/compare_proteome.py work/panel-out/panel-0.1.0 work/Proteome_1011 chromosome1
```
Esperado (15/10): 9.972 linhas; homozigotas que o proteoma confirma 0,940; heterozigotas 0,924; portadoras do proteoma que o painel também chama 0,987. Diferença grande indica regressão na construção.

- [ ] **Passo 3: Registrar** — `docs/build-check/README.md`:

```markdown
# Build check against the public proteome

`compare_proteome.py` compares the missense carriers of a built panel with `1011proteome.tar.gz`
(SNP-only proteins). Run on panel draft built 2026-10-15 (bcftools 1.20, min GQ 20, min DP 5):

| chromosome1, 9,972 single-codon missense rows | agreement |
|---|---|
| homozygous carriers the proteome agrees with | 204,428 / 217,364 = 0.940 |
| heterozygous carriers the proteome agrees with | 50,688 / 54,859 = 0.924 |
| proteome carriers the panel also calls | 255,116 / 258,466 = 0.987 |

Of the 12,936 homozygous disagreements, 9,855 (76%) are codons with a second change (the proteome shows a third
amino acid; `csq -l` annotates each variant alone, a documented limit) and 3,081 (24%) show the reference amino acid
in the proteome. Heterozygous calls agree far more than the 50% a random allele draw would give, so the proteome's
"random allele" is not uniform; we did not investigate further.
```

- [ ] **Passo 4: Commit**

```bash
git add docs/build-check
git commit -m "Check the built panel against the public SNP-only proteome"
```

---

## Fora deste plano

Testes de associação sobre o pacote novo e a nova calibração (subprojeto 2); máscara de regiões difíceis, apelidos RefSeq/GenBank e modo cepa externa (subprojeto 3); publicação no Zenodo (depois da autorização).
