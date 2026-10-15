# Mutation Hunter 2.0: design

Data: 2026-10-15. Nome provisório; a troca é um find-and-replace de `mutation-hunter` / `mutation_hunter`.

## O que é

Um pacote Python, de código aberto, que cruza variantes raras de aminoácido do painel dos 1011 genomas de
*S. cerevisiae* (Peter et al. 2018) com fenótipos. Modelo de referência: o CaeNDR/NemaScan (*C. elegans*) e o
DGRPool (*Drosophila*).

Faz duas perguntas, com o mesmo núcleo:

- **PheWAS**: parte de uma variante (`YNL079C:A123V`) ou de um gene (`YNL079C`) e testa os 35 fenótipos.
- **GWAS**: parte de um fenótipo (uma das 35 condições ou um TSV do usuário) e testa todos os genes ou
  todas as variantes raras.

## Decisões do Eduardo

| Data | Decisão |
|---|---|
| 04/10 | Público: publicação e comunidade. Pacote instalável, CLI, testes. |
| 04/10 | Só *S. cerevisiae*, painel 1011 e matriz de 35 condições. |
| 04/10 | Lógica só no pacote Python; Nextflow é casca fina. |
| 15/10 | Referência de direção: CaeNDR (portal + NemaScan) e DGRPool (PheWAS). |
| 15/10 | Dois pontos de entrada: variante/gene (PheWAS) e fenótipo (GWAS). |
| 15/10 | Formato em camadas: (1) pacote → (2) resultados pré-calculados com DOI → (3) atlas estático, opcional. Este design cobre só a camada 1. |
| 15/10 | Código aberto, organizado e reprodutível. Sem coautoria do Claude nos commits. |
| 15/10 | Após a calibração: modelo misto + 10 PCs de ORFs + INT; GWAS por gene; raridade ≤ 2%. |

## Arquitetura

```
src/mutation_hunter/
  data.py      baixa, confere MD5 e lê os 4 arquivos públicos; cache em ~/.cache/mutation_hunter
  variants.py  matriz cepas × posições; variantes raras; portadoras de uma variante
  stats.py     parentesco, PCs, INT, modelo misto (EMMAX), teste hipergeométrico, BH
  analysis.py  phewas() e gwas(): juntam os três módulos acima
  report.py    TSV + JSON de proveniência + HTML de arquivo único
  cli.py       mutation-hunter fetch | phewas | gwas | fdr
nextflow/main.nf   FETCH → GWAS por fenótipo em paralelo → FDR global
tests/             painel falso de 40 cepas, sem internet
docs/calibration/  evidência e scripts da escolha do modelo
```

Dependências: numpy, pandas, scipy, matplotlib. Python ≥ 3.11. Nada de BLAST, MAFFT ou R.

### Dados (todos de http://1002genomes.u-strasbg.fr/files/, MD5 fixado no código)

| Arquivo | Uso |
|---|---|
| `1011proteome.tar.gz` | Um FASTA por gene (6.575), já alinhado. Cepa ausente = gene ausente. |
| `phenoMatrix_35ConditionsNormalizedByYPD.tab.gz` | 971 cepas × 35 condições, sem NA. |
| `1011DistanceMatrixBasedOnSNPs.tab.gz` | 1011 × 1011; base da matriz de parentesco (efeito aleatório). |
| `1011DistanceMatrixBasedOnORFs.tab.gz` | 1011 × 1011, genoma acessório; 10 PCs entram como covariáveis. |

IDs: o código da cepa é o prefixo do cabeçalho antes de `_<GENE>_`, o que preserva `SACE_YBA` (93 cepas;
era o bug do `seq_id[:3]`). Cruzamento real: 970 cepas têm fenótipo e distância. `ARC` não tem distância.

### Saídas

Cada análise grava `<prefixo>.tsv` (o produto), `<prefixo>.json` (parâmetros, versões, MD5 dos dados,
λGC) e `<prefixo>.html` (gráficos embutidos, abre sem servidor).

Colunas do TSV: `phenotype gene test n n_carriers beta se p p_hypergeom q`.

- `beta`: efeito de ser portadora, em desvios-padrão do fenótipo transformado (INT).
- `p`: modelo misto. `p_hypergeom`: o teste original da IC, sem ajuste, mantido como comparação.
- `q`: Benjamini-Hochberg sobre todos os testes da análise. No Nextflow, `fdr` recalcula sobre tudo.

## Método (decidido em 15/10 pela calibração)

1. **Alelo de referência = alelo majoritário** no painel, porque o proteoma não traz a S288C à parte.
2. **Variante rara**: frequência ≤ 2% entre as cepas que têm o gene e ≥ 3 portadoras.
3. **Burden por gene** (CAST): portadora = tem qualquer variante rara do gene. É o padrão do GWAS.
4. **Modelo misto EMMAX**: parentesco K de SNPs (−½·J·D²·J, MDS clássico) como efeito aleatório; 10 PCs do
   parentesco de ORFs como covariáveis fixas; razão de variâncias estimada uma vez por fenótipo.
5. **Fenótipo transformado por postos (INT)** antes do teste.
6. **Cepa sem o gene**: genótipo ausente, imputado pela média (não entra como referência).
7. **Heterozigose**: aceita o alelo sorteado pelo proteoma público.
8. **Teste original**: hipergeométrico no fenótipo bruto, extremo = média + 1,5 DP, reportado ao lado.

## Calibração (15/10)

Detalhes, tabelas e scripts em `docs/calibration/README.md`. Resumo:

- Só o modelo misto ficou calibrado em fenótipos nulos simulados (λ 0,91–1,00). OLS: 4,1–10,3; PCs de SNPs: 2,1–3,2.
- O GWAS por variante quase não tem poder (0–8% para efeito de 1 DP). Por gene: 85–95%.
- Raridade: o λ real sobe com o limiar (1,64 a 1% → 1,93 a 10%); o poder estabiliza a partir de 2%.
- O parentesco de SNPs não vê a estrutura do genoma acessório (λ 1,73, cauda 61× em nulos simulados por ORFs).
  10 PCs de ORFs corrigem (cauda 1,4–2,2×) sem perder poder (0,98).
- O INT conserta a cauda inflada por fenótipos assimétricos com poucas portadoras.

Efeito nos dados reais (cobre, YPDCUSO410MM, por gene):

| Versão | λGC | q < 0,05 |
|---|---|---|
| Burden 10%, 5 PCs de SNPs | 19,6 | 4.302 |
| Burden 10%, modelo misto | 4,25 | 980 |
| **Final: 2%, modelo misto + 10 PCs de ORFs + INT** | **1,41** | **2** (YOL003C, YPL028W) |

PheWAS do *TPM1* (YNL079C): os p de 10⁻¹³ da primeira versão sumiram (λ 0,88, nada significativo).

## Decisões em aberto (do Eduardo)

1. O λ real residual (~1,4) é poligenia ou estrutura não capturada? Não dá para separar com estes testes.
2. Heterozigose: manter o sorteio ou usar o CDS com códigos de ambiguidade.
3. Licença (proposta: MIT) e idioma da documentação (proposta: inglês).
4. Impacto do bug do `SACE_` nos resultados da IC.

## Fora do escopo desta versão

Relatório de cepa nova, priorização por efeito (ESM/SIFT), CNVs e genes acessórios, frequência por clado
(precisa da tabela S1 do artigo), camadas 2 e 3, bioconda, contêiner.
