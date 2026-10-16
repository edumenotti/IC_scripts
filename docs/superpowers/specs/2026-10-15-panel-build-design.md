# Subprojeto 1: construção do painel — design

Data: 2026-10-15. Nome do pacote provisório (`mutation-hunter`).
Decisões que governam este documento: vault `decisions/2026-10-15-mutation-hunter-redesenho-vcf.md`.

## Objetivo

Transformar uma vez os dados públicos do painel dos 1011 genomas de *S. cerevisiae* (Peter et al. 2018) num
**pacote de dados do painel** pequeno (~150 MB), versionado e verificável. Os subprojetos 2 (motor de análise) e
3 (cepa externa) leem só esse pacote. O usuário nunca precisa do VCF de 5,4 GB.

Critério de sucesso: a construção roda de ponta a ponta a partir das entradas fixadas, produz um pacote cujo
manifesto descreve exatamente como foi feito, e as portadoras de missense raras no `chromosome1` coincidem com
as do branch `baseline/proteome`, exceto nos sítios heterozigotos.

## Contexto verificado (revisão de 15/10)

- `1011Matrix.gvcf.gz` **não é um gVCF**: é saída do GATK 3.3 GenotypeGVCFs com `includeNonVariantSites=false`.
  Só sítios variantes (1.754.866 registros), sem blocos de referência. Fora desses sítios, nada se sabe.
- Nenhum filtro de qualidade aplicado (FILTER é `.` em todos os registros; só QUAL ≥ 30). Há genótipos com
  1 leitura e GQ 3.
- Todos os genótipos foram chamados como diploides, mas o painel tem 132 haploides, 787 diploides e 92
  poliploides (Peter 2018, tabela S16). Haploides mostram "heterozigotos" que são artefato.
- `bcftools csq` com opções padrão para no primeiro heterozigoto sem fase; com `-s -` sem `-l`, combina
  variantes de cepas diferentes em consequências falsas. `-l -s -` (um registro por vez) dá classes limpas.
- Ensembl R64-1-1 release 115: os 16 cromossomos batem em tamanho com o VCF; 0 divergências de REF
  (`bcftools norm -c x`). A release 116 mudou o nome do GFF, por isso fixar a 115.
- As matrizes públicas de presença/ausência e número de cópias têm IDs de gene desfigurados pelo R
  (`X3051.YDR360W`, `YAL016C.A`); 1.348 genes da Ensembl não têm correspondência direta.
- A cópia local `~/ic/vcf/1011Matrix.gvcf.gz` tem o MD5 publicado (`42478e3e9dff4bd46993d82d8eab40d4`; o
  `md5.txt` a chama de `1011_matrix.gvcf.gz`). O servidor de origem entrega ~20 KB/s.

## Entradas (MD5 fixado no código)

| Entrada | Origem | Uso |
|---|---|---|
| `1011Matrix.gvcf.gz` | 1002genomes.u-strasbg.fr (ou cópia local) | variantes e genótipos |
| `Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz` (`5ae49d2d22cd977abdd6b054b94c544d`) | Ensembl release 115 | referência |
| `Saccharomyces_cerevisiae.R64-1-1.115.gff3.gz` (`17ed999c1853e03562f89bbcaefa80b6`) | Ensembl release 115 | anotação |
| `phenoMatrix_35ConditionsNormalizedByYPD.tab.gz` | 1002genomes | fenótipos |
| `genesMatrix_PresenceAbsence.tab.gz`, `genesMatrix_CopyNumber.tab.gz` | 1002genomes | genes acessórios e CNV |
| `allORFs_pangenome.fasta.gz` | 1002genomes | presença de ORFs numa cepa externa |
| `1011GWASMatrix.tar.gz` | 1002genomes | sítios comuns para posicionar cepa externa |
| `41586_2018_30_MOESM3_ESM.xls` (tabelas S1 e S16) | Nature, suplemento de Peter 2018 | clado, ploidia, zigosidade |

## Etapas

1. **Pacote de referência.** Remove a mitocôndria; renomeia `I…XVI` → `chromosome1…16` no FASTA e no GFF;
   indexa; gera `aliases.tsv` (`I`, `chrI`, `NC_001133.9`, `BK006935.2` → `chromosome1`, etc.). O subprojeto 3
   reutiliza esse pacote sem alterá-lo.
2. **Normalização.** `bcftools norm -f ref -c x -m -any -a` (separa multialélicos, alinha indels à esquerda,
   atomiza MNPs). Guarda o POS original.
3. **Consequência.** `bcftools csq -f ref -g gff -l -s - -v 0`. Limite documentado: duas variantes no mesmo
   códon da mesma cepa não são combinadas.
4. **Genótipos.** `bcftools query` extrai GT, GQ e DP. Viram ausentes: GQ < 20, DP < 5, e qualquer
   heterozigoto em cepa haploide (pela tabela S16). Codificação: −1 ausente, 0 referência, 1 het, 2 homozigoto
   alternativo. Não é dose aditiva: o motor usa "portadora" (> 0) e "homozigota" (= 2).
5. **Estatísticas por sítio.** AC, AN, n_het, n_hom, taxa de chamada, frequência sobre genótipos chamados, e
   `ref_is_minor` (a S288C carrega o alelo de menor frequência).
6. **Classes de consequência**, fixadas antes de qualquer resultado:
   - `missense`;
   - `lof`: stop_gained, frameshift, splice_acceptor, splice_donor, start_lost; exclui truncamento nos últimos
     5% do CDS;
   - `inframe`: inframe_insertion, inframe_deletion;
   - `synonymous`: só controle negativo na validação.
7. **Parentesco a partir do VCF.** SNPs bialélicos com frequência ≥ 5% e taxa de chamada ≥ 95%, genótipo
   codificado como portadora (0/1) para não depender da ploidia. Matriz de relação genômica (GRM) padronizada
   (VanRaden), ausentes imputados pela média; uma versão completa e uma LOCO por cromossomo. Grupos de clones: cepas cuja distância
   fica abaixo de um limiar (padrão: 1% da distância mediana; o valor e o número de grupos vão para o
   manifesto).
8. **IDs.** Limpa os nomes das matrizes públicas (`X\d+\.` sai, `.` vira `-`); liga ORFs do pangenoma aos nomes
   sistemáticos; gene sem correspondência fica marcado, nunca descartado em silêncio. A construção falha se a
   sobreposição de cepas entre VCF, fenótipo e matrizes for diferente da esperada.

## Pacote de dados `panel-<versão>/`

| Arquivo | Conteúdo |
|---|---|
| `manifest.json` | versão do formato e do pacote; MD5/SHA256 e origem de cada entrada; versões de bcftools, Nextflow e do código (commit); comandos exatos; limiares; contagens (sítios, cepas, genes, grupos de clones); data |
| `SHA256SUMS` | soma de cada arquivo |
| `strains.tsv` | id, ploidia, zigosidade, clado, grupo de clones, tem_fenótipo |
| `variants.tsv.gz` | id, chrom, pos, pos_original, ref, alt, gene, transcrito, consequência, classe, troca de aminoácido, posição relativa no CDS, AC, AN, n_het, n_hom, taxa de chamada, ref_is_minor |
| `genotypes/{data,indices,indptr}.npy` + `shape.json` | CSR variantes × cepas, int8, ordenado por gene; lido com `mmap` por bloco de gene |
| `genes.tsv` | gene, chrom, início, fim, primeira e última linha na matriz, ORF do pangenoma |
| `kinship.npy`, `kinship_loco/chromosomeN.npy` | GRM float32 1011 × 1011 |
| `phenotypes.tsv` | 971 × 35, como publicado |
| `presence_absence.npy`, `copy_number.npy`, `orfs.tsv` | cepas × ORFs com IDs limpos e mapa para nomes sistemáticos |
| `placement.npy`, `placement_sites.tsv` | doses nos SNPs comuns do `1011GWAS_matrix` (as 925 linhas de CNV, codificadas como cromossomo 17, ficam de fora) |
| `reference/` | FASTA e GFF renomeados, `aliases.tsv`, `allORFs_pangenome.fasta.gz` |

Sem dependências novas: arquivos `.npy`/`.npz` (numpy/scipy) e TSV (pandas). Nenhum `pickle`.

## Execução

- `nextflow run nextflow/build.nf --vcf <1011Matrix.gvcf.gz> --outdir panel/` a partir do repositório.
- Processos com `bcftools` declaram `conda 'bioconda::bcftools=1.20'` (ativado com `-with-conda`); sem ele, usam o
  `bcftools` do `PATH`. Ambiente completo dos mantenedores em `environment.yml`. Passos numéricos ficam em `src/mutation_hunter/build.py`, chamados pelo Nextflow.
- A construção é feita pelos mantenedores; o usuário comum não a executa.

## Distribuição

- Versão do pacote de dados independente da versão do código (`panel-1.0.0`).
- Zenodo, um DOI por versão, **só depois de autorização** do grupo Schacherer/Liti. Até lá, uso local.
- O pacote Python fixa o DOI da versão e o SHA256 do arquivo. `MH_PANEL` aponta para um pacote local.
- Cache em `~/.cache/mutation_hunter/<versão>/`; escrita em diretório temporário, renomeação atômica e marcador
  `.complete` gravado por último.

## Testes

1. **Unitários, sem ferramentas externas:** gerador de pacote sintético em Python (usado também pelo
   subprojeto 2); testa leitura por bloco de gene, IDs, manifesto e cache.
2. **Construção com bcftools** (pulados sem bcftools; na CI, job com micromamba): VCF escrito à mão, 8 cepas,
   cromossomo falso de 2 kb, GFF com 2 genes (um na fita reversa). Casos: multialélico, indel, `./.`, het em
   haploide, GQ baixo, stop nos últimos 5%, splice, sinônima. Saída comparada com um TSV esperado.
3. **Dados reais, manual:** constrói só o `chromosome1` e compara portadoras de missense raras com o branch
   `baseline/proteome`. Contagens e divergências registradas em `docs/build-check/README.md`.

## Fora do escopo

Máscara de regiões difíceis (subtelômeros, transposons) e callable mask para cepa externa → subprojeto 3.
Profundidade de alelos (AD) de poliploides → não usada pela codificação aprovada. Mitocôndria. Testes de
associação → subprojeto 2.
