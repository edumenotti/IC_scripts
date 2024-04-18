Conjunto de scripts em Python, R e Groovy desenvolvidos durante o projeto de iniciação científica do aluno Eduardo Menoti Silva sob financiamento da FAPESP, processo 2022/13111-4. O intuito do projeto foi desenvolver uma ferramenta capaz de identificar mutações raras e verificar seu provável peso na expressão de fenótipos de interesse em leveduras, utilizando resistência a estresse oxidativo como estudo de caso.

### Validade e usabilidade

É importante ressaltar que o script em sua forma atual realiza as análises de cálculo de mutações com base na cepa industrial de *Saccharomyces cerevisiae* Pedra2 quando os arquivos de input estão no formato .fasta. Para alterar, é necessário cria um novo banco de dados BLAST utilizando os dados que se quer ter como base através do comando `makeblastdb -in {dados de interesse} -dbtype {proteínas ou nucl} -out {nome base dos arquivos de saída}`. Caso os arquivos de entrada já estejam no formato XML, o script funcionará normalmente.

### Hyper_script

Está disponível também o código utilizado para análise por distribuição hipergeométrica da significância estatística das mutações encontradas. 

# Arquivos CLUSTAL
Para o bom funcionamento do script desenvolvido, é fundamental o download dos dados de alinhamento CLUSTAL disponibilizados no link:
https://drive.google.com/file/d/1j4gbMhXMoOlZnTipuzuOo7MyvspAB0GX/view?usp=sharing

A pasta contendo os arquivos de alinhamento deve ser corretamente referenciada no parâmetro `--alignment_dir` do script Python.

# Pipeline Nextflow
O pipeline MutationHunter compõe o passo a passo necessário para partir de arquiovs GenBank contendo uma série de genes a serem analisados até a parte final de análise da relevância das mutações por análise de distribuição hipergeométrica. 
Para executar o pipeline, atente-se aos parâmetros necessários

# Script GenBank_extract
Para o bom funcionamento do pipeline, é necessário também a presença de um segundo script .py responsável por extrair a sequência .fa de todos os genes contidos nos arquivos GenBank providenciados.
