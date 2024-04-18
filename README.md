Conjunto de scripts desenvolvidos durante o projeto de iniciação científica do aluno Eduardo Menoti Silva sob financiamento da FAPESP, processo 2022/13111-4.

Para o bom funcionamento do script desenvolvido, é fundamental o download dos dados de alinhamento CLUSTAL disponibilizados no link:
https://drive.google.com/file/d/18zFyeWcrTIstM4N-KX8Dl3aR8p_RPmMY/view?usp=drive_link
#

# Pipeline Nextflow
O pipeline MutationHunter compõe o passo a passo necessário para partir de arquiovs GenBank contendo uma série de genes a serem analisados até a parte final de análise da relevância das mutações por análise de distribuição hipergeométrica. 
Para executar o pipeline, atente-se aos parâmetros necessários

# Script GenBank_extract
Para o bom funcionamento do pipeline, é necessário também a presença de um segundo script .py responsável por extrair a sequência .fa de todos os genes contidos nos arquivos GenBank providenciados.
