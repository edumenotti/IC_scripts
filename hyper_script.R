library(dplyr)
library(readr)
library(readxl)
library(stats) # Para a função phyper()

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 1) {
  stop("No phenotype specified as argument.\nUsage: Rscript hyper_script.R PhenotypeName")
}

phenotype <- args[1]

if(length(args) >= 2) {
  csv_folder <- args[2]
} else {
  csv_folder <- dirname(normalizePath(rstudioapi::getActiveDocumentContext()$path))
}

# Listar todos os arquivos CSV no diretório
csv_files <- list.files(path = csv_folder, pattern = "_mutations.csv$", full.names = TRUE)

# Carregar dados fenotípicos
phenotypic_data <- read_excel("phenoMatrix_35ConditionsNormalizedByYPD.xlsx")
colnames(phenotypic_data)[1] <- "cepas"

phenotype_mean <- mean(phenotypic_data[[phenotype]], na.rm = TRUE)
phenotype_sd <- sd(phenotypic_data[[phenotype]], na.rm = TRUE)

# Calcular o valor de corte
threshold <- phenotype_mean + (1.5 * phenotype_sd)

# Inicializar um vetor para armazenar os resultados
filtered_results <- data.frame(Gene = character(), Mutation = character(), stringsAsFactors = FALSE)

# Função para processar cada arquivo
hyper <- function(file_path) {
  dados_mutacoes <- read.csv(file_path)
  
  if(nrow(dados_mutacoes) <= 0) {
    return(NULL) # Pula o arquivo se ele tiver apenas o cabeçalho ou estiver vazio
  }
  
  temp_list <- list()
  
  for(i in 1:nrow(dados_mutacoes)) {
    mutation_strains <- unlist(strsplit(dados_mutacoes$Strains[i], ";"))
    phenotypic_data_filtered <- phenotypic_data[phenotypic_data$cepas %in% mutation_strains, ]
    
    M <- nrow(phenotypic_data)
    N <- sum(phenotypic_data[[phenotype]] > threshold)
    n <- length(mutation_strains)
    k <- sum(phenotypic_data_filtered[[phenotype]] > threshold)
    
    if(M > 0 && N <= M && n <= M && k <= N && n > 0) {
      p_value <- phyper(k-1, N, M-N, n, lower.tail=FALSE)
      if(p_value < 0.05) {
        gene_name <- tools::file_path_sans_ext(basename(file_path))
        gene_name <- strsplit(gene_name, "_mutations")[[1]][1]
        temp_list[[i]] <- data.frame(Gene = gene_name, Mutation = dados_mutacoes$Mutation[i], P_Value = p_value)
      }
    }
  }
  
  if(length(temp_list) > 0) {
    return(do.call(rbind, temp_list))
  } else {
    return(NULL)
  }
}

# Aplicar a função a cada arquivo e armazenar os resultados
for(file in csv_files) {
  results <- hyper(file)
  if(!is.null(results)) {
    filtered_results <- rbind(filtered_results, results[, c("Gene", "Mutation")])
  }
}

# Remover duplicatas, se houver
filtered_results <- unique(filtered_results)

output_csv_path <- file.path(getwd(), "filtered_genes_mutations.csv")

# Verificar se o arquivo já existe e decidir entre escrever ou fazer append
if(file.exists(output_csv_path)) {
  write_csv(filtered_results, output_csv_path, append = TRUE)
} else {
  write_csv(filtered_results, output_csv_path)
}
# Imprimir os resultados no formato solicitado
print(filtered_results)
