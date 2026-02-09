
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(stringr))
setwd("/home/zhouchi/Project")
file_paths <- list.files("./TLS_project/data/TCGA_pheno_new/")




cancer_type <- list.files("./TLS_project/data/ICB_pheno_uni") %>% strsplit("\\.") %>% sapply("[", 1)
ICB_pheno <- list()


for (i in 1:length(cancer_type)){
    pheno <- read.table(paste0("./TLS_project/data/ICB_pheno_uni/", cancer_type[i], ".tsv"), header = T, sep = "\t")
    #pheno_filter <- subset(pheno, select = c("overall.survival..days.", "vital.status", "response_NR"))
    pheno_filter <- data.frame(sample_id=pheno$sample_id, os=pheno$overall.survival..days., status=pheno$vital.status, response=pheno$response_NR,age=pheno$age_start,sex=pheno$Gender)
    pheno_filter$status <- lapply(pheno_filter$status, function(x){ifelse(x == "Dead", 1, 0)})
    pheno_filter$response <- lapply(pheno_filter$response, function(x){ifelse(x == "R", 1, ifelse(x == "N", 0, NA))})
    #print(head(pheno_filter, n=5))
    pheno_filter$os <- as.numeric(pheno_filter$os)
    pheno_filter$status <- as.numeric(pheno_filter$status)
    pheno_filter$response <- as.numeric(pheno_filter$response)
    if (cancer_type[i]=="Nathanson_Melanoma_2017"){pheno_filter$sample_id <- paste0("X", pheno_filter$sample_id)}
    ICB_pheno[[cancer_type[i]]] <- pheno_filter
    
}

save(ICB_pheno, file = "./TLS_project/data/ICB_pheno_uni.Rdata")


cancer_type <- list.files("./TLS_project/data/ICB_pheno_uni") %>% strsplit("\\.") %>% sapply("[", 1)

ICB_expr <- list()

for (i in 1:length(cancer_type)){
    expr <- read.table(paste0("./TLS_project/data/ICB_expr/", cancer_type[i], ".tsv"), header = T, row.names = 1, sep = "\t") 
    ICB_expr[[cancer_type[i]]] <- expr
}

save(ICB_expr, file = "./TLS_project/data/ICB_expr_uni.Rdata")

