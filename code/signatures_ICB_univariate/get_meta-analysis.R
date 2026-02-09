
suppressMessages(library(tidyverse))

source("./TLS_project/code/meta_analysis/meta-analysis.R")

signature_list <- list.files("./TLS_project/results/signatures_ICB_uni/")

for (signatureID in signature_list){
    prefix <- paste0("./TLS_project/results/signatures_ICB_uni/", signatureID, "/")
    cox_os_Rdata <- paste0(prefix, "COX_OS.Rdata")
    log_RData <- paste0(prefix, "Log_Response.Rdata")
    load(cox_os_Rdata)
    cox_os <- as.data.frame(cox_os)
    cox_os = cox_os[ as.numeric(cox_os$SE) >= 0 & !is.na( as.numeric(cox_os$Pval) ) , ]
    if( nrow( cox_os[ !is.na( cox_os[ , 1 ] ) , ] ) ){
      Get_Cox_Forestplot_ICB( data = cox_os , prefix = prefix , signatureID = signatureID ) 
    }
    
    load( log_RData )
    log_response <- as.data.frame(log_response)
    ## Meta-analysis of the Log Regression models (Response) Continous
    log_response = log_response[ as.numeric(log_response$SE) >= 0 & !is.na( as.numeric(log_response$Pval) ) , ]
    if( nrow( log_response[ !is.na( log_response[ , 1 ] ) , ] ) ){
      Get_LogReg_Forestplot( data = log_response , prefix= prefix , signatureID = signatureID ) 
    }
    
    
}




