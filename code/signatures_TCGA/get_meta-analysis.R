
source("./TLS_project/code/meta_analysis/meta-analysis.R")

signature_list <- list.files("./TLS_project/results/signatures_TCGA/")

hot_tumor=c('BLCA', 'BRCA', 'CESC', 'COAD', 'READ', 'HNSC', 'KIRC', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PRAD', 'SKCM', 'STAD', 'UCEC','LUNG')

cold_tumor=c('ACC', 'CHOL', 'DLBC', 'ESCA', 'KICH', 'KIRP', 'TGCT', 'THCA', 'THYM', 'UCS', 'UVM', 'GBM', 'LGG', 'PCPG', 'SARC')

print(signature_list)
## All 
for (signatureID in signature_list){
    prefix <- paste0("./TLS_project/results/signatures_TCGA/", signatureID, "/")
    print(prefix)
    ###OS
    cox_os_Rdata <- paste0(prefix, "COX_OS.Rdata")
    load(cox_os_Rdata)
    cox_os <- as.data.frame(cox_os)
    #hot
    cox_os_hot=cox_os[which(cox_os$study %in% c(hot_tumor)),]
    cox_os_hot = cox_os_hot[ cox_os_hot$SE <= 10 & !is.na( cox_os_hot$Pval ) , ]
    if( nrow( cox_os_hot[ !is.na( cox_os_hot[ , 1 ] ) , ] ) ){
      # cancer <- Get_Cancer( cancer=cox_os )
      # seq <- Get_Seq( seq=cox_os )
      Get_Cox_Forestplot( data = cox_os_hot , prefix = prefix , signatureID = signatureID, TIL_type = "hot_tumor") 
    }
    #cold
    cox_os_cold=cox_os[which(cox_os$study %in% c(cold_tumor)),]
    cox_os_cold = cox_os_cold[ cox_os_cold$SE <= 10 & !is.na( cox_os_cold$Pval ) , ]
    if( nrow( cox_os_cold[ !is.na( cox_os_cold[ , 1 ] ) , ] ) ){
      # cancer <- Get_Cancer( cancer=cox_os )
      # seq <- Get_Seq( seq=cox_os )
      Get_Cox_Forestplot( data = cox_os_cold , prefix = prefix , signatureID = signatureID, TIL_type = "cold_tumor") 
    }



    ###PFI
    cox_pfi_Rdata <- paste0(prefix, "COX_PFI.Rdata")
    load(cox_pfi_Rdata)
    cox_pfi <- as.data.frame(cox_pfi)
    ##hot
    cox_pfi_hot=cox_pfi[which(cox_pfi$study %in% c(hot_tumor)),]
    cox_pfi_hot = cox_pfi_hot[ cox_pfi_hot$SE <= 10 & !is.na( cox_pfi_hot$Pval ) , ]
    if( nrow( cox_pfi[ !is.na( cox_os[ , 1 ] ) , ] ) ){
      # cancer <- Get_Cancer( cancer=cox_os )
      # seq <- Get_Seq( seq=cox_os )
      Get_Cox_Forestplot_Pfi( data = cox_pfi_hot , prefix = prefix , signatureID = signatureID,  TIL_type = "hot_tumor") 
    }
    ##cold
    cox_pfi_cold=cox_pfi[which(cox_pfi$study %in% c(cold_tumor)),]
    cox_pfi_cold = cox_pfi_cold[ cox_pfi_cold$SE <= 10 & !is.na( cox_pfi_cold$Pval ) , ]
    if( nrow( cox_pfi[ !is.na( cox_os[ , 1 ] ) , ] ) ){
      # cancer <- Get_Cancer( cancer=cox_os )
      # seq <- Get_Seq( seq=cox_os )
      Get_Cox_Forestplot_Pfi( data = cox_pfi_cold , prefix = prefix , signatureID = signatureID,  TIL_type = "cold_tumor") 
    }
    
    
}

