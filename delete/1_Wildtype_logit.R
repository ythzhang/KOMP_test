
#######################################################################################################
#1: SD in metabolomics data (Wildtype mice)
#load raw data file from metabolomics assays
#subset into categorical and continuous datasets

library(wcmc)
library(data.table)
library(RNOmni)

setwd("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
Null_Metabolomics <- wcmc::read_data("Supplementary Data 1.xlsx")

Null_Metabolomics_p <- Null_Metabolomics$p
Null_Metabolomics_f <- Null_Metabolomics$f
Null_Metabolomics_e <- Null_Metabolomics$e_matrix


## 1a: wildtype mice - sex
### deal with missing values: filter missing value > 70%, replace other missing values.

num_missing_male <- apply(Null_Metabolomics_e, 1, function(x){
  sum(is.na(x[Null_Metabolomics_p$Gender %in% "Male"]))
})
num_missing_female <- apply(Null_Metabolomics_e, 1, function(x){
  sum(is.na(x[Null_Metabolomics_p$Gender %in% "Female"]))
})

names(num_missing_male) = names(num_missing_female) = Null_Metabolomics_f$label
missing_index <- ((num_missing_male/ sum(Null_Metabolomics_p$Gender %in% "Male")) >= 0.7) & ((num_missing_female/ sum(Null_Metabolomics_p$Gender %in% "Female")) >= 0.7)
Null_Metabolomics_e <- Null_Metabolomics_e[!missing_index, ]
Null_Metabolomics_f <- Null_Metabolomics_f[!missing_index, ]

for(i in 1 : nrow(Null_Metabolomics_e)){
  Null_Metabolomics_e[i,is.na(Null_Metabolomics_e[i, ])] <- 0.5 * min(Null_Metabolomics_e[i, ], na.rm = TRUE)
}

rownames(Null_Metabolomics_e) <- Null_Metabolomics_f$label

## statistics
########################################################################

test1_Norma = Normality_F = Normality_M = p_val1 = fc1 = p_val_adj1 = F_Count = M_Count = Coef1 = c()

for (i in 1 : nrow(Null_Metabolomics_e)){

     Null_Metabol <- data.table(Intensity = Null_Metabolomics_e[i, ], Gender = Null_Metabolomics_p$Gender)
     Null_Metabol <- Null_Metabol[!is.na(Null_Metabol$Intensity), ]
 
     test1_Norma[i] <- shapiro.test(Null_Metabol$Intensity)$p.value
     Null_Metabol$trans <- rankNorm(Null_Metabol$Intensity)

     Null_Metabolomics_f$CompoundName[i]  
     Null_Metabolomics_f$label[i] 
  
  
     F_Count[i] <- sum((Null_Metabol$Gender %in% "Female"))      
     M_Count[i] <- sum(Null_Metabol$Gender %in% "Male")
  
  
     Null_Metabol_M <-  Null_Metabol[Null_Metabol$Gender %in% "Male"]
     Null_Metabol_F <-  Null_Metabol[Null_Metabol$Gender %in% "Female"]
  
     Normality_M[i] <- shapiro.test(Null_Metabol_M$trans)$p.value
     
      if( length(unique(Null_Metabol_F)) == 1 ){
        Normality_F[i] <- NA
    
        }else{
           Normality_F[i] <- shapiro.test(Null_Metabol_F$trans)$p.value 
    
           }
  
  
     Null_e_M <- Null_Metabol[Null_Metabol$Gender %in% "Male", ]
     Null_e_F <- Null_Metabol[Null_Metabol$Gender %in% "Female", ]
  
  
     lm_result <- summary(gls(Null_Metabol$trans ~ Null_Metabol$Gender, na.action = 'na.exclude'))
     p_val1[i] <- summary(gls(Null_Metabol$trans ~ Null_Metabol$Gender, na.action = 'na.exclude'))$coefficients[2, 4]
     Coef1[i] <- summary(gls(Null_Metabol$trans ~ Null_Metabol$Gender, na.action = 'na.exclude'))$coefficients[2, 1]
     fc1[i] <-  mean(Null_Metabol$Intensity[Null_Metabol$Gender %in% "Male"])/ mean(Null_Metabol$Intensity[Null_Metabol$Gender %in% "Female"])
     
     
     }

    p_val_adj1 <- p.adjust(p_val1, method = "fdr")

fwrite(data.table(label = Null_Metabolomics_f$label, Metbaolites = Null_Metabolomics_f$CompoundName,Platform= Null_Metabolomics_f$Assay, F_Count = F_Count, M_Count = M_Count, Normality_before = test1_Norma, Normality_M = Normality_M, Normality_F = Normality_F, p_value = p_val1, Coefficient = Coef1, fold_change = fc1, adjusted_p_value = p_val_adj1),"xxxxxxxxxxxxxxxxxxxxxxxxx")


#######################################################################################################
##1b: wildtype mice - sex adjust by body weight
#load bodyweight data from phenotype results 

Phenotype_bw <- wcmc::read_data("xxxxxxxxxxxxxxxxxxxxxxxxx.xlsx")
Phenotype_bw_p <- Phenotype_bw$p
Phenotype_bw_f <- Phenotype_bw$f
Phenotype_bw_e <- Phenotype_bw$e_matrix

body_weight_all <- Phenotype_bw_e[,Phenotype_bw_p[["phenotype"]] %in% "Body weight"]
names(body_weight_all) <- Phenotype_bw_f$label
body_weight <- body_weight_all[Null_Metabolomics_p$label]
names(body_weight) <- Null_Metabolomics_p$label

## statistics
p_val2 = fc2 = p_val_adj2 = F_Count = M_Count = Coef2 =  c()

 for (i in 1 : nrow(Null_Metabolomics_e)){

      Null_Metabolomics_e_trans <- rankNorm(Null_Metabolomics_e[i, ])

      Null_Metabolomics_e_Stats <- data.table(Null_Metabolomics_p$Gender, Null_Metabolomics_e[i, ], Null_Metabolomics_e_trans, body_weight)
      # Null_Metabolomics_f$CompoundName[i]  
      # Null_Metabolomics_f$label[i]  


      colnames(Null_Metabolomics_e_Stats) <- c("Gender", "MetaboVAlue_raw", "MetaboVAlue_trans", "Weight")
  
      Null_Metabolomics_e_Stats <- Null_Metabolomics_e_Stats[!is.na(Null_Metabolomics_e_Stats$Weight), ]
      Null_Metabolomics_e_Stats$weight_trans <- rankNorm(Null_Metabolomics_e_Stats$Weight)
  
      F_Count[i] <- sum((Null_Metabolomics_e_Stats$Gender %in% "Female"))      
      M_Count[i] <- sum(Null_Metabolomics_e_Stats$Gender %in% "Male")
      # body_weight_stats<-body_weight[names(body_weight)%in%Null_Metabolomics_p_New$label]

      p_val2[i] <- summary(gls(Null_Metabolomics_e_Stats$MetaboVAlue_trans ~ Null_Metabolomics_e_Stats$Gender + Null_Metabolomics_e_Stats$weight_trans, na.action = 'na.exclude'))$coefficients[2, 4]
      Coef2[i] <- summary(gls(Null_Metabolomics_e_Stats$MetaboVAlue_trans ~ Null_Metabolomics_e_Stats$Gender + Null_Metabolomics_e_Stats$weight_trans, na.action = 'na.exclude'))$coefficients[2, 1]
  
      fc2 [i] <-  mean(Null_Metabolomics_e_Stats$MetaboVAlue_raw[Null_Metabolomics_e_Stats$Gender %in% "Male"])/ mean(Null_Metabolomics_e_Stats$MetaboVAlue_raw[Null_Metabolomics_e_Stats$Gender %in% "Female"])

  
  
  
    }

p_val_adj2 <- p.adjust(p_val2, method = "fdr")

fwrite(data.table(label = Null_Metabolomics_f$label, Metbaolites = Null_Metabolomics_f$CompoundName, Platform = Null_Metabolomics_f$Assay, F_Count = F_Count, M_Count = M_Count, p_value = p_val2, Coefficient = Coef2, fold_change = fc2, adjusted_p_value = p_val_adj2), "XXXXXXXXXXXXXXXXXXXXXXXXXX")




#######################################################################################################
#2: SD in phenotype data (Wildtype mice)
#load raw data file from metabolomics assays

## 2a: wildtype mice phenotype - sex
setwd("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
Raw_Phenotype <- wcmc::read_data("Supplementary Data 6.xlsx")

Raw_Phenotype_UCDavis_p <- Raw_Phenotype_UCDavis$p
Raw_Phenotype_UCDavis_f <- Raw_Phenotype_UCDavis$f
Raw_Phenotype_UCDavis_e <- Raw_Phenotype_UCDavis$e_cat_matrix
# nrow(Raw_Phenotype_UCDavis_e)

Raw_Phenotype_UCDavis_f$label_procedure <- sapply(strsplit(Raw_Phenotype_UCDavis_f$label,"___"), function(x){x[1]})
Raw_Phenotype_UCDavis_f$label_phenotype <- sapply(strsplit(Raw_Phenotype_UCDavis_f$label,"___"), function(x){x[2]})
# rownames(Raw_Phenotype_UCDavis_e) = Raw_Phenotype_UCDavis_f$label


# combine the phenotypes because one phenotype may have multiple procedures.
# length(Raw_Phenotype_UCDavis_f$label_phenotype)
unique_phenotype <- unique(Raw_Phenotype_UCDavis_f$label_phenotype)

for(u in 1 : length(unique_phenotype)){
  if(sum(Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u]) > 1){
    if(length(table(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],])) > 1){
      duplicated_value_index = which(apply(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u], 1 : 2], 2, function(x){
        if(length(unique(x)) == 1){
          if(!unique(x) == "NA"){
            return(TRUE)
          }else{
            return(FALSE)
          }
        }else{
          return(FALSE)
        }
      }))
      if(length(duplicated_value_index) > 0){
        print(unique_phenotype[u])
        print(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u], 1 : 10])
        Sys.sleep(20)
      }
    }
  }
}


Raw_Phenotype_UCDavis_e_merge <- matrix(NA, nrow = length(unique(Raw_Phenotype_UCDavis_f$label_phenotype)), ncol = ncol(Raw_Phenotype_UCDavis_e))
rownames(Raw_Phenotype_UCDavis_e_merge) <- unique(Raw_Phenotype_UCDavis_f$label_phenotype)


i <- 1
for(u in 1 : length(unique_phenotype)){
  
  if(sum(Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u]) > 1){
    # stop(u)
    for(j in 1 : ncol(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u], ])){
      the_j_th_col = Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u], j]
      if(!length(unique(the_j_th_col)) == 1){
        Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u], ][the_j_th_col == "NA", j] = Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],][!the_j_th_col == "NA", j] 
      } 
    }
    Raw_Phenotype_UCDavis_e_merge[i, ] <- unique(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u], ])
    i <- i + 1
  }else{
    Raw_Phenotype_UCDavis_e_merge[i, ] <- Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u], ]
    i <- i + 1
  }
}


Raw_Phenotype_UCDavis_e_merge[Raw_Phenotype_UCDavis_e_merge == "NA"] <- NA


Raw_Phenotype_UCDavis_e_null <- Raw_Phenotype_UCDavis_e_merge[,Raw_Phenotype_UCDavis_p$Genotype %in% "null"]
Raw_Phenotype_UCDavis_p_null <- Raw_Phenotype_UCDavis_p[Raw_Phenotype_UCDavis_p$Genotype %in% "null", ]
Raw_Phenotype_UCDavis_e_null <- as.matrix(Raw_Phenotype_UCDavis_e_null)

Raw_Phenotype_UCDavis_e_null <- apply(Raw_Phenotype_UCDavis_e_null, 2, as.numeric)
rownames(Raw_Phenotype_UCDavis_e_null) <- rownames(Raw_Phenotype_UCDavis_e_merge)

### deal with missing values: filter missing value > 70%, replace other missing values.
num_missing_male <- apply(Raw_Phenotype_UCDavis_e_null, 1, function(x){
  sum(is.na(x[Raw_Phenotype_UCDavis_p_null$Gender %in% "Male"]))
})
num_missing_female <- apply(Raw_Phenotype_UCDavis_e_null, 1, function(x){
  sum(is.na(x[Raw_Phenotype_UCDavis_p_null$Gender %in% "Female"]))
})


missing_index <- ((num_missing_male/ sum(Raw_Phenotype_UCDavis_p_null$Gender %in% "Male")) >= 0.7) & ((num_missing_female/ sum(Raw_Phenotype_UCDavis_p_null$Gender %in% "Female")) >= 0.7)
Phenotype_UCDavis_e_null <- Raw_Phenotype_UCDavis_e_null[!missing_index, ]
rownames(Phenotype_UCDavis_e_null) <- rownames(Raw_Phenotype_UCDavis_e_merge)[!missing_index]


p_val3 = fc3 = p_val_adj3 = Coef3 = F_Count = M_Count = continuous_index = c()

for (i in 1 : nrow(Phenotype_UCDavis_e_null)){
    
  
  #subset into categorical and continuous datasets
  continuous <- sum(is.na(as.numeric(Phenotype_UCDavis_e_null[i, !is.na(Phenotype_UCDavis_e_null[i, ])]))) == 0
  
  
  
  if(continuous){
      
    continuous_index[i] <- TRUE
  
    Phenotype_UCDavis_e_Null_Stats <- data.table(Phenotype_UCDavis_e_Null, Phenotype_UCDavis_p_Null$Gender)
    colnames(Phenotype_UCDavis_e_Null_Stats) <- c("PhenoVAlue", "Gender")
    Phenotype_UCDavis_e_Null_Stats <- Phenotype_UCDavis_e_Null_Stats[!is.na(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue), ]
    F_Count[i] <- sum(Phenotype_UCDavis_e_Null_Stats$Gender %in% "Female")      
    M_Count[i] <- sum(Phenotype_UCDavis_e_Null_Stats$Gender %in% "Male")
  
  
    p_val3[i] <- summary(gls(Phenotype_UCDavis_e_Null_Stats ~ Phenotype_UCDavis_p_Null_Stats$Gender, na.action = 'na.exclude'))$coefficients[2, 4]
    Coef3[i] <- summary(gls(Phenotype_UCDavis_e_Null_Stats ~ Phenotype_UCDavis_p_Null_Stats$Gender, na.action = 'na.exclude'))$coefficients[2, 1]
  
    fc3[i] <- mean(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue[Phenotype_UCDavis_e_Null_Stats$Gender %in% "Male"])/ mean(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue[Phenotype_UCDavis_e_Null_Stats$Gender %in% "Female"])
    
      
      }else{
          continuous_index[i] <- FALSE
          F_Count[i] <- NA
          M_Count[i] <- NA
  
          p_val3[i] <- NA
          Coef3[i] <- NA
          fc3[i] <- NA
          
          }
    
    
    
    }



p_val_adj3 <- p.adjust(p_val3, method = "fdr")


fwrite(data.table(label = rownames(Phenotype_UCDavis_e_null),F_Count = F_Count, M_Count = M_Count, p_value = p_val3, coefficient = Coef3, fold_change = fc3, adjusted_p_value = p.adjust(p_val3,'fdr')),"XXXXXXXXXXXXXXXXXXXXXXXXXX.csv")


####################################################################################################
## 2b: wildtype mice phenotype -  sex adjusted by body weight
#load bodyweight data from phenotype results 

Phenotype_UCDavisSelect <- wcmc::read_data("XXXXXXXXXXXXXXXXXXXXXXXXXX.xlsx")
Phenotype_UCDavisSelect_p <- Phenotype_UCDavisSelect$p
Phenotype_UCDavisSelect_f <- Phenotype_UCDavisSelect$f
Phenotype_UCDavisSelect_e <- Phenotype_UCDavisSelect$e_matrix

body_weight_all <- Phenotype_UCDavisSelect_e[,Phenotype_UCDavisSelect_p[["phenotype"]] %in% "Body weight"]
names(body_weight_all) <- Phenotype_UCDavisSelect_f$label

body_weight <- body_weight_all[Raw_Phenotype_UCDavis_p_null$label]

#statistics
p_val4 = Coef4 = fc4 = p_val_adj4 = F_Count = M_Count = c()

for (i in 1:nrow(Phenotype_UCDavis_e_null)){
  
  continuous <- sum(is.na(as.numeric(Phenotype_UCDavis_e_null[i, !is.na(Phenotype_UCDavis_e_null[i, ])]))) == 0
  
  
  if(continuous){
    continuous_index[i] <- TRUE
    
    
    body_weight_stats <- body_weight[names(body_weight) %in% Phenotype_UCDavis_p_Null$label]
    
    Phenotype_UCDavis_e_Null_Stats <- data.table(Phenotype_UCDavis_e_Null, Phenotype_UCDavis_p_Null$Gender, body_weight[names(body_weight) %in% Phenotype_UCDavis_p_Null$label])
    colnames(Phenotype_UCDavis_e_Null_Stats) <- c("PhenoVAlue", "Gender", "Weight")
    Phenotype_UCDavis_e_Null_Stats <- Phenotype_UCDavis_e_Null_Stats[!(is.na(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue) | is.na(Phenotype_UCDavis_e_Null_Stats$Weight)), ]
    
    F_Count[i] <- sum(Phenotype_UCDavis_e_Null_Stats$Gender %in% "Female")      
    M_Count[i] <- sum(Phenotype_UCDavis_e_Null_Stats$Gender %in% "Male")
    
    
    p_val4[i] <- summary(gls(Phenotype_UCDavis_e_Null_Stats ~ Phenotype_UCDavis_p_Null_Stats$Gender + body_weight_stats, na.action = 'na.exclude'))$coefficients[2, 4]
    Coef4[i] <- summary(gls(Phenotype_UCDavis_e_Null_Stats ~ Phenotype_UCDavis_p_Null_Stats$Gender + body_weight_stats, na.action = 'na.exclude'))$coefficients[2, 1]
    
    fc4[i] <- mean(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue[Phenotype_UCDavis_e_Null_Stats$Gender %in% "Male"])/ mean(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue[Phenotype_UCDavis_e_Null_Stats$Gender %in% "Female"])
    
    
    
  }else{
    continuous_index[i] <- FALSE
    
    
    F_Count[i] <- NA     
    M_Count[i] <- NA
    p_val4[i] <- NA
    Coef4[i] <- NA
    fc4[i] <- NA
  }
  
}


p_val_adj4 <- p.adjust(p_val4, method = "fdr")
# sum(p_val4[continuous_index] < 0.05,na.rm = TRUE)/sum(continuous_index)


fwrite(data.table(label = rownames(Phenotype_UCDavis_e_null), F_Count = F_Count, M_Count = M_Count, p_value = p_val4, Coefficient = Coef4, fold_change = fc4, adjusted_p_value = p.adjust(p_val4, 'fdr')), "XXXXXXXXXXXXXXXXXXXXXXXXXX.csv")