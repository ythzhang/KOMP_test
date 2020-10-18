# TODO: Add comment

# R version 3.6.0 (2019-04-26) -- "Planting of a Tree"
# Copyright (C) 2019 The R Foundation for Statistical Computing
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# 
# R is free software and comes with ABSOLUTELY NO WARRANTY.
# You are welcome to redistribute it under certain conditions.
# Type 'license()' or 'licence()' for distribution details.
# 
# R is a collaborative project with many contributors.
# Type 'contributors()' for more information and
# 'citation()' on how to cite R or R packages in publications.
# 
# Type 'demo()' for some demos, 'help()' for on-line help, or
# 'help.start()' for an HTML browser interface to help.
# Type 'q()' to quit R.
# 
# Author: YZ
###############################################################################

#Load packages
library(wcmc)
library(data.table)
library(RNOmni)

# Read data file Supplementary Data 1 for the analysis of Metabolomics data in wildtype mice
setwd("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
Null_Metabolomics <- wcmc::read_data("Supplementary Data 1.xlsx")

Null_Metabolomics_p <- Null_Metabolomics$p
Null_Metabolomics_f <- Null_Metabolomics$f
Null_Metabolomics_e <- Null_Metabolomics$e_matrix

# replace any non numeric text with NA
Null_Metabolomics_e[Null_Metabolomics_e == 0] <- NA
Null_Metabolomics_e[Null_Metabolomics_e == "NA"] <- NA
Null_Metabolomics_e[Null_Metabolomics_e == "NaN"] <- NA
Null_Metabolomics_e[Null_Metabolomics_e == "Inf"] <- NA
Null_Metabolomics_e[!is.finite(Null_Metabolomics_e)] <- NA
sum(Null_Metabolomics_e == "Inf")

## 1 wild - sex
## 1a: wildtype mice - sex
### deal with missing values: filter missing value > 70%, replace other missing values.

### subset data for wildtype male mice and calculate percentage of missing values
num_missing_male <- apply(Null_Metabolomics_e, 1, function(x){
  sum(is.na(x[Null_Metabolomics_p$Gender %in% "Male"]))
})

### subset data for wildtype female mice and calculate percentage of missing values
num_missing_female <- apply(Null_Metabolomics_e, 1, function(x){
  sum(is.na(x[Null_Metabolomics_p$Gender %in% "Female"]))
})

### remove metabolite variables that have greater than 70% missing values. 
names(num_missing_male) = names(num_missing_female) = Null_Metabolomics_f$label
missing_index <- ((num_missing_male/ sum(Null_Metabolomics_p$Gender %in% "Male")) >= 0.7) & ((num_missing_female/ sum(Null_Metabolomics_p$Gender %in% "Female")) >= 0.7)
Null_Metabolomics_e <- Null_Metabolomics_e[!missing_index, ]
Null_Metabolomics_f <- Null_Metabolomics_f[!missing_index, ]


### replace missing values by half minimum of the values in each metabolite variable
for(i in 1:nrow(Null_Metabolomics_e)){
  Null_Metabolomics_e[i,is.na(Null_Metabolomics_e[i,])] = 0.5 * min(Null_Metabolomics_e[i,], na.rm = TRUE)
}

rownames(Null_Metabolomics_e) = Null_Metabolomics_f$label

## statistics
########################################################################
#Load packages used for statistics
library(nlme)
test1_Norma = Normality_F = Normality_M = F_Count = M_Count  = Coef1 = p_val1 = fc1 = p_val_adj1 = c()

for (i in 1:nrow(Null_Metabolomics_e)){
  #construct dataframe
  Null_Metabol <- data.table(Intensity = Null_Metabolomics_e[i,], Gender=Null_Metabolomics_p$Gender)
  #remove missing vaues
  Null_Metabol <- Null_Metabol[!is.na(Null_Metabol$Intensity), ]
 
  ### transform data
  test1_Norma[i] <- shapiro.test(Null_Metabol$Intensity)$p.value
  Null_Metabol$Intensity_trans <- rankNorm(Null_Metabol$Intensity)
  Null_Metabolomics_f$CompoundName[i]  
  Null_Metabolomics_f$label[i] 
  
  ### find outlier
  outliers_nullF<-boxplot.stats(Null_Metabol$Intensity_trans[Null_Metabol$Gender%in%'Female'])$out
  outliers_nullM<-boxplot.stats(Null_Metabol$Intensity_trans[Null_Metabol$Gender%in%'Male'])$out
  ### remove outliers if necessary
  if(length(outliers_nullF) > 0 & length(outliers_nullM) > 0){
  
  Label_outlier <- c(which(Null_Metabol$Intensity_trans %in% outliers_nullF[[1]]), which(Null_Metabol$Intensity_trans %in% outliers_nullM[[1]]))
  Null_Metabol_new <- Null_Metabol [-Label_outlier,]
  }else if(length(outliers_nullF) > 0 & length(outliers_nullM) == 0){
    Label_outlier <- c(which(Null_Metabol$Intensity_trans %in% outliers_nullF[[1]]))
    Null_Metabol_new <- Null_Metabol [-Label_outlier, ]
  }else if(length(outliers_nullF) == 0 & length(outliers_nullM) > 0){
    Label_outlier <- c(which(Null_Metabol$Intensity_trans %in% outliers_nullM[[1]]))
    Null_Metabol_new <- Null_Metabol[-Label_outlier, ]
  }else{
    Null_Metabol_new <- Null_Metabol
  }
  # count number for female and male mice
  F_Count[i] <- sum((Null_Metabol_new$Gender %in% "Female"))      
  M_Count[i] <- sum(Null_Metabol_new$Gender %in% "Male")
  
  ### subset data by sex and test normality for WT male mice and female mice
  Null_Metabol_new_M <-  Null_Metabol_new[Null_Metabol_new$Gender %in% "Male"]
  Normality_M[i]<-shapiro.test(Null_Metabol_new_M$Intensity_trans)$p.value
  Null_Metabol_new_F <-  Null_Metabol_new[Null_Metabol_new$Gender %in% "Female"]
  if( nrow(unique(Null_Metabol_new_F)) == 1 ){
    Normality_F[i] <- NA
  }else{
    Normality_F[i] <- shapiro.test(Null_Metabol_new_F$Intensity_trans)$p.value 
  }
  
  Null_e_M <- Null_Metabol_new[Null_Metabol_new$Gender %in% "Male", ]
  Null_e_F <- Null_Metabol_new[Null_Metabol_new$Gender %in% "Female", ]
  
  # test sex as biological variable using gls() model
  WT_sex_stats <- tryCatch({
    
    gls(Intensity_trans ~ Gender, Null_Metabol_new, na.action = 'na.exclude')
  }, error = function(er){
    NA
  })
  
  if(length(WT_sex_stats) <= 1){
    p_val1[i] = NA
    Coef1[i] = NA
  }else{
    
    p_val1[i] = nlme:::summary.gls( WT_sex_stats)$tTable[2,4]
    Coef1[i] = nlme:::summary.gls( WT_sex_stats)$tTable[2,1]
  }
  
  ### calculate fold-changes for each metabolite variable using female mice as reference
  fc1[i]=  mean(Null_Metabol_new$Intensity[Null_Metabol_new$Gender %in% "Male"])/mean(Null_Metabol_new$Intensity[Null_Metabol_new$Gender %in% "Female"])
  
}

### adjustment method using ("BH" or its alias "fdr")
p_val_adj1 = p.adjust(p_val1, method = "fdr")

### calculate significant ratio
sum(p_val1<0.05)/length(p_val1)
sum(p_val_adj1<0.05)/length(p_val1)

#Output result
setwd("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
fwrite(data.table(label = Null_Metabolomics_f$label, Metbaolites = Null_Metabolomics_f$CompoundName, Platform= Null_Metabolomics_f$Assay, InChiKey = Null_Metabolomics_f$InChiKey, PubChemID = Null_Metabolomics_f$PubChemID, SMILES=Null_Metabolomics_f$SMILES, F_Count = F_Count, M_Count = M_Count, Normality_before = test1_Norma, Normality_M = Normality_M, Normality_F = Normality_F, p_value = p_val1, Coefficient= Coef1, fold_change = fc1, adjusted_p_value = p_val_adj1),"1 Metabo_WT_sex effect.csv")




#######################################################################################################
## 1b: wildtype mice - sex as covariate by adjusting to body weight
### load bodyweight data from phenotype results, body weight recorded as "IMPC_HWT_007_001"

Phenotype_weight <- wcmc::read_data("xxxxxxxxxxxxxxxxxxxxxxxxx.xlsx")
Phenotype_weight_p = Phenotype_weight$p
Phenotype_weight_f = Phenotype_weight$f
Phenotype_weight_e = Phenotype_weight$e_matrix

# subset body weight from IMPC phenotype data
body_weight_all = Phenotype_weight_e[Phenotype_weight_f$label %in% "IMPC_HWT_007_001", Phenotype_weight_p$label %in% Null_Metabolomics_p$label]

test1_Norma = Normality_F = Normality_M = F_Count = M_Count = p_val2 = Coef2 = fc2 = p_val_adj2 = c()

for (i in 1 : nrow(Null_Metabolomics_e)){
  
  ###test normality
  test1_Norma[i] <- shapiro.test(Null_Metabolomics_e[i, ])$p.value
  ###transform data
  Null_Metabolomics_e_trans <- rankNorm(Null_Metabolomics_e[i, ])
  
  ### find outlier
  outliers_nullF <- boxplot.stats(Null_Metabolomics_e_trans[Null_Metabolomics_p$Gender %in% 'Female'])$out
  outliers_nullM <- boxplot.stats(Null_Metabolomics_e_trans[Null_Metabolomics_p$Gender %in% 'Male'])$out
  ### remove outliers if necessary
  if( length(outliers_nullF) > 0 & length(outliers_nullM) > 0){
    
    Label_outlier <- c(which(Null_Metabolomics_e_trans %in% outliers_nullF[[1]]), which(Null_Metabolomics_e_trans %in% outliers_nullM[[1]]))
    Null_Metabolomics_e_trans_New <- Null_Metabolomics_e_trans[-Label_outlier]
    Null_raw_e_new <- Null_Metabolomics_e[i,][-Label_outlier]
    Null_Metabolomics_p_New <- Null_Metabolomics_p[-Label_outlier]
  }else if(length(outliers_nullF) > 0 & length(outliers_nullM) == 0){
    Label_outlier <- c(which(Null_Metabolomics_e_trans %in% outliers_nullF[[1]]))
    Null_Metabolomics_e_trans_New <- Null_Metabolomics_e_trans[-Label_outlier]
    Null_raw_e_new <- Null_Metabolomics_e[i,][-Label_outlier]
    Null_Metabolomics_p_New <- Null_Metabolomics_p[-Label_outlier]
  }else if(length(outliers_nullF) == 0 & length(outliers_nullM) > 0){
    Label_outlier <- c(which(Null_Metabolomics_e_trans %in% outliers_nullM[[1]]))
    Null_Metabolomics_e_trans_New <- Null_Metabolomics_e_trans[-Label_outlier]
    Null_raw_e_new <- Null_Metabolomics_e[i,][-Label_outlier]
    Null_Metabolomics_p_New <- Null_Metabolomics_p[-Label_outlier]
  }else{
    Null_Metabolomics_e_trans_New <- Null_Metabolomics_e_trans
    Null_raw_e_new <- Null_Metabolomics_e[i,]
    Null_Metabolomics_p_New <- Null_Metabolomics_p
  }
  
  Null_Metabolomics_e_Stats <- data.table(Null_Metabolomics_p_New$Gender, Null_raw_e_new, Null_Metabolomics_e_trans_New, body_weight_all[names(body_weight_all) %in% Null_Metabolomics_p_New$label])
  colnames(Null_Metabolomics_e_Stats) <- c("Gender","MetaboVAlue_raw","MetaboVAlue","Weight")
  ### need to remove data for which body weight was not available
  Null_Metabolomics_e_Stats <- Null_Metabolomics_e_Stats[!is.na(Null_Metabolomics_e_Stats$Weight),]

  ### Count number of female and male mice if metabolite values and body weight are both available
  F_Count[i] <- sum(Null_Metabolomics_e_Stats$Gender %in% "Female")      
  M_Count[i] <- sum(Null_Metabolomics_e_Stats$Gender %in% "Male")
  
  ### test normality for WT male mice and female mice
  if( length(unique(Null_Metabolomics_e_Stats$MetaboVAlue[Null_Metabolomics_e_Stats$Gender%in%"Female"])) == 1 ){
    Normality_M[i] <- shapiro.test(Null_Metabolomics_e_Stats$MetaboVAlue[Null_Metabolomics_e_Stats$Gender %in% "Male"])$p.value
    Normality_F[i] <- NA
  }else{
    Normality_M[i] <- shapiro.test(Null_Metabolomics_e_Stats$MetaboVAlue[Null_Metabolomics_e_Stats$Gender %in% "Male"])$p.value
    Normality_F[i] <- shapiro.test(Null_Metabolomics_e_Stats$MetaboVAlue[Null_Metabolomics_e_Stats$Gender %in% "Female"])$p.value
  }
  
  ### test sex as covariable and adjusted to bosy weight using gls() model
  WT_sex_stats <- tryCatch({
    
    gls(MetaboVAlue ~ Gender + Weight, Null_Metabolomics_e_Stats, na.action = 'na.exclude')
  }, error = function(er){
    NA
  })
  
  if(length(WT_sex_stats) <= 1){
    p_val1[i] = NA
    Coef1[i] = NA
    
  }else{
    
    p_val2[i] = nlme:::summary.gls(WT_sex_stats)$tTable[2,4]
    Coef2[i] = nlme:::summary.gls(WT_sex_stats)$tTable[2,1]
  }
  
  ### calculate fold-changes for each metabolite variable using female mice as reference
  fc2[i]=  mean(Null_Metabolomics_e_Stats$MetaboVAlue_raw[Null_Metabolomics_e_Stats$Gender %in% "Male"])/mean(Null_Metabolomics_e_Stats$MetaboVAlue_raw[Null_Metabolomics_e_Stats$Gender %in% "Female"])
}

### adjustment method using ("BH" or its alias "fdr")
p_val_adj2 = p.adjust(p_val2, method ="fdr")

### calculate significant ratio
sum(p_val2<0.05)/length(p_val2)
sum(p_val_adj2<0.05)/length(p_val2)

#Output result
setwd("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
fwrite(data.table(label = Null_Metabolomics_f$label, Metbaolites = Null_Metabolomics_f$CompoundName, Platform= Null_Metabolomics_f$Assay, InChiKey = Null_Metabolomics_f$InChiKey, PubChemID = Null_Metabolomics_f$PubChemID, SMILES=Null_Metabolomics_f$SMILES, F_Count = F_Count, M_Count = M_Count, Normality_before = test1_Norma, Normality_M = Normality_M, Normality_F = Normality_F, p_value = p_val2, Coefficient = Coef2, fold_change = fc2, adjusted_p_value = p_val_adj2),"1b Metabol_adjust byBW.csv")

