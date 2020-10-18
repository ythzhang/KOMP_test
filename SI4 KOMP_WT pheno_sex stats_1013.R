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
# Load packages
library(wcmc)
library(data.table)
library(RNOmni)
# Read data file Supplementary Data 2 for the analysis of IMPC phenotype data in wildtype mice
setwd("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
Phenotype_UCDavis <- wcmc::read_data("Supplementary Data 2.xlsx")
Raw_Phenotype_UCDavis_p_null <- Phenotype_UCDavis$p
Raw_Phenotype_UCDavisS_f_null <- Phenotype_UCDavis$f
Raw_Phenotype_UCDavis_e_null <- Phenotype_UCDavis$e_matrix
rownames(Raw_Phenotype_UCDavis_e_null) <- Raw_Phenotype_UCDavisS_f_null$label
### check data structure and set as numeric
Raw_Phenotype_UCDavis_e_null = apply(Raw_Phenotype_UCDavis_e_null,2,as.numeric)
## 1 wild - sex
## 1a: wildtype mice - sex
### deal with missing values: filter missing value > 70%, replace other missing values.
### subset data for wildtype male mice and calculate percentage of missing values
num_missing_male = apply(Raw_Phenotype_UCDavis_e_null,1,function(x){
  sum(is.na(x[Raw_Phenotype_UCDavis_p_null$Gender %in% "Male"]))
})
### subset data for wildtype female mice and calculate percentage of missing values
num_missing_female = apply(Raw_Phenotype_UCDavis_e_null,1,function(x){
  sum(is.na(x[Raw_Phenotype_UCDavis_p_null$Gender %in% "Female"]))
})
### remove metabolite variables that have greater than 70% missing values.
missing_index = ((num_missing_male/sum(Raw_Phenotype_UCDavis_p_null$Gender %in% "Male")) >= 0.7) & ((num_missing_female/sum(Raw_Phenotype_UCDavis_p_null$Gender %in% "Female")) >= 0.7)
Phenotype_UCDavis_e_null = Raw_Phenotype_UCDavis_e_null[!missing_index,]
rownames(Phenotype_UCDavis_e_null) = rownames(Raw_Phenotype_UCDavis_e_null)[!missing_index]
nrow(Phenotype_UCDavis_e_null)
## statistics
########################################################################
#Load packages used for statistics
library(nlme)
F_Count = M_Count = Norma_F = Norma_M = p_val3 = fc3 = p_val_adj3 = Coef3 = continuous_index = c()
for (i in 1:nrow(Phenotype_UCDavis_e_null)){
  #confirm continuous trait
  continuous = sum(is.na(as.numeric(Phenotype_UCDavis_e_null[i,!is.na(Phenotype_UCDavis_e_null[i,])]))) == 0
  if(continuous){
    continuous_index[i] = TRUE
    ### find outlier
    outliers_nullF <- boxplot.stats(Phenotype_UCDavis_e_null[i,][Raw_Phenotype_UCDavis_p_null$Gender %in% 'Female'])$out
    outliers_nullM <- boxplot.stats(Phenotype_UCDavis_e_null[i,][Raw_Phenotype_UCDavis_p_null$Gender %in% 'Male'])$out
    ### remove outliers if necessary
    if( length(outliers_nullF) > 0 & length(outliers_nullM) > 0){
      Label_outlier <- c(which(Phenotype_UCDavis_e_null[i,] %in% outliers_nullF), which(Phenotype_UCDavis_e_null[i,] %in% outliers_nullM[[1]]))
      Phenotype_UCDavis_e_Null_New <- Phenotype_UCDavis_e_null[i,][-Label_outlier]
      Phenotype_UCDavis_p_Null_New <- Raw_Phenotype_UCDavis_p_null[-Label_outlier]
    }else if(length(outliers_nullF) > 0 & length(outliers_nullM) == 0){
      Label_outlier <- c(which(Phenotype_UCDavis_e_null[i,] %in% outliers_nullF))
      Phenotype_UCDavis_e_Null_New <- Phenotype_UCDavis_e_null[i,][-Label_outlier]
      Phenotype_UCDavis_p_Null_New <- Raw_Phenotype_UCDavis_p_null[-Label_outlier]
    }else if(length(outliers_nullF) == 0 & length(outliers_nullM) > 0){
      Label_outlier <- c(which(Phenotype_UCDavis_e_null[i,] %in% outliers_nullM))
      Phenotype_UCDavis_e_Null_New <- Phenotype_UCDavis_e_null[i,][-Label_outlier]
      Phenotype_UCDavis_p_Null_New <- Raw_Phenotype_UCDavis_p_null[-Label_outlier]
    }else{
      Phenotype_UCDavis_e_Null_New <- Phenotype_UCDavis_e_null[i,]
      Phenotype_UCDavis_p_Null_New <- Raw_Phenotype_UCDavis_p_null
    }
    # construct dataframe
    Phenotype_UCDavis_e_Null_Stats <- data.table(Phenotype_UCDavis_e_Null_New,Phenotype_UCDavis_p_Null_New$Gender)
    colnames(Phenotype_UCDavis_e_Null_Stats)<-c("PhenoVAlue","Gender")
    Phenotype_UCDavis_e_Null_Stats <- Phenotype_UCDavis_e_Null_Stats[!is.na(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue),]
    ### calculate number of available values for female and male mice
    F_Count[i] <- sum(Phenotype_UCDavis_e_Null_Stats$Gender%in%"Female")      
    M_Count[i] <- sum(Phenotype_UCDavis_e_Null_Stats$Gender%in%"Male")
    ### transform data and test normality for male and female mice
    Phenotype_UCDavis_e_Null_Stats$PhenoVAlue_trans <- rankNorm(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue)
    if (length(unique(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue_trans)) == 1){
      Norma_F[i] <- NA
      Norma_M[i] <- NA
    }else{
      Norma_F[i] <- shapiro.test(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue_trans[Phenotype_UCDavis_e_Null_Stats$Gender %in% "Female"])$p.value
      Norma_M[i] <- shapiro.test(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue_trans[Phenotype_UCDavis_e_Null_Stats$Gender %in% "Male"])$p.value
    }
    # test sex as biological variable using gls() model
    WT_sex_stats <- tryCatch({
      gls(PhenoVAlue_trans ~ Gender, Phenotype_UCDavis_e_Null_Stats, na.action = 'na.exclude')
    }, error = function(er){
      NA
      })
    #p-values
      if(length(WT_sex_stats) <= 1){
        p_val3[i] = NA
        Coef3[i] = NA
      }else{
        p_val3[i]= nlme:::summary.gls( WT_sex_stats)$tTable[2,4]
        Coef3[i]= nlme:::summary.gls( WT_sex_stats)$tTable[2,1]
      }
    ### calculate fold-changes for each metabolite variable using female mice as reference
      fc3[i]=  mean(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue[Phenotype_UCDavis_e_Null_Stats$Gender %in% "Male"])/mean(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue[Phenotype_UCDavis_e_Null_Stats$Gender %in% "Female"])
    }else{
      continuous_index[i] = FALSE
      F_Count[i] = NA
      M_Count[i] = NA
      p_val3[i] = NA
      Coef3[i] = NA
      fc3[i] = NA
    }
}
### adjustment method using ("BH" or its alias "fdr")
p_val_adj3 = p.adjust(p_val3,method = "fdr")
### calculate significant ratio
sum(p_val3[continuous_index]<0.05,na.rm = TRUE)/sum(continuous_index)
# sum(fdr3[continuous_index]<0.05,na.rm = TRUE)/sum(continuous_index)
#Output result
setwd("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
fwrite(data.table(label = rownames(Phenotype_UCDavis_e_null), F_Count=F_Count, M_Count=M_Count, Normality_F = Norma_F, Normality_M = Norma_M, p_value = p_val3, coefficient=Coef3, fold_change = fc3, adjusted_p_value = p.adjust(p_val3,'fdr')),"2 Pheno WT_sex effect.csv")


####################################################################################################
## 2b phenotype -  sex effect by adjusting to body weight
### load bodyweight data from phenotype results, body weight measured at the same time point as each parameter was measured
Phenotype_UCDavis <- wcmc::read_data("Supplementary Data 2.xlsx")
Phenotype_UCDavis_p_null = Phenotype_UCDavis$p
Phenotype_UCDavisS_f_null = Phenotype_UCDavis$f
Phenotype_UCDavis_e_null = Phenotype_UCDavis$e_matrix
# subset body weight from IMPC phenotype data
Phenotype_weight = wcmc::read_data("xxxxxxxxxxxxxxxxxxxxxxxxx.xlsx", sheet = "WT_weight")
Phenotype_weight_p = Phenotype_weight$p
Phenotype_weight_f = Phenotype_weight$f
Phenotype_weight_e = Phenotype_weight$e_matrix
continuous_index = Norma_F = Norma_M = p_val4 = Coef4 = fc4 = p_val_adj4 = F_Count = M_Count = c()
for (i in 1:nrow(Phenotype_UCDavis_e_null)){
  #confirm if parameter is continuous tratits or not
  continuous = sum(is.na(as.numeric(Phenotype_UCDavis_e_null[i,!is.na(Phenotype_UCDavis_e_null[i,])]))) == 0
  if(continuous){
    continuous_index[i] = TRUE
    ### find outlier
    outliers_nullF <- boxplot.stats(Phenotype_UCDavis_e_null[i,][Phenotype_UCDavis_p_null$Gender %in% 'Female'])$out
    outliers_nullM <- boxplot.stats(Phenotype_UCDavis_e_null[i,][Phenotype_UCDavis_p_null$Gender %in% 'Male'])$out
    ### remove outliers if necessary
    if( length(outliers_nullF) > 0 & length(outliers_nullM)>0){
      Label_outlier <- c(which(Phenotype_UCDavis_e_null[i,] %in% outliers_nullF),which(Phenotype_UCDavis_e_null[i,] %in% outliers_nullM[[1]]))
      Phenotype_UCDavis_e_Null_New <- Phenotype_UCDavis_e_null[i,][-Label_outlier]
      Phenotype_UCDavis_p_Null_New <- Phenotype_UCDavis_p_null[-Label_outlier]
    }else if(length(outliers_nullF) > 0 & length(outliers_nullM) == 0){
      Label_outlier <- c(which(Phenotype_UCDavis_e_null[i,] %in% outliers_nullF))
      Phenotype_UCDavis_e_Null_New <- Phenotype_UCDavis_e_null[i,][-Label_outlier]
      Phenotype_UCDavis_p_Null_New <- Phenotype_UCDavis_p_null[-Label_outlier]
    }else if(length(outliers_nullF) == 0 & length(outliers_nullM) > 0){
      Label_outlier <- c(which(Phenotype_UCDavis_e_null[i,] %in% outliers_nullM))
      Phenotype_UCDavis_e_Null_New <- Phenotype_UCDavis_e_null[i,][-Label_outlier]
      Phenotype_UCDavis_p_Null_New <- Phenotype_UCDavis_p_null[-Label_outlier]
    }else{
      Phenotype_UCDavis_e_Null_New <- Phenotype_UCDavis_e_null[i,]
      Phenotype_UCDavis_p_Null_New <- Phenotype_UCDavis_p_null
    }
    #construct dataframe
    Phenotype_UCDavis_e_Null_Stats <- data.table(Phenotype_UCDavis_e_Null_New,Phenotype_UCDavis_p_Null_New$Gender,Phenotype_weight_e[i,Phenotype_weight_p$label %in% Phenotype_UCDavis_p_Null_New$label])
    colnames(Phenotype_UCDavis_e_Null_Stats) <- c("PhenoVAlue","Gender","Weight")
    ### need to remove data for which phenotype value or body weight was not available
    Phenotype_UCDavis_e_Null_Stats <- Phenotype_UCDavis_e_Null_Stats[!(is.na(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue)|is.na(Phenotype_UCDavis_e_Null_Stats$Weight)),]
    ### Count number of female and male mice if metabolite values and body weight are both available
    F_Count[i] <- sum(Phenotype_UCDavis_e_Null_Stats$Gender %in% "Female")      
    M_Count[i] <- sum(Phenotype_UCDavis_e_Null_Stats$Gender %in% "Male")
    ### transform data and check if data points are numeric  
    Phenotype_UCDavis_e_Null_Stats$PhenoVAlue_trans <- as.numeric(rankNorm(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue))
    class(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue_trans)
    ### test normality for WT male mice and female mice
    if (length(unique(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue_trans)) == 1){
      Norma_F[i] <- NA
      Norma_M[i] <- NA
    }else{
      Norma_F[i] <- shapiro.test(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue_trans[Phenotype_UCDavis_e_Null_Stats$Gender %in% "Female"])$p.value
      Norma_M[i] <- shapiro.test(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue_trans[Phenotype_UCDavis_e_Null_Stats$Gender %in% "Male"])$p.value
    }
    ### test sex as covariable and adjusted to bosy weight using gls() model
    WT_sex_stats <- tryCatch({
      gls(PhenoVAlue_trans ~ Gender + Weight, Phenotype_UCDavis_e_Null_Stats, na.action = 'na.exclude')
    }, error = function(er){
      NA
    })
    #p-values
    if(length(WT_sex_stats) <= 1){
      p_val4[i] = NA
      Coef4[i] = NA
    }else{
      p_val4[i]= nlme:::summary.gls( WT_sex_stats)$tTable[2,4]
      Coef4[i]= nlme:::summary.gls( WT_sex_stats)$tTable[2,1]
    }
    ### calculate fold-changes for each metabolite variable using female mice as reference
    fc4[i]=  mean(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue[Phenotype_UCDavis_e_Null_Stats$Gender %in% "Male"])/mean(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue[Phenotype_UCDavis_e_Null_Stats$Gender %in% "Female"])
  }else{
    continuous_index[i] = FALSE
    F_Count[i] <- NA     
    M_Count[i] <- NA
    Norma_F[i] <- NA
    Norma_M[i] <- NA
    p_val4[i] <- NA
    Coef4[i] <- NA
    fc4[i] <- NA
  }
}
### adjustment method using ("BH" or its alias "fdr")
p_val_adj4 = p.adjust(p_val4,method ="fdr")
### calculate significant ratio
sum(p_val4[continuous_index]<0.05,na.rm = TRUE)/sum(continuous_index)
#Output result
setwd("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
fwrite(data.table(label = Phenotype_UCDavisS_f_null$label, F_Count = F_Count, M_Count = M_Count, Normality_F = Norma_F, Normality_M = Norma_M, p_value = p_val4, Coefficient = Coef4, fold_change = fc4, adjusted_p_value = p.adjust(p_val4,'fdr')),"3 Pheno WT_sex effect_adjust weight.csv")
