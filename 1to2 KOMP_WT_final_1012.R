# Null_Metabolomics = read.csv("C:\\Users\\yzhan\\Documents\\Fiehn\\3_KOMP\\7_Sexual D_20190605\\Jennly-Zhang-MetaboliteSD_ver03-master\\Jennly-Zhang-MetaboliteSD_ver03-master\\Null_Metabolomics.csv")
library(wcmc)
setwd('C:\\Users\\yzhan\\Documents\\Fiehn\\3_KOMP\\7_Sexual D_20190605\\Jennly Zhang KOMP SD_ver10\\Input\\ONLINE_data')
Null_Metabolomics = wcmc::read_data("Null Metabolomics_online.xlsx")

# Null_Metabolomics<-as.data.frame(Null_Metabolomics)
Null_Metabolomics_p = Null_Metabolomics$p
Null_Metabolomics_f = Null_Metabolomics$f
Null_Metabolomics_e = Null_Metabolomics$e_matrix
Null_Metabolomics_e[Null_Metabolomics_e==0]<-NA
Null_Metabolomics_e[Null_Metabolomics_e=="NA"]<-NA
Null_Metabolomics_e[Null_Metabolomics_e=="NaN"]<-NA
Null_Metabolomics_e[Null_Metabolomics_e=="Inf"]<-NA
Null_Metabolomics_e[!is.finite(Null_Metabolomics_e)] = NA
sum(Null_Metabolomics_e=="Inf")

## 1 wild - sex
### deal with missing value. 
#### filter missing value > 80%


num_missing_male = apply(Null_Metabolomics_e,1,function(x){
  sum(is.na(x[Null_Metabolomics_p$Gender %in% "Male"]))
})
num_missing_female = apply(Null_Metabolomics_e,1,function(x){
  sum(is.na(x[Null_Metabolomics_p$Gender %in% "Female"]))
})
names(num_missing_male) = names(num_missing_female) = Null_Metabolomics_f$label
missing_index = ((num_missing_male/sum(Null_Metabolomics_p$Gender %in% "Male")) >= 0.7)&((num_missing_female/sum(Null_Metabolomics_p$Gender %in% "Female")) >= 0.7)
Null_Metabolomics_e = Null_Metabolomics_e[!missing_index,]
Null_Metabolomics_f = Null_Metabolomics_f[!missing_index,]


#### filter missing value > 80%
for(i in 1:nrow(Null_Metabolomics_e)){
  Null_Metabolomics_e[i,is.na(Null_Metabolomics_e[i,])] = 0.5 * min(Null_Metabolomics_e[i,], na.rm = TRUE)
}

rownames(Null_Metabolomics_e) = Null_Metabolomics_f$label

## statistics
########################################################################








library(RNOmni)
library(nlme)
test1_Norma=test2_Norma=Normality_F=Normality_M= p_val1=fc1=p_val_adj1=F_Count=M_Count=Coef1=Mean.diff =Std.error=DF=Residual=Effect.size = c()

for (i in 1:nrow(Null_Metabolomics_e)){

  Null_Metabol <- data.table(Intensity = Null_Metabolomics_e[i,], Gender=Null_Metabolomics_p$Gender)
  
  Null_Metabol <- Null_Metabol[!is.na(Null_Metabol$Intensity), ]
  # Null_Metabolomics_e <- Null_Metabolomics_e[i,!is.na(Null_Metabolomics_e[i,])]
  test1_Norma[i] <- shapiro.test(Null_Metabol$Intensity)$p.value
  Null_Metabol$Intensity_trans <- rankNorm(Null_Metabol$Intensity)
  # plotNormalHistogram(Null_Metabolomics_e_Blom)
  test2_Norma[i] <- shapiro.test(Null_Metabol$Intensity_trans)$p.value

  Null_Metabolomics_f$CompoundName[i]  
  Null_Metabolomics_f$label[i] 
  
  # Null_Metabolomics_e_Blom <- Null_Metabolomics_e_Blom[!is.na(Null_Metabolomics_e_Blom)]
  outliers_nullF<-boxplot.stats(Null_Metabol$Intensity_trans[Null_Metabol$Gender%in%'Female'])$out
  outliers_nullM<-boxplot.stats(Null_Metabol$Intensity_trans[Null_Metabol$Gender%in%'Male'])$out
  
  if( length(outliers_nullF)>0 & length(outliers_nullM)>0){
    
  Label_outlier<-c(which(Null_Metabol$Intensity_trans %in% outliers_nullF[[1]]),which(Null_Metabol$Intensity_trans %in% outliers_nullM[[1]]))
  Null_Metabol_new <- Null_Metabol [-Label_outlier,]
  
  }else if(length(outliers_nullF)>0 & length(outliers_nullM)==0){
    
    Label_outlier<-c(which(Null_Metabol$Intensity_trans %in% outliers_nullF[[1]]))
    Null_Metabol_new<-Null_Metabol [-Label_outlier, ]

    
  }else if(length(outliers_nullF)==0 & length(outliers_nullM)>0){
    
    Label_outlier<-c(which(Null_Metabol$Intensity_trans %in% outliers_nullM[[1]]))
    Null_Metabol_new <- Null_Metabol[-Label_outlier, ]
    
  }else{
    
    Null_Metabol_new <- Null_Metabol
  
  }
  
  
  F_Count[i] <- sum((Null_Metabol_new$Gender %in% "Female"))      
  M_Count[i] <- sum(Null_Metabol_new$Gender %in% "Male")
  
  
  Null_Metabol_new_M <-  Null_Metabol_new[Null_Metabol_new$Gender %in% "Male"]
  Normality_M[i]<-shapiro.test(Null_Metabol_new_M$Intensity_trans)$p.value
  Null_Metabol_new_F <-  Null_Metabol_new[Null_Metabol_new$Gender %in% "Female"]
  
  if( nrow(unique(Null_Metabol_new_F))==1 ){
    Normality_F[i] <- NA
    
  }else{
    Normality_F[i]<-shapiro.test(Null_Metabol_new_F$Intensity_trans)$p.value 
  }
  
  
  Null_e_M <- Null_Metabol_new[Null_Metabol_new$Gender %in% "Male",]
  Null_e_F <- Null_Metabol_new[Null_Metabol_new$Gender %in% "Female",]
  
  
  
  
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
  
  
  
  # lm_result<- summary(lm(Null_Metabol_new$Blom ~ Null_Metabol_new$Gender, na.action = 'na.exclude'))
  # p_val1[i]= summary(lm(Null_Metabol_new$Blom ~ Null_Metabol_new$Gender, na.action = 'na.exclude'))$coefficients[2,4]
  # Coef1[i]= summary(lm(Null_Metabol_new$Blom ~ Null_Metabol_new$Gender, na.action = 'na.exclude'))$coefficients[2,1]
  fc1 [i]=  mean(Null_Metabol_new$Intensity[Null_Metabol_new$Gender %in% "Male"])/mean(Null_Metabol_new$Intensity[Null_Metabol_new$Gender %in% "Female"])
  
  # fc1_2 [i]=  mean(Null_raw_e_new_M)/mean(Null_raw_e_new_F)
  # Mean.diff[i] <- mean(Null_Metabol_new$Intensity[Null_Metabol_new$Gender %in% "Male"])-mean(Null_Metabol_new$Intensity[Null_Metabol_new$Gender %in% "Female"])
  # Std.error[i]<-lm_result$coefficients[2,2]
  # Residual[i]<-lm_result$sigma
  # DF[i]<-lm_result$df[2]
  # Effect.size[i]<-abs(Mean.diff[i])/ Residual[i]

}

p_val_adj1 = p.adjust(p_val1, method = "fdr")


sum(p_val1<0.05)/length(p_val1)
sum(p_val_adj1<0.05)/length(p_val1)

library(data.table)

fwrite(data.table(label = Null_Metabolomics_f$label, Metbaolites = Null_Metabolomics_f$CompoundName,Platform= Null_Metabolomics_f$Assay, InChiKey = Null_Metabolomics_f$InChiKey, PubChemID = Null_Metabolomics_f$PubChemID, SMILES=Null_Metabolomics_f$SMILES, F_Count=F_Count, M_Count =M_Count, Normality_before = test1_Norma, Normality_before = test2_Norma, Normality_M = Normality_M, Normality_F = Normality_F, p_value = p_val1, Coefficient=Coef1, fold_change = fc1, adjusted_p_value = p_val_adj1,Mean.diff=Mean.diff,Std.error=Std.error,Residual=Residual,D_freedom=DF,Effect.size=Effect.size),"C:\\Users\\yzhan\\Documents\\Fiehn\\3_KOMP\\7_Sexual D_20190605\\Jennly Zhang KOMP SD_ver10\\Output\\1 Metabo_wild_sex_outlier_Norm_202001012.csv")




# ----------------------------------------------------------------------------------------------------
## 1b adjust by body weight

Phenotype_weight = wcmc::read_data("C:\\Users\\yzhan\\Documents\\Fiehn\\3_KOMP\\_01 New extracted phenotype_1008\\9_final_pehnotype dataset\\2_Pheno_WT_206parameter.xlsx", sheet = "WT_weight")
Phenotype_weight_p = Phenotype_weight$p
Phenotype_weight_f = Phenotype_weight$f
Phenotype_weight_e = Phenotype_weight$e_matrix


body_weight_all = Phenotype_weight_e[Phenotype_weight_f$label %in% "IMPC_HWT_007_001",Phenotype_weight_p$label %in% Null_Metabolomics_p$label]
# names(body_weight_all) = Phenotype_UCDavisSelect_f$label
# body_weight = body_weight_all[Null_Metabolomics_p$label]

# names(body_weight)<-Null_Metabolomics_p$label



## statistics
library(RNOmni)
test1_Norma=test2_Norma=Normality_F=Normality_M= p_val2=fc2=p_val_adj2=Coef2 = F_Count=M_Count=Coef1=Mean.diff =Std.error=DF=Residual=Effect.size = c()


for (i in 1:nrow(Null_Metabolomics_e)){
  
  test1_Norma[i] <- shapiro.test(Null_Metabolomics_e[i,])$p.value
  Null_Metabolomics_e_trans <- rankNorm(Null_Metabolomics_e[i,])
  # plotNormalHistogram(Null_Metabolomics_e_Blom)
  test2_Norma[i] <- shapiro.test(Null_Metabolomics_e_trans)$p.value
  
  outliers_nullF<-boxplot.stats(Null_Metabolomics_e_trans[Null_Metabolomics_p$Gender%in%'Female'])$out
  outliers_nullM<-boxplot.stats(Null_Metabolomics_e_trans[Null_Metabolomics_p$Gender%in%'Male'])$out
  
  if( length(outliers_nullF)>0 & length(outliers_nullM)>0){
    
    Label_outlier<-c(which(Null_Metabolomics_e_trans %in% outliers_nullF[[1]]),which(Null_Metabolomics_e_trans %in% outliers_nullM[[1]]))
    Null_Metabolomics_e_trans_New<-Null_Metabolomics_e_trans[-Label_outlier]
    Null_raw_e_new <-Null_Metabolomics_e[i,][-Label_outlier]
    Null_Metabolomics_p_New<-Null_Metabolomics_p[-Label_outlier]
    
  }else if(length(outliers_nullF)>0 & length(outliers_nullM)==0){
    
    Label_outlier <- c(which(Null_Metabolomics_e_trans %in% outliers_nullF[[1]]))
    Null_Metabolomics_e_trans_New <- Null_Metabolomics_e_trans[-Label_outlier]
    Null_raw_e_new <-Null_Metabolomics_e[i,][-Label_outlier]
    Null_Metabolomics_p_New<-Null_Metabolomics_p[-Label_outlier]
    
  }else if(length(outliers_nullF)==0 & length(outliers_nullM)>0){
    
    Label_outlier <- c(which(Null_Metabolomics_e_trans %in% outliers_nullM[[1]]))
    Null_Metabolomics_e_trans_New<-Null_Metabolomics_e_trans[-Label_outlier]
    Null_raw_e_new <-Null_Metabolomics_e[i,][-Label_outlier]
    Null_Metabolomics_p_New<-Null_Metabolomics_p[-Label_outlier]
    
  }else{
    Null_Metabolomics_e_trans_New <- Null_Metabolomics_e_trans
    Null_raw_e_new <-Null_Metabolomics_e[i,]
    Null_Metabolomics_p_New<-Null_Metabolomics_p
  }
  
  # 
  # F_Count[i] <- sum(Null_Metabolomics_p_New$Gender%in%"Female")      
  # M_Count[i] <- sum(Null_Metabolomics_p_New$Gender%in%"Male")
  # 
  
  Null_Metabolomics_e_Stats <- data.table(Null_Metabolomics_p_New$Gender, Null_raw_e_new, Null_Metabolomics_e_trans_New, body_weight_all[names(body_weight_all) %in% Null_Metabolomics_p_New$label])
  
  
  colnames(Null_Metabolomics_e_Stats) <- c("Gender","MetaboVAlue_raw","MetaboVAlue","Weight")
  
  
  Null_Metabolomics_e_Stats<-Null_Metabolomics_e_Stats[!is.na(Null_Metabolomics_e_Stats$Weight),]
  # Null_Metabolomics_e_Stats$weight_blom <-rankNorm(Null_Metabolomics_e_Stats$Weight)
  
  # body_weight_stats<-body_weight[names(body_weight) %in% Null_Metabolomics_p_New$label]
  F_Count[i] <- sum(Null_Metabolomics_e_Stats$Gender%in%"Female")      
  M_Count[i] <- sum(Null_Metabolomics_e_Stats$Gender%in%"Male")
  
  
  # Normality_M[i]<- shapiro.test(Null_Metabolomics_e_Stats$MetaboVAlue[Null_Metabolomics_e_Stats$Gender %in% "Male"])$p.value
  
  if( length(unique(Null_Metabolomics_e_Stats$MetaboVAlue[Null_Metabolomics_e_Stats$Gender%in%"Female"])) == 1 ){
    Normality_M[i] <- shapiro.test(Null_Metabolomics_e_Stats$MetaboVAlue[Null_Metabolomics_e_Stats$Gender %in% "Male"])$p.value
    Normality_F[i] <- NA
    
  }else{
    Normality_M[i] <- shapiro.test(Null_Metabolomics_e_Stats$MetaboVAlue[Null_Metabolomics_e_Stats$Gender %in% "Male"])$p.value
    Normality_F[i] <- shapiro.test(Null_Metabolomics_e_Stats$MetaboVAlue[Null_Metabolomics_e_Stats$Gender %in% "Female"])$p.value
  }
  
  
  
  
  
  WT_sex_stats <- tryCatch({
    
    # gls(MetaboVAlue_raw ~ Gender + Weight, Null_Metabolomics_e_Stats, na.action = 'na.exclude')
    gls(MetaboVAlue ~ Gender + Weight, Null_Metabolomics_e_Stats, na.action = 'na.exclude')
    
  }, error = function(er){
    NA
  })
  
  if(length(WT_sex_stats) <= 1){
    p_val1[i] = NA
    Coef1[i] = NA
    
  }else{
    
    p_val2[i] = nlme:::summary.gls( WT_sex_stats)$tTable[2,4]
    Coef2[i] = nlme:::summary.gls( WT_sex_stats)$tTable[2,1]
    
  }
  

  
  # p_val2[i]= summary(lm(Null_Metabolomics_e_Stats$MetaboVAlue ~ Null_Metabolomics_e_Stats$Gender + Null_Metabolomics_e_Stats$weight_blom,na.action = 'na.exclude'))$coefficients[2,4]
  # Coef2[i]= summary(lm(Null_Metabolomics_e_Stats$MetaboVAlue ~ Null_Metabolomics_e_Stats$Gender + Null_Metabolomics_e_Stats$weight_blom,na.action = 'na.exclude'))$coefficients[2,1]
  
  fc2 [i]=  mean(Null_Metabolomics_e_Stats$MetaboVAlue_raw[Null_Metabolomics_e_Stats$Gender %in% "Male"])/mean(Null_Metabolomics_e_Stats$MetaboVAlue_raw[Null_Metabolomics_e_Stats$Gender %in% "Female"])

  
  
  
    
}


p_val_adj2 = p.adjust(p_val2,method ="fdr")

sum(p_val2<0.05)/length(p_val2)
sum(p_val_adj2<0.05)/length(p_val2)

setwd('C:\\Users\\yzhan\\Documents\\Fiehn\\3_KOMP\\7_Sexual D_20190605\\Jennly Zhang KOMP SD_ver10\\Output')
fwrite(data.table(label = Null_Metabolomics_f$label, Metbaolites = Null_Metabolomics_f$CompoundName,Platform= Null_Metabolomics_f$Assay, InChiKey = Null_Metabolomics_f$InChiKey, PubChemID = Null_Metabolomics_f$PubChemID, SMILES=Null_Metabolomics_f$SMILES, F_Count=F_Count, M_Count =M_Count, Normality_before = test1_Norma, Normality_before = test2_Norma, Normality_M = Normality_M, Normality_F = Normality_F,p_value = p_val2, Coefficient=Coef2,fold_change = fc2, adjusted_p_value = p_val_adj2),"1b Metabol_raw intensity_adjust byBW(norm)_1014.csv")


## 2 adjust by body weight
# ------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------



# ===================================================================================================================================================
## 2 phenotype -  sex
library(wcmc)
setwd('C:\\Users\\yzhan\\Documents\\Fiehn\\3_KOMP\\7_Sexual D_20190605\\Jennly Zhang KOMP SD_ver09\\Input')
Raw_Phenotype_UCDavis = wcmc::read_data("Raw_Phenotype_UCDavis.xlsx")
Raw_Phenotype_UCDavis_p = Raw_Phenotype_UCDavis$p
Raw_Phenotype_UCDavis_f = Raw_Phenotype_UCDavis$f
Raw_Phenotype_UCDavis_e = Raw_Phenotype_UCDavis$e_cat_matrix
nrow(Raw_Phenotype_UCDavis_e)


Raw_Phenotype_UCDavis_f$label_procedure = sapply(strsplit(Raw_Phenotype_UCDavis_f$label,"___"),function(x){x[1]})
Raw_Phenotype_UCDavis_f$label_phenotype = sapply(strsplit(Raw_Phenotype_UCDavis_f$label,"___"),function(x){x[2]})

rownames(Raw_Phenotype_UCDavis_e) = Raw_Phenotype_UCDavis_f$label



# combine the phenotypes because one phenotype may have multiple procedures.
# length(Raw_Phenotype_UCDavis_f$label_phenotype)
unique_phenotype = unique(Raw_Phenotype_UCDavis_f$label_phenotype)
# length(unique_phenotype)

for(u in 1:length(unique_phenotype)){
  if(sum(Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u])>1){
    if(length(table(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],]))>1){
      duplicated_value_index = which(apply(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],1:2],2,function(x){
        if(length(unique(x)) == 1){
          if(!unique(x)=="NA"){
            return(TRUE)
          }else{
            return(FALSE)
          }
        }else{
          return(FALSE)
        }
      }))
      if(length(duplicated_value_index)>0){
        print(unique_phenotype[u])
        print(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],1:10])
        Sys.sleep(20)
      }
    }
  }
}




Raw_Phenotype_UCDavis_e_merge = matrix(NA, nrow = length(unique(Raw_Phenotype_UCDavis_f$label_phenotype)), ncol = ncol(Raw_Phenotype_UCDavis_e))
rownames(Raw_Phenotype_UCDavis_e_merge) = unique(Raw_Phenotype_UCDavis_f$label_phenotype)
# Raw_Phenotype_UCDavis_f_merge = Raw_Phenotype_UCDavis_f[which(unique(Raw_Phenotype_UCDavis_f$label_phenotype)]

i = 1
for(u in 1:length(unique_phenotype)){
  
  if(sum(Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u])>1){
    # stop(u)
    for(j in 1:ncol(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],])){
      the_j_th_col = Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],j]
      if(!length(unique(the_j_th_col)) == 1){
        Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],][the_j_th_col == "NA",j] = Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],][!the_j_th_col == "NA",j] 
      } 
    }
    Raw_Phenotype_UCDavis_e_merge[i,] = unique(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],])
    i = i + 1
  }else{
    Raw_Phenotype_UCDavis_e_merge[i,] = Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],]
    i = i + 1
  }
}


Raw_Phenotype_UCDavis_e_merge[Raw_Phenotype_UCDavis_e_merge == "NA"] = NA


Raw_Phenotype_UCDavis_e_null = Raw_Phenotype_UCDavis_e_merge[,Raw_Phenotype_UCDavis_p$Genotype %in% "null"]
Raw_Phenotype_UCDavis_p_null = Raw_Phenotype_UCDavis_p[Raw_Phenotype_UCDavis_p$Genotype %in% "null",]
Raw_Phenotype_UCDavis_e_null<-as.matrix(Raw_Phenotype_UCDavis_e_null)

Raw_Phenotype_UCDavis_e_null = apply(Raw_Phenotype_UCDavis_e_null,2,as.numeric)
rownames(Raw_Phenotype_UCDavis_e_null) = rownames(Raw_Phenotype_UCDavis_e_merge)




num_missing_male = apply(Raw_Phenotype_UCDavis_e_null,1,function(x){
  sum(is.na(x[Raw_Phenotype_UCDavis_p_null$Gender %in% "Male"]))
})
num_missing_female = apply(Raw_Phenotype_UCDavis_e_null,1,function(x){
  sum(is.na(x[Raw_Phenotype_UCDavis_p_null$Gender %in% "Female"]))
})
# names(num_missing_male) = names(num_missing_female) = Raw_Phenotype_UCDavis_f_merge$label


missing_index = ((num_missing_male/sum(Raw_Phenotype_UCDavis_p_null$Gender %in% "Male")) >= 0.7)&((num_missing_female/sum(Raw_Phenotype_UCDavis_p_null$Gender %in% "Female")) >= 0.7)
Phenotype_UCDavis_e_null = Raw_Phenotype_UCDavis_e_null[!missing_index,]
rownames(Phenotype_UCDavis_e_null) = rownames(Raw_Phenotype_UCDavis_e_merge)[!missing_index]




p_val3=fc3=p_val_adj3=Coef3=F_Count=M_Count=continuous_index=c()

for (i in 1:nrow(Phenotype_UCDavis_e_null)){
  
    continuous = sum(is.na(as.numeric(Phenotype_UCDavis_e_null[i,!is.na(Phenotype_UCDavis_e_null[i,])]))) == 0
    
    
    if(continuous){
      continuous_index[i] = TRUE
      
  
  outliers_nullF<-boxplot.stats(Phenotype_UCDavis_e_null[i,][Raw_Phenotype_UCDavis_p_null$Gender%in%'Female'])$out
  outliers_nullM<-boxplot.stats(Phenotype_UCDavis_e_null[i,][Raw_Phenotype_UCDavis_p_null$Gender%in%'Male'])$out
  
  if( length(outliers_nullF)>0 & length(outliers_nullM)>0){
    
    Label_outlier<-c(which(Phenotype_UCDavis_e_null[i,]%in%outliers_nullF),which(Phenotype_UCDavis_e_null[i,]%in%outliers_nullM[[1]]))
    Phenotype_UCDavis_e_Null_New<-Phenotype_UCDavis_e_null[i,][-Label_outlier]
    Phenotype_UCDavis_p_Null_New<-Raw_Phenotype_UCDavis_p_null[-Label_outlier]
    
  }else if(length(outliers_nullF)>0 & length(outliers_nullM)==0){
    
    Label_outlier<-c(which(Phenotype_UCDavis_e_null[i,]%in%outliers_nullF))
    Phenotype_UCDavis_e_Null_New<-Phenotype_UCDavis_e_null[i,][-Label_outlier]
    Phenotype_UCDavis_p_Null_New<-Raw_Phenotype_UCDavis_p_null[-Label_outlier]
    
  }else if(length(outliers_nullF)==0 & length(outliers_nullM)>0){
    
    Label_outlier<-c(which(Phenotype_UCDavis_e_null[i,]%in%outliers_nullM))
    Phenotype_UCDavis_e_Null_New<-Phenotype_UCDavis_e_null[i,][-Label_outlier]
    Phenotype_UCDavis_p_Null_New<-Raw_Phenotype_UCDavis_p_null[-Label_outlier]
    
  }else{
    Phenotype_UCDavis_e_Null_New<-Phenotype_UCDavis_e_null[i,]
    Phenotype_UCDavis_p_Null_New<-Raw_Phenotype_UCDavis_p_null
  }
  
  
  Phenotype_UCDavis_e_Null_Stats<-data.table(Phenotype_UCDavis_e_Null_New,Phenotype_UCDavis_p_Null_New$Gender)
  colnames(Phenotype_UCDavis_e_Null_Stats)<-c("PhenoVAlue","Gender")
  Phenotype_UCDavis_e_Null_Stats<-Phenotype_UCDavis_e_Null_Stats[!is.na(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue),]
  F_Count[i] <- sum(Phenotype_UCDavis_e_Null_Stats$Gender%in%"Female")      
  M_Count[i] <- sum(Phenotype_UCDavis_e_Null_Stats$Gender%in%"Male")
  
  
  
  p_val3[i]= summary(lm(Phenotype_UCDavis_e_Null_New ~ Phenotype_UCDavis_p_Null_New$Gender,na.action = 'na.exclude'))$coefficients[2,4]
  Coef3[i]= summary(lm(Phenotype_UCDavis_e_Null_New ~ Phenotype_UCDavis_p_Null_New$Gender ,na.action = 'na.exclude'))$coefficients[2,1]
  
  fc3[i]=  mean(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue[Phenotype_UCDavis_e_Null_Stats$Gender %in% "Male"])/mean(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue[Phenotype_UCDavis_e_Null_Stats$Gender %in% "Female"])
  
}else{
  continuous_index[i] = FALSE
  F_Count[i]=NA
  M_Count[i]=NA
  
  p_val3[i]=NA
  Coef3[i] =NA
  # fdr3[length(fdr3)+1] = NA
  fc3[i] =NA
}
    
}



p_val_adj3 = p.adjust(p_val3,method ="fdr")

sum(p_val3[continuous_index]<0.05,na.rm = TRUE)/sum(continuous_index)
# sum(fdr3[continuous_index]<0.05,na.rm = TRUE)/sum(continuous_index)

fwrite(data.table(label = rownames(Phenotype_UCDavis_e_null),F_Count=F_Count, M_Count=M_Count,p_value = p_val3, coefficient=Coef3, fold_change = fc3, adjusted_p_value = p.adjust(p_val3,'fdr')),"C:\\Users\\yzhan\\Documents\\Fiehn\\3_KOMP\\7_Sexual D_20190605\\Jennly Zhang KOMP SD_ver09\\Output\\3 Pheno WT_continuous_remove outlier_20200204.csv")




####################################################################################################


## 2b phenotype -  sex adjusted by body weight
# body_weight_all
# sum(Raw_Phenotype_UCDavis_p_null$label%in%names(body_weight_all))/length(Raw_Phenotype_UCDavis_p_null$label)

Phenotype_UCDavisSelect = wcmc::read_data("C:\\Users\\yzhan\\Documents\\Fiehn\\3_KOMP\\7_Sexual D_20190605\\Jennly Zhang KOMP SD_ver08\\Input\\Phenotype_UCDavisSelect.xlsx")
Phenotype_UCDavisSelect_p = Phenotype_UCDavisSelect$p
Phenotype_UCDavisSelect_f = Phenotype_UCDavisSelect$f
Phenotype_UCDavisSelect_e = Phenotype_UCDavisSelect$e_matrix

body_weight_all = Phenotype_UCDavisSelect_e[,Phenotype_UCDavisSelect_p[["phenotype"]] %in% "Body weight"]
names(body_weight_all) = Phenotype_UCDavisSelect_f$label

body_weight = body_weight_all[Raw_Phenotype_UCDavis_p_null$label]




p_val4=Coef4=fc4=p_val_adj4=F_Count=M_Count=c()

for (i in 1:nrow(Phenotype_UCDavis_e_null)){
  
  continuous = sum(is.na(as.numeric(Phenotype_UCDavis_e_null[i,!is.na(Phenotype_UCDavis_e_null[i,])]))) == 0
  
  
  if(continuous){
    continuous_index[i] = TRUE
    
    
    outliers_nullF<-boxplot.stats(Phenotype_UCDavis_e_null[i,][Raw_Phenotype_UCDavis_p_null$Gender%in%'Female'])$out
    outliers_nullM<-boxplot.stats(Phenotype_UCDavis_e_null[i,][Raw_Phenotype_UCDavis_p_null$Gender%in%'Male'])$out
    
    if( length(outliers_nullF)>0 & length(outliers_nullM)>0){
      
      Label_outlier<-c(which(Phenotype_UCDavis_e_null[i,]%in%outliers_nullF),which(Phenotype_UCDavis_e_null[i,]%in%outliers_nullM[[1]]))
      Phenotype_UCDavis_e_Null_New<-Phenotype_UCDavis_e_null[i,][-Label_outlier]
      Phenotype_UCDavis_p_Null_New<-Raw_Phenotype_UCDavis_p_null[-Label_outlier]
      
    }else if(length(outliers_nullF)>0 & length(outliers_nullM)==0){
      
      Label_outlier<-c(which(Phenotype_UCDavis_e_null[i,]%in%outliers_nullF))
      Phenotype_UCDavis_e_Null_New<-Phenotype_UCDavis_e_null[i,][-Label_outlier]
      Phenotype_UCDavis_p_Null_New<-Raw_Phenotype_UCDavis_p_null[-Label_outlier]
      
    }else if(length(outliers_nullF)==0 & length(outliers_nullM)>0){
      
      Label_outlier<-c(which(Phenotype_UCDavis_e_null[i,]%in%outliers_nullM))
      Phenotype_UCDavis_e_Null_New<-Phenotype_UCDavis_e_null[i,][-Label_outlier]
      Phenotype_UCDavis_p_Null_New<-Raw_Phenotype_UCDavis_p_null[-Label_outlier]
      
    }else{
      Phenotype_UCDavis_e_Null_New<-Phenotype_UCDavis_e_null[i,]
      Phenotype_UCDavis_p_Null_New<-Raw_Phenotype_UCDavis_p_null
    }
    
    
    body_weight_stats<-body_weight[names(body_weight)%in%Phenotype_UCDavis_p_Null_New$label]
    
    Phenotype_UCDavis_e_Null_Stats<-data.table(Phenotype_UCDavis_e_Null_New,Phenotype_UCDavis_p_Null_New$Gender,body_weight[names(body_weight)%in%Phenotype_UCDavis_p_Null_New$label])
    colnames(Phenotype_UCDavis_e_Null_Stats)<-c("PhenoVAlue","Gender","Weight")
    Phenotype_UCDavis_e_Null_Stats<-Phenotype_UCDavis_e_Null_Stats[!(is.na(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue)|is.na(Phenotype_UCDavis_e_Null_Stats$Weight)),]
    
    F_Count[i] <- sum(Phenotype_UCDavis_e_Null_Stats$Gender%in%"Female")      
    M_Count[i] <- sum(Phenotype_UCDavis_e_Null_Stats$Gender%in%"Male")
    
    
    p_val4[i]= summary(lm(Phenotype_UCDavis_e_Null_New ~ Phenotype_UCDavis_p_Null_New$Gender+body_weight_stats,na.action = 'na.exclude'))$coefficients[2,4]
    Coef4[i]= summary(lm(Phenotype_UCDavis_e_Null_New ~ Phenotype_UCDavis_p_Null_New$Gender+body_weight_stats ,na.action = 'na.exclude'))$coefficients[2,1]
    
    fc4[i]=  mean(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue[Phenotype_UCDavis_e_Null_Stats$Gender %in% "Male"])/mean(Phenotype_UCDavis_e_Null_Stats$PhenoVAlue[Phenotype_UCDavis_e_Null_Stats$Gender %in% "Female"])
    
    
    
  }else{
    continuous_index[i] = FALSE
    
    
    F_Count[i] <- NA     
    M_Count[i] <- NA
    p_val4[i]=NA
    Coef4[i] =NA
    # fdr3[length(fdr3)+1] = NA
    fc4[i] =NA
  }
  
}



p_val_adj4 = p.adjust(p_val4,method ="fdr")

sum(p_val4[continuous_index]<0.05,na.rm = TRUE)/sum(continuous_index)


fwrite(data.table(label = rownames(Phenotype_UCDavis_e_null),F_Count=F_Count, M_Count=M_Count,p_value = p_val4, Coefficient= Coef4,fold_change = fc4, adjusted_p_value = p.adjust(p_val4,'fdr')),"C:\\Users\\yzhan\\Documents\\Fiehn\\3_KOMP\\7_Sexual D_20190605\\Jennly Zhang KOMP SD_ver09\\Output\\3 Pheno WT_adjust byBW_continuous_remove outlier.csv")




############################################################################################################
############################################################################################################
  

# 6 in null, the correlation between each phenotype and each metabolite.

setwd("C:\\Users\\yzhan\\Documents\\Fiehn\\3_KOMP\\7_Sexual D_20190605\\Jennly Zhang KOMP SD_ver08\\Input\\ONLINE_data")
Null_Metabolomics = wcmc::read_data("Null Metabolomics_online.xlsx")

Null_Metabolomics_p = Null_Metabolomics$p
Null_Metabolomics_f = Null_Metabolomics$f
Null_Metabolomics_e = Null_Metabolomics$e_matrix

Null_Metabolomics_e[Null_Metabolomics_e==0]<-NA
Null_Metabolomics_e[126,1]



#### filter missing value > 80%
num_missing_male = apply(Null_Metabolomics_e,1,function(x){
  sum(is.na(x[Null_Metabolomics_p$Gender %in% "Male"]))
})
num_missing_female = apply(Null_Metabolomics_e,1,function(x){
  sum(is.na(x[Null_Metabolomics_p$Gender %in% "Female"]))
})
names(num_missing_male) = names(num_missing_female) = Null_Metabolomics_f$label

missing_index = ((num_missing_male/sum(Null_Metabolomics_p$Gender %in% "Male")) > 0.8)|((num_missing_female/sum(Null_Metabolomics_p$Gender %in% "Female")) > 0.8)
Null_Metabolomics_e = Null_Metabolomics_e[!missing_index,]
Null_Metabolomics_f = Null_Metabolomics_f[!missing_index,]
nrow(Null_Metabolomics_e)



#### replace NA value
for(i in 1:nrow(Null_Metabolomics_e)){
  Null_Metabolomics_e[i,is.na(Null_Metabolomics_e[i,])] = 0.5 * min(Null_Metabolomics_e[i,], na.rm = TRUE)
}
rownames(Null_Metabolomics_e) = Null_Metabolomics_f$label


Raw_Phenotype_UCDavis = wcmc::read_data("C:\\Users\\yzhan\\Documents\\Fiehn\\3_KOMP\\7_Sexual D_20190605\\Jennly Zhang KOMP SD_ver08\\Input\\Raw_Phenotype_UCDavis.xlsx")
Raw_Phenotype_UCDavis_p = Raw_Phenotype_UCDavis$p
Raw_Phenotype_UCDavis_f = Raw_Phenotype_UCDavis$f
Raw_Phenotype_UCDavis_e = Raw_Phenotype_UCDavis$e_cat_matrix



Raw_Phenotype_UCDavis_f$label_procedure = sapply(strsplit(Raw_Phenotype_UCDavis_f$label,"___"),function(x){x[1]})
Raw_Phenotype_UCDavis_f$label_phenotype = sapply(strsplit(Raw_Phenotype_UCDavis_f$label,"___"),function(x){x[2]})

rownames(Raw_Phenotype_UCDavis_e) = Raw_Phenotype_UCDavis_f$label

# combine the phenotypes because one phenotype may have multiple procedures.
unique_phenotype = unique(Raw_Phenotype_UCDavis_f$label_phenotype)
for(u in 1:length(unique_phenotype)){
  if(sum(Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u])>1){
    if(length(table(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],]))>1){
      duplicated_value_index = which(apply(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],],2,function(x){
        if(length(unique(x)) == 1){
          if(!unique(x)=="NA"){
            return(TRUE)
          }else{
            return(FALSE)
          }
        }else{
          return(FALSE)
        }
      }))
      if(length(duplicated_value_index)>0){
        print(unique_phenotype[u])
        print(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],1:10])
        Sys.sleep(20)
      }
    }
  }
}




Raw_Phenotype_UCDavis_e_merge = matrix(NA, nrow = length(unique(Raw_Phenotype_UCDavis_f$label_phenotype)), ncol = ncol(Raw_Phenotype_UCDavis_e))

colnames(Raw_Phenotype_UCDavis_e_merge) = colnames(Raw_Phenotype_UCDavis_e)

rownames(Raw_Phenotype_UCDavis_e_merge) = unique(Raw_Phenotype_UCDavis_f$label_phenotype)

i = 1
for(u in 1:length(unique_phenotype)){
  
  if(sum(Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u])>1){
    # stop(u)
    for(j in 1:ncol(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],])){
      the_j_th_col = Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],j]
      if(!length(unique(the_j_th_col)) == 1){
        Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],][the_j_th_col == "NA",j] = Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],][!the_j_th_col == "NA",j] 
      } 
    }
    Raw_Phenotype_UCDavis_e_merge[i,] = unique(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],])
    i = i + 1
  }else{
    Raw_Phenotype_UCDavis_e_merge[i,] = Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],]
    i = i + 1
  }
}


nrow(Raw_Phenotype_UCDavis_e_merge)
nrow(Raw_Phenotype_UCDavis_e)
# 
Raw_Phenotype_UCDavis_e_merge[Raw_Phenotype_UCDavis_e_merge == "NA"] = NA


Raw_Phenotype_UCDavis_e_null = Raw_Phenotype_UCDavis_e_merge[,Raw_Phenotype_UCDavis_p$Genotype %in% "null"]
Raw_Phenotype_UCDavis_p_null = Raw_Phenotype_UCDavis_p[Raw_Phenotype_UCDavis_p$Genotype %in% "null",]
rownames(Raw_Phenotype_UCDavis_e_null) = rownames(Raw_Phenotype_UCDavis_e_merge)


Raw_Phenotype_UCDavis_e_null = Raw_Phenotype_UCDavis_e_null[,colnames(Null_Metabolomics_e)]
Raw_Phenotype_UCDavis_e_null = apply(Raw_Phenotype_UCDavis_e_null,2,as.numeric)
rownames(Raw_Phenotype_UCDavis_e_null) = rownames(Raw_Phenotype_UCDavis_e_merge)

# all gender
cor = cor(t(Null_Metabolomics_e), t(Raw_Phenotype_UCDavis_e_null),use = "pairwise.complete.obs", method = "spearman")

rownames(cor) = rownames(Null_Metabolomics_e)
colnames(cor) = rownames(Raw_Phenotype_UCDavis_e_null)

cor_test = cor
for(i in 1:nrow(cor_test)){
  print(i)
  for(j in 1:ncol(cor_test)){
    if(is.na(cor[i,j])){
      cor_test[i,j] = NA
    }else{
      cor_test[i,j] = cor.test(Null_Metabolomics_e[i,], Raw_Phenotype_UCDavis_e_null[j,], method = "spearman")$p.value
    }
  }
}
write.csv(cor,paste0("9 in null, the correlation pheno-metabo_corr both sex_ver08.csv"))
write.csv(cor_test,paste0("9 in null, the correlation pheno-metabo_(corrp value both sex)_ver08.csv"))


# female only
female_labels = Raw_Phenotype_UCDavis_p$label[Raw_Phenotype_UCDavis_p$Gender %in% "Female"]
cor_female = cor(t(Null_Metabolomics_e[,colnames(Null_Metabolomics_e) %in% female_labels]), t(Raw_Phenotype_UCDavis_e_null[,colnames(Raw_Phenotype_UCDavis_e_null) %in% female_labels]),use = "pairwise.complete.obs", method = "spearman")

rownames(cor_female) = rownames(Null_Metabolomics_e[,colnames(Null_Metabolomics_e) %in% female_labels])
colnames(cor_female) = rownames(Raw_Phenotype_UCDavis_e_null[,colnames(Raw_Phenotype_UCDavis_e_null) %in% female_labels])

cor_test_female = cor_female
for(i in 1:nrow(cor_test_female)){
  print(i)
  for(j in 1:ncol(cor_test_female)){
    if(is.na(cor_female[i,j])){
      cor_test_female[i,j] = NA
    }else{
      cor_test_female[i,j] = cor.test(Null_Metabolomics_e[i,colnames(Null_Metabolomics_e) %in% female_labels], 
                                      Raw_Phenotype_UCDavis_e_null[j,colnames(Raw_Phenotype_UCDavis_e_null) %in% female_labels], method = "spearman")$p.value
    }
  }
}

write.csv(cor_female,paste0("9 in null female, the correlation pheno-metabo_(coeff female)_ver08.csv"))
write.csv(cor_test_female,paste0("9 in null female, the correlation pheno-metabo_(corrp value female)_ver08.csv"))


# male only
male_labels = Raw_Phenotype_UCDavis_p$label[Raw_Phenotype_UCDavis_p$Gender %in% "Male"]
cor_male = cor(t(Null_Metabolomics_e[,colnames(Null_Metabolomics_e) %in% male_labels]), t(Raw_Phenotype_UCDavis_e_null[,colnames(Raw_Phenotype_UCDavis_e_null) %in% male_labels]),use = "pairwise.complete.obs", method = "spearman")

rownames(cor_male) = rownames(Null_Metabolomics_e[,colnames(Null_Metabolomics_e) %in% male_labels])
colnames(cor_male) = rownames(Raw_Phenotype_UCDavis_e_null[,colnames(Raw_Phenotype_UCDavis_e_null) %in% male_labels])

cor_test_male = cor_male
for(i in 1:nrow(cor_test_male)){
  print(i)
  for(j in 1:ncol(cor_test_male)){
    if(is.na(cor_male[i,j])){
      cor_test_male[i,j] = NA
    }else{
      cor_test_male[i,j] =  cor_test_male[i,j] = cor.test(Null_Metabolomics_e[i,colnames(Null_Metabolomics_e) %in% male_labels], 
                                                            Raw_Phenotype_UCDavis_e_null[j,colnames(Raw_Phenotype_UCDavis_e_null) %in% male_labels], method = "spearman")$p.value
    }
  }
}

write.csv(cor_male,paste0("9 in null male, the correlation pheno-metabo_(coeff male)_ver08.csv"))
write.csv(cor_test_male,paste0("9 in null male, the correlation pheno-metabo_(corrp value male)_ver08.csv"))























































