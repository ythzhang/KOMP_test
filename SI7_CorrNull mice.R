# SI7 in null mice, the correlation between each phenotype and each metabolite.
# setweed
setwd("")

# read metabolomics dataset and extract the matrix
Null_Metabolomics = wcmc::read_data("Supplementary Data 1_raw metabolomics_v2_null.xlsx")

Null_Metabolomics_p = Null_Metabolomics$p
Null_Metabolomics_f = Null_Metabolomics$f
Null_Metabolomics_e = Null_Metabolomics$e_matrix

Null_Metabolomics_e[Null_Metabolomics_e==0] <- NA
Null_Metabolomics_e[Null_Metabolomics_e=="NA"] <- NA
Null_Metabolomics_e[Null_Metabolomics_e=="NaN"] <- NA
Null_Metabolomics_e[Null_Metabolomics_e=="Inf"] <- NA
Null_Metabolomics_e[!is.finite(Null_Metabolomics_e)] <- NA
sum(Null_Metabolomics_e=="Inf")

# Null_Metabolomics_e[100,1]
#### filter missing value >= 70%
num_missing_male = apply(Null_Metabolomics_e,1,function(x){
  sum(is.na(x[Null_Metabolomics_p$Gender %in% "Male"]))
})
num_missing_female = apply(Null_Metabolomics_e,1,function(x){
  sum(is.na(x[Null_Metabolomics_p$Gender %in% "Female"]))
})
names(num_missing_male) = names(num_missing_female) = Null_Metabolomics_f$label

missing_index = ((num_missing_male/sum(Null_Metabolomics_p$Gender %in% "Male")) >= 0.7) & ((num_missing_female/sum(Null_Metabolomics_p$Gender %in% "Female")) >= 0.7)
Null_Metabolomics_e = Null_Metabolomics_e[!missing_index,]
Null_Metabolomics_f = Null_Metabolomics_f[!missing_index,]
# nrow(Null_Metabolomics_e)
# write.table(Null_Metabolomics_f, "Null_Metabolomics_f.txt")
rownames(Null_Metabolomics_e) = Null_Metabolomics_f$label
#### replace NA value
# for(i in 1:nrow(Null_Metabolomics_e)){
#   Null_Metabolomics_e[i,is.na(Null_Metabolomics_e[i,])] = 0.5 * min(Null_Metabolomics_e[i,], na.rm = TRUE)
# }
# read phenotype data and extrax matrix
Phenotype_UCDavis = wcmc::read_data("Supplementary Data 2_raw phenotype_null.xlsx")
Raw_Phenotype_UCDavis_p_null = Phenotype_UCDavis$p
Raw_Phenotype_UCDavis_f_null = Phenotype_UCDavis$f
Raw_Phenotype_UCDavis_e_null = Phenotype_UCDavis$e_matrix

rownames(Raw_Phenotype_UCDavis_e_null) = Raw_Phenotype_UCDavis_f_null$label
Raw_Phenotype_UCDavis_e_null = apply(Raw_Phenotype_UCDavis_e_null,2,as.numeric)

# filter missingvalue >=70%
num_missing_male = apply(Raw_Phenotype_UCDavis_e_null,1,function(x){
  sum(is.na(x[Raw_Phenotype_UCDavis_p_null$Gender %in% "Male"]))
})
num_missing_female = apply(Raw_Phenotype_UCDavis_e_null,1,function(x){
  sum(is.na(x[Raw_Phenotype_UCDavis_p_null$Gender %in% "Female"]))
})
# names(num_missing_male) = names(num_missing_female) = Raw_Phenotype_UCDavis_f_merge$label

missing_index = ((num_missing_male/sum(Raw_Phenotype_UCDavis_p_null$Gender %in% "Male")) >= 0.7)&((num_missing_female/sum(Raw_Phenotype_UCDavis_p_null$Gender %in% "Female")) >= 0.7)
Phenotype_UCDavis_e_null = Raw_Phenotype_UCDavis_e_null[!missing_index,]
rownames(Phenotype_UCDavis_e_null) = rownames(Raw_Phenotype_UCDavis_e_null)[!missing_index]

# nrow(Phenotype_UCDavis_e_null)
# ==================================================================================================================================
# Correlation analysis, sexes combined
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

write.csv(cor,paste0("1_Null_correlation,Pheno-Metabo_(coeff both sex).csv"))
write.csv(cor_test,paste0("2_Null_correlation,Pheno-Metabo_(pvalue both sex).csv"))

# ============================================================================================
# correlation analysis, female only
female_labels = Raw_Phenotype_UCDavis_p_null$label[Raw_Phenotype_UCDavis_p_null$Gender %in% "Female"]
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
write.csv(cor_female,paste0("3_Null_correlation_female,Corr Pheno-Metabo_(coeff female).csv"))
write.csv(cor_test_female,paste0("4_Null_correlation_female,the Corr Pheno-Metabo_(pvalue female).csv"))

# correlation analysis, male only
male_labels = Raw_Phenotype_UCDavis_p_null$label[Raw_Phenotype_UCDavis_p_null$Gender %in% "Male"]
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
write.csv(cor_male,paste0("5_Null_correlation_male,Corr Pheno-Metabo_(coeff male).csv"))
write.csv(cor_test_male,paste0("6_Null_correlation_male,the Corr Pheno-Metabo_(pvalue male).csv"))

# ========================================================================================================================================================
# Recombine correlation results
# setwd("")

cor_bothsex<-read.csv("1_Null_correlation,Pheno-Metabo_(coeff both sex).csv",header = T )

pval_bothsex<-read.csv("2_Null_correlation,Pheno-Metabo_(pvalue both sex).csv",header = T)

cor_female<-read.csv("3_Null_correlation_female,Corr Pheno-Metabo_(coeff female).csv", header = T)

pval_female<-read.csv("4_Null_correlation_female,the Corr Pheno-Metabo_(pvalue female).csv",header = T)

cor_male<-read.csv("5_Null_correlation_male,Corr Pheno-Metabo_(coeff male).csv",header = T)

pval_male<-read.csv("6_Null_correlation_male,the Corr Pheno-Metabo_(pvalue male).csv",header = T)

# ===========================================================================================================================================
#report sample number for metabolomics data
Metabo_count_female = apply(Null_Metabolomics_e, 1, function(x){
  sum(!is.na(x[Null_Metabolomics_p$Gender %in% "Female"]))
}) 
Metabo_count_male = apply(Null_Metabolomics_e, 1, function(x){
  sum(!is.na(x[Null_Metabolomics_p$Gender %in% "Male"]))
}) 

#report sample number for phenotype data
Pheno_count_female = apply(Raw_Phenotype_UCDavis_e_null,1,function(x){
  sum(!is.na(x[Raw_Phenotype_UCDavis_p_null$Gender %in% "Female"]))
})
Pheno_count_male = apply(Raw_Phenotype_UCDavis_e_null,1,function(x){
  sum(!is.na(x[Raw_Phenotype_UCDavis_p_null$Gender %in% "Male"]))
})


# ========================================================================================
library(data.table)
# confirm that rownames and colnames of files match to each other, exclude row 1 and colume 1.
colmatch<-c()
for (i in 2:length(colnames(pval_female))){
  if(colnames(pval_female)[i]%in%colnames(pval_male)[i]){
    print('TRUE')
    colmatch[i] <- "true"
  }else{
    print('FALSE')
    colmatch[i] <- "FALSE" 
  }
}

library(data.table)
colmatch<-c()
for (i in 2:length(colnames(pval_female))){
  if(colnames(pval_female)[i]%in%colnames(pval_bothsex)[i]){
    print('TRUE')
    colmatch[i] <- "true"
  }else{
    print('FALSE')
    colmatch[i] <- "FALSE"
  }
}


# ==================================================================================================================================
#merge pval_both, pval_f, pval_m, cor_both, cor_f, cor_m based on each phenotype with all metabolites
list<-sprintf('%0.3d', 1:300)
colnames(pval_female)[1] <- 'Label'
colnames(pval_male)[1] <- 'Label'
colnames(pval_bothsex)[1] <- 'Label'

i=2
for (i in 2:ncol(pval_female)){
  
  Pheno_name <- colnames(pval_female)[i]
  
  Corr_pval_both <- pval_bothsex[,i]
  Corr_pval_f <- pval_female[,i]
  Corr_pval_m <- pval_male[,i]
  
  Corr_coeff_both <- cor_bothsex[,i]
  Corr_coeff_f <- cor_female[,i]
  Corr_coeff_m <- cor_male[,i]
  
  Metabolite_label <- pval_female$Label
  Assay <- Null_Metabolomics_f$Assay
  Metabolites <- Null_Metabolomics_f$CompoundName
  Phenotype_label <- Raw_Phenotype_UCDavis_f_null[i-1,]$Parameter
  Procedure_name <- Raw_Phenotype_UCDavis_f_null[i-1,]$`Procedure name`
  Phenotype_name <- Raw_Phenotype_UCDavis_f_null[i-1,]$`Parameter name`
  Pheno_current_count_female <- Pheno_count_female[i-1]
  Pheno_current_count_male <- Pheno_count_male[i-1]
  fwrite(data.table(Metabolite_label, Assay, Metabolites, Corr_coeff_both, Corr_coeff_f, Corr_coeff_m, Corr_pval_both, Corr_pval_f, Corr_pval_m, Phenotype_label, Procedure_name, Phenotype_name, Pheno_name, Metabo_count_female, Metabo_count_male, Pheno_current_count_female, Pheno_current_count_male), paste0("9_COMBINED correlation result (",list[i],") (",colnames(pval_female)[i],") in Null.csv"))
}


# ==================================================================================================================
### classify the correlation result based on the p value and cor value

files_cor <- list.files(pattern="^9_COMBINED correlation result (.*) in Null",full.names = T)
length(files_cor)
list<-sprintf('%0.3d', 1:300)
for(i in 1:length(files_cor)){
  corr <- read.csv(files_cor[i])
  corr$classify[(corr$Corr_pval_f < 0.05)|(corr$Corr_pval_m < 0.05)] <- 'significant'
  corr$classify[(corr$Corr_pval_f >= 0.05)&(corr$Corr_pval_m >= 0.05)] <- 'non-significant'
  
  corr$class <- 'Non-significant'
  corr$class[(corr$Corr_pval_f < 0.05 & corr$Corr_pval_m < 0.05) & (corr$Corr_coeff_f/corr$Corr_coeff_m < 0)] = 'Different direction'
  corr$class[(corr$Corr_pval_f < 0.05 & corr$Corr_pval_m < 0.05 )& (corr$Corr_coeff_f/corr$Corr_coeff_m > 0)] = 'Different size'
  corr$class[(corr$Corr_pval_f < 0.05 & corr$Corr_pval_m >= 0.05) | (corr$Corr_pval_f >= 0.05 & corr$Corr_pval_m < 0.05)] = 'One sex'
  corr$class[(corr$Corr_pval_f >= 0.05 & corr$Corr_pval_m >= 0.05) & (corr$Corr_coeff_f/corr$Corr_coeff_m > 0) & (corr$Corr_coeff_f/corr$Corr_coeff_both > 0) & corr$Corr_pval_both < 0.05] = 'Significant only if sexes combined'
  k = i + 1
  #corr1<-corr[corr$class=='significant',]
  write.csv(corr,paste0("9 Extracted correlation result (",list[k],") (",colnames(pval_female)[k],") significant in Null_both sex.csv"))
}


# -------------------------------------------------------------------------------------------------------------
###Summary classes for each phenotype
files_cor <- list.files(pattern="^9_COMBINED correlation result (.*) in Null", full.names = T)
length(files_cor)

#list<-sprintf('%0.3d', 1:300)
summ<-matrix(c('Metabolite label', 'Different direction', 'Different size', 'Significant only in female', 'Significant only in male', 'Significant only if sexes combined', 'Non-significant'), byrow = FALSE)
#summ[,1]=c('Differentdirection','Differentsize','Onesex','Nonsignificant')
for(i in 1:length(files_cor)){
  corr <- read.csv(files_cor[i])
  corr$classify[(corr$Corr_pval_f < 0.05) | (corr$Corr_pval_m < 0.05)] <- 'significant'
  corr$classify[(corr$Corr_pval_f >= 0.05) & (corr$Corr_pval_m >= 0.05)] <- 'non-significant'
  
  corr$class <-'Non-significant'
  corr$class[(corr$Corr_pval_f < 0.05 & corr$Corr_pval_m < 0.05) & (corr$Corr_coeff_f/corr$Corr_coeff_m < 0)] = 'Different direction'
  corr$class[(corr$Corr_pval_f < 0.05 & corr$Corr_pval_m < 0.05) & (corr$Corr_coeff_f/corr$Corr_coeff_m > 0)] = 'Different size'
  corr$class[corr$Corr_pval_f < 0.05 & corr$Corr_pval_m > 0.05] = 'Significant only in female'
  corr$class[corr$Corr_pval_f > 0.05 & corr$Corr_pval_m < 0.05] = 'Significant only in male'
  corr$class[(corr$Corr_pval_f >= 0.05 & corr$Corr_pval_m >= 0.05) & (corr$Corr_coeff_f/corr$Corr_coeff_m > 0) & (corr$Corr_coeff_f/corr$Corr_coeff_both > 0) & corr$Corr_pval_both < 0.05] = 'Significant only if sexes combined'
  
  Nonsignificant <- sum(corr$class == 'Non-significant')
  Differentdirection <- sum(corr$class == 'Different direction')
  Differentsize <- sum(corr$class == 'Different size')
  Sig_only_in_female <- sum(corr$class == 'Significant only in female')
  Sig_only_in_male <- sum(corr$class == 'Significant only in male')
  Significant_only_if_sexes_combined <- sum(corr$class == 'Significant only if sexes combined')
  
  k = i + 1
  a <- c(colnames(pval_female[k]), Differentdirection, Differentsize, Sig_only_in_female, Sig_only_in_male, Significant_only_if_sexes_combined, Nonsignificant)
  a <- as.matrix(a)
  summ <- cbind(summ, a)
}
#colnames(summ)=colnames(pval_female)
fwrite(data.table(summ), "9_Summary correlation in each class.csv")


# -----------------------------------------------------------------------------------------------------------
#### Merge all significant corr result in  phenotype and metabolite

files_cor <- list.files(pattern  ="^9_COMBINED correlation result (.*) in Null", full.names = T)
length(files_cor)
corrmerge <- matrix(nrow = 0, ncol = 17)
for(i in 1:length(files_cor)){
  corr <- read.csv(files_cor[i])
  corr$classify[(corr$Corr_pval_f < 0.05) | (corr$Corr_pval_m < 0.05)] <- 'significant'
  corr$classify[(corr$Corr_pval_f >= 0.05) & (corr$Corr_pval_m >= 0.05)] <- 'non-significant'
  
  corr$class <-'Non-significant'
  corr$class[(corr$Corr_pval_f < 0.05 & corr$Corr_pval_m < 0.05) & (corr$Corr_coeff_f/corr$Corr_coeff_m < 0)] = 'Different direction'
  corr$class[(corr$Corr_pval_f < 0.05 & corr$Corr_pval_m < 0.05) & (corr$Corr_coeff_f/corr$Corr_coeff_m > 0)] = 'Different size'
  corr$class[corr$Corr_pval_f < 0.05 & corr$Corr_pval_m > 0.05] = 'Significant only in female'
  corr$class[corr$Corr_pval_f > 0.05 & corr$Corr_pval_m < 0.05] = 'Significant only in male'
  corr$class[(corr$Corr_pval_f >= 0.05 & corr$Corr_pval_m >= 0.05) & (corr$Corr_coeff_f/corr$Corr_coeff_m > 0) & (corr$Corr_coeff_f/corr$Corr_coeff_both > 0) & corr$Corr_pval_both < 0.05] = 'Significant only if sexes combined'
  
  index <- corr$class == 'Non-significant'
  corr1 <- corr[!index, ]
  k = i + 1
  if(nrow(corr1) > 0){
    corr1$pheno <- colnames(pval_female[k])
  }
  
  corrmerge <- rbind(corrmerge,corr1)
  head(corrmerge)
  # corrmerge <- na.omit(corrmerge)
}

write.csv(corrmerge,paste0("9_Heatmap Merged Table corr significant in Null.csv"))