#######################################################################################################
# Statistics methods are adapted from : Supporting Data of Prevalence of Sexual Dimorphism in Mammalian Phenotypic Traits (2017).

#######################################################################################################
# Load packages
setwd("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
library(nlme)
library(RNOmni)
# for each genotype, check which metabolite is associated with gene * gender, gene, gender, or not significant at all.
# library(car)
# Function 1: build model for testing fixed effects
# Assumption: assume no batch effect and weight is not a independet variance 
model_forFIXEDtest<-function(dataset, Intensity){
  model.formula <- as.formula(paste(Intensity, "~", paste("genotype", "sex", "genotype:sex", sep= "+")))
  model=gls(model.formula, dataset, na.action="na.omit")
}
# -------------------------------------------------------------------------------------------------------------------------------
#Function 2: testing the fixed effects and building genotype model formula
# Goal: Genotype is always include in this model. The goal is to test genotype and sex effect.
# anova function: null hypothesis, the regression coefficients are equal to zero  
#                 alternative hypothesis that the regression coefficient are not equal to zero.
# the null hypothesis was rejected when p-values < 0.05 (e.g.,accept the alternative hypothesis, the components of the model should be included in later analysis)
# Note a complexity surrounds the interaction term  - if it is significant but gender is excluded it is excluded. 
final_genotype_model<-function(dataset, Intensity){
  model_afterFIXED=model_forFIXEDtest(dataset, Intensity)
  anova_results = anova(model_afterFIXED, type="marginal")$"p-value" < 0.05
  keepSex = anova_results[3]
  keepInteraction = anova_results[4]
  if(keepSex && keepInteraction){
    return(model.formula <- as.formula(paste(Intensity, "~", paste("genotype", "sex", "genotype:sex", sep= "+"))))
  }else if(keepSex &&  !keepInteraction){
    return(model.formula <- as.formula(paste(Intensity, "~", paste("genotype", "sex", sep= "+"))))
  }else if(!keepSex &&  !keepInteraction){
    return(model.formula <- as.formula(paste(Intensity, "~", "genotype")))
  }else if(!keepSex  && keepInteraction){
    return(model.formula <- as.formula(paste(Intensity, "~", paste("genotype", "sex", "genotype:sex", sep= "+"))))
  }
}
#---------------------------------------------------------------------------------------------------------------------------
#Function 3: testing the fixed effects and building final Genotype effect null model
# Goal:  to test fixed effects of the model and based on the output build the final null model formula for later testing - as a null model it automatically excludes genotype and interaction term. 
# The anova function tests the fixed effects associated by treatment with a null hypothesis that the regression coefficients are equal to zero and an alternative hypothesis that the regression coefficient are not equal to zero.
# If the p-values of these tests are less than 0.05 we reject the null and accept the alternative that the are significant components of the model and should be included. 
# If no terms are significant a model can be build with just an intercept element this is specified as  "model.formula <- as.formula(paste(depVariable, "~", "1"))"
null_model_genotype<-function(dataset, Intensity){ 
  model_afterFIXED=model_forFIXEDtest(dataset, Intensity)
  anova_results = anova(model_afterFIXED, type="marginal")$"p-value" < 0.05
  # anova_results = anova.lm(model_afterFIXED, type="marginal")[p-value] < 0.05
  keepSex = anova_results[3]
  keepInteraction = anova_results[4]     
  if(!keepSex && !keepInteraction){
    return(model.formula <- as.formula(paste(Intensity, "~", "1")))
  }else{
    return(model.formula <- as.formula(paste(Intensity, "~", "sex")))
  } 
}
#----------------------------------------------------------------------------------
#Function 3-b: testing the fixed effects and building final Interaction effect null model
null_model_Interaction<-function(dataset, Intensity){ 
  model_afterFIXED=model_forFIXEDtest(dataset, Intensity)
  anova_results = anova(model_afterFIXED, type="marginal")$"p-value" < 0.05
  keepSex = anova_results[3]
  keepGenotype = anova_results[2]     
  if(!keepSex && !keepGenotype){
    return(model.formula <- as.formula(paste(Intensity, "~", "1")))
  }else{
    return(model.formula <- as.formula(paste(Intensity, "~", paste("genotype","sex",sep="+"))))
  } 
}
# ------------------------------------------------------------------------------------------------
#Function 4: testing the genotype effect
# The goal of the function is to give a pvalue of whether genotype effect is significant by compare the genotype model with the null model with the anova function
# The model_formula_null and model_formula_genotype are called to define the models for comparison.
# For each possible combination,  then an anova model is used to report the pvalue.
# For testing the genotype effect we use method= ML
testing_genotype_effect<-function(dataset,Intensity){
  #call functions to determine model.formula
  model_formula_null = null_model_genotype(dataset,Intensity)
  model_formula_genotype = final_genotype_model(dataset,Intensity)
  model_genotype=gls(model_formula_genotype,dataset,method='ML', na.action="na.omit")
  model_null=gls(model_formula_null, dataset,method='ML',na.action="na.omit")
  return(pvalue_genotype=(anova(model_genotype,model_null)$`p-value`[2]))
}

# ------------------------------------------------------------------------------------------------
#Function 4-b: testing the interaction effect
# The goal of the function is to give a pvalue of whether genotype effect is significant by compare the genotype model with the null model with the anova function
# The model_formula_null and model_formula_genotype are called to define the models for comparison.
# For each possible combination,  then an anova model is used to report the pvalue.
# For testing the genotype effect we use method= ML
testing_Interaction_effect<-function(dataset, Intensity){
  #call functions to determine model.formula
  formula_interaction_null = null_model_Interaction(dataset,Intensity)
  model_interaction_full=gls(Intensity~genotype+sex+genotype:sex,dataset, method='ML',na.action="na.omit")
  model_null_interaction =gls(formula_interaction_null, dataset,method='ML', na.action="na.omit")
  return(pvalue_Stage2=(anova(model_interaction_full, model_null_interaction)$`p-value`[2]))
}

# ------------------------------------------------------------------------------------------------
#Function 4-d: testing the sex effect
null_model_sex<-function(dataset, Intensity){ 
  model_afterFIXED=model_forFIXEDtest(dataset, Intensity)
  anova_results = anova(model_afterFIXED, type="marginal")$"p-value" < 0.05
  keepInteraction = anova_results[4]
  keepGenotype = anova_results[2]     
  if(!keepInteraction && !keepGenotype){
    return(model.formula <- as.formula(paste(Intensity, "~", "1")))
  }else{
    return(model.formula <- as.formula(paste(Intensity, "~", "genotype")))
  } 
}
# ------------------------------------------------------------------------------------------------
final_sex_model<-function(dataset, Intensity){
  model_afterFIXED=model_forFIXEDtest(dataset, Intensity)
  anova_results = anova(model_afterFIXED, type="marginal")$"p-value" < 0.05
  keepGenotype = anova_results[2]
  keepInteraction = anova_results[4]
  if( keepGenotype && keepInteraction){
    return(model.formula <- as.formula(paste(Intensity, "~", paste("genotype", "sex", "genotype:sex", sep= "+"))))
  }else if(keepGenotype &&  !keepInteraction){
    return(model.formula <- as.formula(paste(Intensity, "~", paste("genotype", "sex", sep= "+"))))
  }else if(!keepGenotype &&  !keepInteraction){
    return(model.formula <- as.formula(paste(Intensity, "~", "sex" )))
  }else if(!keepGenotype  && keepInteraction){
    return(model.formula <- as.formula(paste(Intensity, "~", paste("genotype", "sex", "genotype:sex", sep= "+"))))
  }
}

testing_sex_effect<-function(dataset, Intensity){
  #call functions to determine model.formula
  formula_sex_null = null_model_sex(dataset, Intensity)
  formula_sex_full = final_sex_model(dataset,Intensity)
  model_null_sex =gls(formula_sex_null, dataset,method='ML', na.action="na.omit")
  model_full_sex=gls(formula_sex_full, dataset, method='ML',na.action="na.omit")
  return(pvalue_Stage1.5=(anova(model_full_sex, model_null_sex)$`p-value`[2]))
}
# ---------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------
#######################################################################################################
#load statistical output from large scale assessment of role of sex across IMPC data
setwd("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
All_Metabolomics = wcmc::read_data("XXXXXXXXXXXXXXXXXXXXXXXXXX.xlsx")
All_Metabolomics_p = All_Metabolomics$p ##(first row plus label row)
All_Metabolomics_f = All_Metabolomics$f  ##(first column plus label column)
All_Metabolomics_e = All_Metabolomics$e_matrix ###(just matrix values)
# count the total genotypes
unique_genes = unique(All_Metabolomics_p$Genotype)
# replace any non numeric text with NA
All_Metabolomics_e_no_mising = All_Metabolomics_e
All_Metabolomics_e[!is.finite(All_Metabolomics_e)] <- NA
All_Metabolomics_e[All_Metabolomics_e ==0] <- NA
All_Metabolomics_e[All_Metabolomics_e =="NA"] <- NA
All_Metabolomics_e[All_Metabolomics_e =="NaN"] <- NA
All_Metabolomics_e[All_Metabolomics_e =="Inf"] <- NA
sum(All_Metabolomics_e == "Inf")
# --------------------------------------------------------------------------------------------------------------------------------------------
# missing values replacement
for(g in 1:length(unique_genes)){
  #read gene
  current_gene = unique_genes[g]
  #subset data for current genotype
  current_label = All_Metabolomics_p$label[All_Metabolomics_p$Genotype %in% c(current_gene)]
  current_e = All_Metabolomics_e[, current_label]
  #subset data for female mice in current genotype 
  current_label_F = All_Metabolomics_p$label[All_Metabolomics_p$Genotype %in% c(current_gene)&All_Metabolomics_p$Gender %in% 'Female']
  current_e_F = All_Metabolomics_e[, current_label_F]
  #subset data for male mice in current genotype 
  current_label_M = All_Metabolomics_p$label[All_Metabolomics_p$Genotype %in% c(current_gene)&All_Metabolomics_p$Gender %in% 'Male']
  current_e_M =All_Metabolomics_e[,current_label_M]
  # replace missing values in WT mice if missing values > 70%
  if (current_gene=='null'){
    for(i in 1 : nrow(current_e)){
      #Replace female NA
      if((sum(is.na(current_e_F[i, ])) >= 14) & (sum(is.na(current_e_M[i, ])) >= 14)){
        current_e_F[i, ] = NA
        current_e_M[i, ] = NA
      }else if((sum(is.na(current_e_F[i, ])) >= 14) & (sum(is.na(current_e_M[i,])) < 14)){
        current_e_F[i, is.na(current_e_F[i, ])] = runif(sum(is.na(current_e_F[i,])), min = 0.49, max = 0.51) * min(current_e[i, !is.na(current_e[i, ])])
        current_e_M[i, ] = current_e_M[i, ]
      }else if((sum(is.na(current_e_F[i, ])) < 14) & (sum(is.na(current_e_M[i, ])) >= 14)){
        current_e_F[i, ] = current_e_F[i, ]
        current_e_M[i, is.na(current_e_M[i, ])] = runif(sum(is.na(current_e_M[i, ])), min = 0.49, max = 0.51) * min(current_e[i, !is.na(current_e[i, ])])
      }else{
        current_e_F[i, ] = current_e_F[i, ]
        current_e_M[i, ] = current_e_M[i, ]
      }
    }
  }else{
    # replace missing values in each genotype if missing values > 70%
    for(i in 1 : nrow(current_e)){
      if(sum(is.na(current_e[i, ])) >= 4){
        current_e_F[i, ] = NA
        current_e_M[i, ] = NA
      }else if(sum(is.na(current_e_F[i, ])) == 3 & sum(is.na(current_e_M[i,])) == 0){
        current_e_F[i, is.na(current_e_F[i, ])] = runif(sum(is.na(current_e_F[i, ])), min = 0.49, max = 0.51) * min(current_e[i, !is.na(current_e[i, ])])
        current_e_M[i, ] = current_e_M[i, ]
      }else if(sum(is.na(current_e_F[i, ])) == 0 & sum(is.na(current_e_M[i, ])) == 3){
        current_e_F[i, ] = current_e_F[i, ]
        current_e_M[i, is.na(current_e_M[i, ])] = runif(sum(is.na(current_e_M[i, ])), min = 0.49, max = 0.51) * min(current_e[i, !is.na(current_e[i, ])])
      }else if(sum(is.na(current_e_F[i, ])) == 2 & sum(is.na(current_e_M[i, ])) == 1){
        current_e_F[i, is.na(current_e_F[i, ])] = runif(sum(is.na(current_e_F[i, ])), min = 0.49, max = 0.51) * min(current_e[i, !is.na(current_e[i, ])])
        current_e_M[i, is.na(current_e_M[i, ])] = runif(sum(is.na(current_e_M[i, ])), min=0.98, max = 1.02) * min(current_e[i, !is.na(current_e[i, ])])
      }else if(sum(is.na(current_e_F[i, ])) == 1 & sum(is.na(current_e_M[i,])) == 2){
        current_e_F[i, is.na(current_e_F[i, ])] = runif(sum(is.na(current_e_F[i, ])), min = 0.98, max = 1.02) * min(current_e[i, !is.na(current_e[i, ])])
        current_e_M[i, is.na(current_e_M[i, ])] = runif(sum(is.na(current_e_M[i, ])), min = 0.49, max = 0.51) * min(current_e[i, !is.na(current_e[i, ])])
      }else if(sum(is.na(current_e[i, ])) == 2){
        if(sum(is.na(current_e_F[i, ])) == 2 & sum(is.na(current_e_M[i, ])) == 0){
          current_e_F[i, is.na(current_e_F[i, ])] = runif(sum(is.na(current_e_F[i, ])), min = 0.98, max = 1.02) * min(current_e[i, !is.na(current_e[i, ])])
          current_e_M[i, ] = current_e_M[i, ]
        }else if(sum(is.na(current_e_F[i, ])) == 0 & sum(is.na(current_e_M[i, ])) == 2){
          current_e_F[i, ] = current_e_F[i, ]
          current_e_M[i, is.na(current_e_M[i, ])] = runif(sum(is.na(current_e_M[i, ])),min = 0.98, max = 1.02) * min(current_e[i, !is.na(current_e[i, ])])
        }else if(sum(is.na(current_e_F[i, ])) == 1 & sum(is.na(current_e_M[i, ])) == 1){
          current_e_F[i, is.na(current_e_F[i, ])] = runif(sum(is.na(current_e_F[i, ])), min = 0.98, max = 1.02) * min(current_e[i, !is.na(current_e[i, ])])
          current_e_M[i, is.na(current_e_M[i, ])] = runif(sum(is.na(current_e_M[i, ])), min = 0.98, max = 1.02) * min(current_e[i, !is.na(current_e[i, ])])
        }
      }else if(sum(is.na(current_e[i, ])) <= 1){
        current_e_F[i, is.na(current_e_F[i, ])] = runif(sum(is.na(current_e_F[i, ])), min = 0.98, max = 1.02) * min(current_e[i, !is.na(current_e[i, ])])
        current_e_M[i, is.na(current_e_M[i, ])] = runif(sum(is.na(current_e_M[i, ])), min = 0.98, max = 1.02) * min(current_e[i, !is.na(current_e[i, ])])
      }
    }
  }
  All_Metabolomics_e_no_mising[,current_label_F] = current_e_F
  All_Metabolomics_e_no_mising[,current_label_M] = current_e_M
}

# ---------------------------------------------------------
# Testing and deal with missing values
for(i in 1: nrow(All_Metabolomics_e_no_mising)){
  All_Metabolomics_f$CompoundName[i]
  # ---------------------------------------------------------
  # FEMALE
  outliers_F_raw_null <- boxplot.stats(All_Metabolomics_e_no_mising[i, All_Metabolomics_p$Gender %in% "Female" & All_Metabolomics_p$Genotype %in% "null"])$out
  outliers_F_label_null <- All_Metabolomics_p$label[which(All_Metabolomics_p$Gender == "Female" & All_Metabolomics_p$Genotype %in% "null" & All_Metabolomics_e_no_mising[i,] %in% outliers_F_raw_null)]
  outliers_F_genotype_null <- All_Metabolomics_p$Genotype[which(All_Metabolomics_p$Gender == "Female" & All_Metabolomics_p$Genotype %in% "null" & All_Metabolomics_e_no_mising[i,] %in% outliers_F_raw_null)]
  if (length(outliers_F_raw_null) == 0){
    All_Metabolomics_e_no_mising[i,All_Metabolomics_p$Gender %in% "Female"  & All_Metabolomics_p$Genotype %in% "null"] <- All_Metabolomics_e_no_mising[i,All_Metabolomics_p$Gender %in% "Female" & All_Metabolomics_p$Genotype %in% "null"]
  }else{
    if(length(outliers_F_raw_null) > (sum(!is.na(All_Metabolomics_e_no_mising[i, All_Metabolomics_p$Gender %in% "Female" & All_Metabolomics_p$Genotype %in% "null"])) /6)){
      outliers_F_raw_null <- NA
      outliers_F_genotype_null <- NA
    }else{
      outliers_F_raw_null <- outliers_F_raw_null
      outliers_F_label_null <- outliers_F_label_null
      outliers_F_genotype_null <- outliers_F_genotype_null
    }
    All_Metabolomics_e_no_mising[i, All_Metabolomics_p$label %in% outliers_F_label_null] <- NA
  }
  outliers_F_raw_KO <- boxplot.stats(All_Metabolomics_e_no_mising[i, All_Metabolomics_p$Gender %in% "Female"])$out
  outliers_F_label_KO <- All_Metabolomics_p$label[which(All_Metabolomics_p$Gender == "Female" & All_Metabolomics_e_no_mising[i,] %in% outliers_F_raw_KO)]
  outliers_F_genotype_KO <- All_Metabolomics_p$Genotype[which(All_Metabolomics_p$label %in% outliers_F_label_KO)]
  ref_F_gene <- unique(outliers_F_genotype_KO)
  if(length(outliers_F_raw_KO) == 0){
    All_Metabolomics_e_no_mising[i,All_Metabolomics_p$Gender %in% "Female"] <- All_Metabolomics_e_no_mising[i,All_Metabolomics_p$Gender %in% "Female"]
  }else{
    for (l in 1: length(ref_F_gene)){
      if (ref_F_gene[l] == "null"){
        null_kp_F_lable <- which(outliers_F_genotype_KO == "null")
        outliers_F_raw_KO <- outliers_F_raw_KO[-null_kp_F_lable]
        outliers_F_genotype_KO <- outliers_F_genotype_KO[-null_kp_F_lable]
        outliers_F_label_KO <- outliers_F_label_KO[-null_kp_F_lable]
      }else{
        if(sum(outliers_F_genotype_KO %in% ref_F_gene[l]) > 1 ){
          outliers_F_label_kp <- outliers_F_label_KO[outliers_F_genotype_KO %in% ref_F_gene[l]]
          outliers_F_label_ad <- All_Metabolomics_p$label[All_Metabolomics_p$Gender == "Female" & All_Metabolomics_p$Genotype %in% ref_F_gene[l] & (!All_Metabolomics_p$label %in% outliers_F_label_kp)]
          outliers_F_Intensity_ad <- All_Metabolomics_e_no_mising[i,All_Metabolomics_p$label %in% outliers_F_label_ad]
          KO_F_kp <- which(outliers_F_genotype_KO %in% ref_F_gene[l])
          outliers_F_raw_KO <- outliers_F_raw_KO[-KO_F_kp]
          outliers_F_raw_KO <- c(outliers_F_Intensity_ad, outliers_F_raw_KO)
          outliers_F_genotype_KO <- outliers_F_genotype_KO[-KO_F_kp]
          outliers_F_genotype_KO <- c(ref_F_gene[l], outliers_F_genotype_KO)
          outliers_F_label_KO <- outliers_F_label_KO[-KO_F_kp]
          outliers_F_label_KO <- c(outliers_F_label_ad ,outliers_F_label_KO)
        }else{
          outliers_F_raw_KO <- outliers_F_raw_KO
          outliers_F_genotype_KO <- outliers_F_genotype_KO
          outliers_F_label_KO <- outliers_F_label_KO
        }
      }
    }
    All_Metabolomics_e_no_mising[i, which(All_Metabolomics_p$label %in% outliers_F_label_KO) ] <- NA
  }
  # ---------------------------------------------------------------------------------------------------
  # MALE
  outliers_M_raw_null <- boxplot.stats(All_Metabolomics_e_no_mising[i, All_Metabolomics_p$Gender %in% "Male" & All_Metabolomics_p$Genotype %in% "null"])$out
  outliers_M_label_null <- All_Metabolomics_p$label[which(All_Metabolomics_p$Gender == "Male" & All_Metabolomics_p$Genotype %in% "null" & All_Metabolomics_e_no_mising[i,] %in% outliers_M_raw_null)]
  outliers_M_genotype_null <- All_Metabolomics_p$Genotype[which(All_Metabolomics_p$Gender == "Male" & All_Metabolomics_p$Genotype %in% "null" & All_Metabolomics_e_no_mising[i,] %in% outliers_M_raw_null)]
  
  if (length(outliers_M_raw_null) == 0){
    All_Metabolomics_e_no_mising[i,All_Metabolomics_p$Gender %in% "Male" & All_Metabolomics_p$Genotype %in% "null"] <- All_Metabolomics_e_no_mising[i,All_Metabolomics_p$Gender %in% "Male" & All_Metabolomics_p$Genotype %in% "null"]
  }else{
    if(length(outliers_M_raw_null) > (sum(!is.na(All_Metabolomics_e_no_mising[i, All_Metabolomics_p$Gender %in% "Male" & All_Metabolomics_p$Genotype %in% "null"])) /6)){
      outliers_M_raw_null <- NA
      outliers_M_genotype_null <- NA
      outliers_M_genotype_null <- NA
    }else{
      outliers_M_raw_null <- outliers_M_raw_null
      outliers_M_label_null <- outliers_M_label_null
      outliers_M_genotype_null <- outliers_M_genotype_null
    }
    All_Metabolomics_e_no_mising[i, All_Metabolomics_p$label %in% outliers_M_label_null] <- NA
  }
  # ------------------------------------------------------------------
  outliers_M_raw_KO <- boxplot.stats(All_Metabolomics_e_no_mising[i, All_Metabolomics_p$Gender %in% "Male"])$out
  outliers_M_label_KO <- All_Metabolomics_p$label[which(All_Metabolomics_p$Gender == "Male" & All_Metabolomics_e_no_mising[i,] %in% outliers_M_raw_KO)]
  outliers_M_genotype_KO <- All_Metabolomics_p$Genotype[which(All_Metabolomics_p$label %in% outliers_M_label_KO)]
  ref_M_gene <- unique(outliers_M_genotype_KO)
  if(length(outliers_M_raw_KO)==0){
    All_Metabolomics_e_no_mising[i, All_Metabolomics_p$Gender %in% "Male"] <- All_Metabolomics_e_no_mising[i, All_Metabolomics_p$Gender %in% "Male"]
  }else{
    for (m in 1: length(ref_M_gene)){
      if (ref_M_gene[m] == "null"){
        null_kp_M_lable <- which(outliers_M_genotype_KO == "null")
        outliers_M_raw_KO <- outliers_M_raw_KO[-null_kp_M_lable]
        outliers_M_genotype_KO <- outliers_M_genotype_KO[-null_kp_M_lable]
        outliers_M_label_KO <- outliers_M_label_KO[-null_kp_M_lable]
      }else{
        if(sum(outliers_M_genotype_KO %in% ref_M_gene[m]) > 1 ){
          outliers_M_label_kp <- outliers_M_label_KO[outliers_M_genotype_KO %in% ref_M_gene[m]]
          outliers_M_label_ad <- All_Metabolomics_p$label[All_Metabolomics_p$Gender == "Male" & All_Metabolomics_p$Genotype %in% ref_M_gene[m] & (!All_Metabolomics_p$label %in% outliers_M_label_kp)]
          outliers_M_Intensity_ad <- All_Metabolomics_e_no_mising[i,All_Metabolomics_p$label %in% outliers_M_label_ad]
          KO_M_kp <- which(outliers_M_genotype_KO %in% ref_M_gene[m])
          outliers_M_raw_KO <- outliers_M_raw_KO[-KO_M_kp]
          outliers_M_raw_KO <- c(outliers_M_Intensity_ad, outliers_M_raw_KO)
          outliers_M_genotype_KO <- outliers_M_genotype_KO[-KO_M_kp]
          outliers_M_genotype_KO <- c(ref_M_gene[m], outliers_M_genotype_KO)
          outliers_M_label_KO <- outliers_M_label_KO[-KO_M_kp]
          outliers_M_label_KO <- c(outliers_M_label_ad ,outliers_M_label_KO)
        }else{
          outliers_M_raw_KO <- outliers_M_raw_KO
          outliers_M_genotype_KO <- outliers_M_genotype_KO
          outliers_M_label_KO <- outliers_M_label_KO
        }
      }
    }
    All_Metabolomics_e_no_mising[i, which(All_Metabolomics_p$label %in% outliers_M_label_KO) ] <- NA
  }
}

#subset data by sex for wildtype mice
null_label = All_Metabolomics_p$label[All_Metabolomics_p$Genotype %in% c("null")]
null_e = All_Metabolomics_e_no_mising[, null_label]
null_label = All_Metabolomics_p$label[All_Metabolomics_p$Genotype %in% c("null") & All_Metabolomics_p$Gender %in% "Female"]
null_e_female = All_Metabolomics_e_no_mising[, null_label]
null_label = All_Metabolomics_p$label[All_Metabolomics_p$Genotype %in% c("null") & All_Metabolomics_p$Gender %in% "Male"]
null_e_male = All_Metabolomics_e_no_mising[, null_label]
# ------------------------------------------------------------------------------------------------------------------------------------------
# Starting TwoStage significance assessment using Two-way ANOVA
Normality_before = Normality_F = Normality_M = list()
numb_null = numb_gene = list()
pvalue_Stage1 = pvalue_Stage1b = pvalue_Stage2 = list()
sep_allKO_pval = sep_allKO_estimate =list()
sep_FvKO_pval = sep_FvKO_estimate = sep_MvKO_pval = sep_MvKO_estimate = list()
foldchange_FvKO = foldchange_MvKO = list()
ctrl_sex_pval = ctrl_sex_estimate  = list() 
KO_sex_pval = KO_sex_estimate = list()
# statistical analysis for each genotype
for(g in 2:length(unique_genes)){
  # subset data for each genotype
  current_gene = unique_genes[g]
  # print(g)
  current_label <- All_Metabolomics_p$label[All_Metabolomics_p$Genotype %in% c(current_gene)]
  current_e <- All_Metabolomics_e_no_mising[, current_label]
  current_p <- All_Metabolomics_p[All_Metabolomics_p$label %in% current_label, ]
  # subset data for WT mice
  # null_p <- All_Metabolomics_p[All_Metabolomics_p$Genotype %in% 'null', ]
  Normality_before[[current_gene]] = Normality_F[[current_gene]] = Normality_M[[current_gene]] = c()
  numb_null[[current_gene]] = numb_gene[[current_gene]] = c()
  pvalue_Stage1[[current_gene]] = pvalue_Stage1b[[current_gene]] = pvalue_Stage2[[current_gene]] = c()
  sep_allKO_pval[[current_gene]] = sep_allKO_estimate[[current_gene]] = c() 
  sep_FvKO_pval[[current_gene]] = sep_FvKO_estimate[[current_gene]] = sep_MvKO_pval[[current_gene]] = sep_MvKO_estimate[[current_gene]] = c()
  foldchange_FvKO[[current_gene]] = foldchange_MvKO[[current_gene]] = c()
  ctrl_sex_pval[[current_gene]] = ctrl_sex_estimate[[current_gene]] = c() 
  KO_sex_pval[[current_gene]] = KO_sex_estimate[[current_gene]] = c()
  # statistical test for each metabolite variable for each genotyope
  for(i in 1 : nrow(All_Metabolomics_e_no_mising)){
    # if missing values > 70%, then not significant 
    if((sum(is.na(current_e[i,])) >= (ncol(current_e)-2) ) | (sum(is.na(null_e[i,])) >= (ncol(null_e)-6) )){
      numb_null[[current_gene]][i] = numb_gene[[current_gene]][i] <- NA
      Normality_before[[current_gene]][i] <- NA
      Normality_F[[current_gene]][i] <- NA
      Normality_M[[current_gene]][i] <- NA
      #p-values not available, NA
      pvalue_Stage1[[current_gene]][i] <- NA
      pvalue_Stage1b[[current_gene]][i] <- NA
      pvalue_Stage2[[current_gene]][i] <- NA
      ctrl_sex_pval[[current_gene]][i] = ctrl_sex_estimate[[current_gene]][i] <- NA
      KO_sex_pval[[current_gene]][i] = KO_sex_estimate[[current_gene]][i] <- NA
      sep_allKO_pval[[current_gene]][i] = sep_allKO_estimate[[current_gene]][i] <- NA
      sep_FvKO_pval[[current_gene]][i] = sep_FvKO_estimate[[current_gene]][i] <- NA
      sep_MvKO_pval[[current_gene]][i] = sep_MvKO_estimate[[current_gene]][i] <- NA
      #foldchange NA
      foldchange_FvKO[[current_gene]][i] <- NA
      foldchange_MvKO[[current_gene]][i] <- NA
    }else{
      # ==================================================================================
      # do statistics for metabolite variables that pass missing values criteria
      data_raw <- data.table(y = c(current_e[i, ], null_e[i, ]), group = rep(c("Gene", "Control"), c(ncol(current_e), ncol(null_e))), sex = c(current_p$Gender, null_p$Gender))
      colnames(data_raw) <- c("Intensity_raw", "genotype", 'sex')
      ### need to remove data for which metabolite value was not available
      data_raw <- data_raw[!is.na(data_raw$Intensity_raw), ]
      # test normality before data tranformation
      Normality_before[[current_gene]][i] <- shapiro.test(data_raw$Intensity_raw)$p.value
      # tranform data
      data_raw$Intensity_trans <- rankNorm(data_raw$Intensity_raw)
      dataset <- data_raw[!is.na(data_raw$Intensity), ]
      # =========================================================================
      #statistical analysis
      model_afterFIXED=tryCatch({
        model_forFIXEDtest(dataset, "Intensity")
      }, error = function(er){
        NA
      })
      if(length(model_afterFIXED) <= 1){
        #count number for WT mice and KO mice
        numb_null[[current_gene]][i] <- sum(dataset$genotype %in% "Control")
        numb_gene[[current_gene]][i] <- sum(dataset$genotype %in% "Gene")
        #normality test is NA
        Normality_before[[current_gene]][i] <- NA
        Normality_F[[current_gene]][i] <- NA
        Normality_M[[current_gene]][i] <- NA
        #p-values NA
        pvalue_Stage1[[current_gene]][i] <- NA
        pvalue_Stage1b[[current_gene]][i] <- NA
        pvalue_Stage2[[current_gene]][i] <- NA
        ctrl_sex_pval[[current_gene]][i] = ctrl_sex_estimate[[current_gene]][i] <- NA 
        KO_sex_pval[[current_gene]][i] = KO_sex_estimate[[current_gene]][i] <- NA
        sep_allKO_pval[[current_gene]][i] = sep_allKO_estimate[[current_gene]][i] <- NA
        sep_FvKO_pval[[current_gene]][i] = sep_FvKO_estimate[[current_gene]][i] <- NA
        sep_MvKO_pval[[current_gene]][i] = sep_MvKO_estimate [[current_gene]][i] <- NA
        
        # fold-change calculation for female KO mice to female WT mice
        foldchange_FvKO[[current_gene]][i] <- mean(dataset$Intensity[dataset$genotype %in% "Gene" & dataset$sex %in% "Female"], na.rm = T)/ mean(dataset$Intensity[dataset$genotype %in% "Control" & dataset$sex %in% "Female"], na.rm = T)
        # fold-change calculation for male KO mice to male WT mice
        foldchange_MvKO[[current_gene]][i] <- mean(dataset$Intensity[dataset$genotype %in% "Gene" & dataset$sex %in% "Male"], na.rm = T)/ mean(dataset$Intensity[dataset$genotype %in% "Control" & dataset$sex %in% "Male"], na.rm = T)
      }else{
        # if statistical analysis can be performed, test normality before data tranformation
        Normality_before[[current_gene]][i] <- shapiro.test(data_raw$Intensity_raw)$p.value
        # test normality after data transformation for male and female mice
        Normality_F[[current_gene]][i] <- shapiro.test(data_raw$Intensity[data_raw$sex %in% "Female"])$p.value
        Normality_M[[current_gene]][i] <- shapiro.test(data_raw$Intensity[data_raw$sex %in% "Male"])$p.value
        #count number of WT mice and KO mice
        numb_null[[current_gene]][i] <- sum(dataset$genotype %in% "Control")
        numb_gene[[current_gene]][i] <- sum(dataset$genotype %in% "Gene")
        # statistical test for each model
        model_formula_null <- null_model_genotype(dataset, "Intensity")
        model_formula_genotype <- final_genotype_model(dataset, "Intensity")
        #test genotyp effect
        geno_effect <- testing_genotype_effect(dataset, "Intensity")
        #test sex effect
        formula_sex_null = null_model_sex(dataset,"Intensity")
        formula_sex_full = final_sex_model(dataset,"Intensity")
        sex_effect= testing_sex_effect(dataset,"Intensity")
        #test genotype-sex interaction effect
        formula_interaction_null <- null_model_Interaction(dataset, "Intensity")
        Interct_effect <- testing_Interaction_effect(dataset, "Intensity")
        result <- gls(Intensity ~ genotype + sex + genotype: sex, data = dataset, na.action = 'na.exclude')
        summary <- nlme:::summary.gls(result)$tTable 
        #p-values
        pvalue_Stage1[[current_gene]][i] <- geno_effect
        pvalue_Stage1[[current_gene]][i] <- sex_effect
        pvalue_Stage2[[current_gene]][i] <- Interct_effect
    
        # Individual test for classying the SD effect without grouping by sex
        result_all<- tryCatch({
          gls(Intensity ~ genotype, data = dataset, na.action = 'na.exclude')
        }, error = function(er){
          NA
        }) 
        #p-values without grouping by sex
        if(length(result_all) <= 1){
          sep_allKO_pval[[current_gene]][i] <- NA
          sep_allKO_estimate[[current_gene]][i] <- NA
        }else{
          sep_all <- nlme:::summary.gls(result_all)$tTable 
          sep_allKO_pval[[current_gene]][i] <- sep_all[2, 4]
          sep_allKO_estimate[[current_gene]][i] <- sep_all[2, 1]
        }
        
        # subset data for female mice an male mice
        sub_F <- dataset[dataset$sex %in% "Female", ]
        sub_M <- dataset[dataset$sex %in% "Male", ]
        # individual comparison for female KO to WT mice
        result_F <- tryCatch({
          gls(Intensity ~ genotype, data = sub_F, na.action = 'na.exclude')
        }, error = function(er){
          NA
        })
        if(length(result_F) <= 1){
          sep_FvKO_pval[[current_gene]][i] <- NA
          sep_FvKO_estimate[[current_gene]][i] <- NA
        }else{
          sep_FvKO <- nlme:::summary.gls(result_F)$tTable 
          sep_FvKO_pval[[current_gene]][i] <- sep_FvKO[2, 4]
          sep_FvKO_estimate[[current_gene]][i] <- sep_FvKO[2, 1]
        }
        # individual comparison for male KO to WT mice
        result_M <- tryCatch({
          gls(Intensity ~ genotype, data = sub_M, na.action = 'na.exclude')
        }, error = function(er){
          NA
        }) 
        if(length(result_M) <= 1){
          sep_MvKO_pval[[current_gene]][i] <- NA
          sep_MvKO_estimate[[current_gene]][i] <- NA
        }else{
          sep_MvKO <- nlme:::summary.gls(result_M)$tTable 
          sep_MvKO_pval[[current_gene]][i] <- sep_MvKO[2, 4]
          sep_MvKO_estimate[[current_gene]][i] <- sep_MvKO[2, 1]
        }
        
        # fold change calculation
        foldchange_FvKO[[current_gene]][i] <- mean(dataset$Intensity_raw[dataset$genotype %in% "Gene" & dataset$sex %in% "Female"], na.rm = T)/ mean(dataset$Intensity_raw[dataset$genotype %in% "Control" & dataset$sex %in% "Female"], na.rm = T)
        foldchange_MvKO[[current_gene]][i] <- mean(dataset$Intensity_raw[dataset$genotype %in% "Gene" & dataset$sex %in% "Male"], na.rm = T)/ mean(dataset$Intensity_raw[dataset$genotype %in% "Control" & dataset$sex %in% "Male"], na.rm = T)
    
        # subset data for wildtype mice by sex, assess sex effect in WT group    
        ctrl_data <- dataset[dataset$genotype%in%"Control"]
        ctrl_test <- tryCatch({
          gls(Intensity~sex, data = ctrl_data, na.action = 'na.exclude')
        }, error = function(er){
          NA
        }) 
        if(length(ctrl_test) <= 1){
          KO_sex_pval[[current_gene]][i] <- NA
          KO_sex_estimate[[current_gene]][i] <- NA
        }else{
          ctrl_result <- nlme:::summary.gls(ctrl_test)$tTable 
          ctrl_sex_pval[[current_gene]][i] <- ctrl_result[2, 4]
          ctrl_sex_estimate[[current_gene]][i] <- ctrl_result[2, 1]
        }
        
        # subset data for KO strain mice by sex, assess sex effect in KO group 
        KO_data = dataset[dataset$genotype %in% "Gene"]
        KO_test<- tryCatch({
          gls(Intensity~sex, data = KO_data, na.action = 'na.exclude')
        }, error = function(er){
          NA
        }) 
        if(length(KO_test) <= 1){
          KO_sex_pval[[current_gene]][i] <- NA
          KO_sex_estimate[[current_gene]][i] <- NA
        }else{
          KO_result <- nlme:::summary.gls(KO_test)$tTable 
          KO_sex_pval[[current_gene]][i] <- KO_result[2, 4]
          KO_sex_estimate[[current_gene]][i] <- KO_result[2, 1]
        }
      }
    }
  }
}
# calculate significance ratio
sapply(pvalue_Stage1,function(x){sum(x<0.05, na.rm = TRUE)})  

#Output result
setwd("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX") 
for(i in 1 : length(names(pvalue_Stage1))){
  fwrite(data.table(label = All_Metabolomics_f$label, assay = All_Metabolomics_f$Assay, Metabolite=All_Metabolomics_f$CompoundName, Normality_before = Normality_before[[i]], Normality_F = Normality_F[[i]],  Normality_M = Normality_M[[i]] , numb_null = numb_null[[i]], numb_gene = numb_gene[[i]], pvalue_Stage1 = pvalue_Stage1[[i]], pvalue_stage1b = pvalue_Stage1b[[i]], pvalue_stage2 = pvalue_Stage2[[i]], fdr_Stage1=fdr_Stage1[[i]], fdr_Stage1b=fdr_Stage1b[[i]], fdr_Stage2=fdr_Stage2[[i]], sep_allKO_pval = sep_allKO_pval[[i]], sep_allKO_estimate = sep_allKO_estimate[[i]], sep_FvKO_pval = sep_FvKO_pval[[i]], sep_FvKO_estimate = sep_FvKO_estimate[[i]], sep_MvKO_pval = sep_MvKO_pval[[i]], sep_MvKO_estimate = sep_MvKO_estimate[[i]], sep_FvKO_fdr=p.adjust(sep_FvKO_pval[[i]],'fdr'), sep_MvKO_fdr=p.adjust(sep_MvKO_pval[[i]],'fdr'), foldchange_FvKO = foldchange_FvKO[[i]], foldchange_MvKO = foldchange_MvKO[[i]], ctrl_sex_pval = ctrl_sex_pval[[i]], fdr_ctrl_sex = fdr_ctrl_sex[[i]], ctrl_sex_estimate = ctrl_sex_estimate[[i]], KO_sex_pval= KO_sex_pval[[i]], fdr_KO_sex = fdr_KO_sex[[i]], KO_sex_estimate = KO_sex_estimate[[i]]), paste0("5 TwoStage METAbol_KOSEX,(",names(pvalue_Stage1)[i],").csv"))
}


#5 classification for pie chart and bar plot: Metabolite variable / genptype: sex interaction effect and direction
setwd("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
# read filelist
files_metmerge<-list.files(pattern="^5 TwoStage PHENOtype (.*),genoSD",full.names = T)
#load genelist
GeneList<-read.csv("GeneList.csv",header = TRUE)
GeneList<-(GeneList$Genotype)
Genelist<-as.character(GeneList)

PieChart_Diff_Dirct<-data.frame(matrix(ncol = 14, nrow = 0))
coltitle<-c('Genotype','Different direction',"Different direction2",'Different size',"Different size2","One sex_Female only", "One sex_Male only","One sex_Female only2","One sex_Male only2", "cannot classify", "cannot classify2","Genotype effect with no sex effect(real)",'Not Significant',"Metabolites count")
colnames(PieChart_Diff_Dirct)<-coltitle
#genotype effect classification for each genotype
for (i in 1:length(files_metmerge)){
  # read file
  met<-read.csv(files_metmerge[i])  
  g<-i
  print(i)
  met1<-met[!is.na(met$pvalue_Stage1),]
  met1$est_dirct<-"Not significant"
  t<-length(met1$pvalue_Stage1)
  met1$est_dirct[(met1$pvalue_Stage1 < 0.05&met1$pvalue_stage2 < 0.05)&(met1$sep_FvKO_pval < 0.05&met1$sep_MvKO_pval < 0.05 )& ((met1$sep_FvKO_estimate>0& met1$sep_MvKO_estimate<0) | (met1$sep_FvKO_estimate<0 & met1$sep_MvKO_estimate>0))] = "Different direction"
  met1$est_dirct[(met1$pvalue_Stage1 < 0.05&met1$pvalue_stage2 < 0.05)&(met1$sep_FvKO_pval < 0.05 & met1$sep_MvKO_pval < 0.05) & ((met1$sep_FvKO_estimate>0 & met1$sep_MvKO_estimate>0)  | (met1$sep_FvKO_estimate<0 & met1$sep_MvKO_estimate<0))] = "Different size"
  met1$est_dirct[(met1$pvalue_Stage1 < 0.05 & met1$pvalue_stage2 < 0.05)&(met1$sep_FvKO_pval < 0.05&met1$sep_MvKO_pval >= 0.05)] = "One sex_Female only"
  met1$est_dirct[(met1$pvalue_Stage1 < 0.05 & met1$pvalue_stage2 < 0.05)&(met1$sep_FvKO_pval >= 0.05&met1$sep_MvKO_pval < 0.05)] = "One sex_Male only"
  met1$est_dirct[(met1$pvalue_Stage1 >= 0.05 & met1$pvalue_stage2 < 0.05)&(met1$sep_FvKO_pval < 0.05&met1$sep_MvKO_pval < 0.05 )& ((met1$sep_FvKO_estimate>0& met1$sep_MvKO_estimate<0) | (met1$sep_FvKO_estimate<0 & met1$sep_MvKO_estimate>0))] = "Different direction2"
  met1$est_dirct[(met1$pvalue_Stage1 >= 0.05 & met1$pvalue_stage2 < 0.05)&(met1$sep_FvKO_pval < 0.05 & met1$sep_MvKO_pval < 0.05) & ((met1$sep_FvKO_estimate>0 & met1$sep_MvKO_estimate>0)  | (met1$sep_FvKO_estimate<0 & met1$sep_MvKO_estimate<0))] = "Different size2"
  met1$est_dirct[(met1$pvalue_Stage1 >= 0.05 & met1$pvalue_stage2 < 0.05)&(met1$sep_FvKO_pval < 0.05&met1$sep_MvKO_pval >= 0.05)] = "One sex_Female only2"
  met1$est_dirct[(met1$pvalue_Stage1 >= 0.05 & met1$pvalue_stage2 < 0.05)&(met1$sep_FvKO_pval >= 0.05&met1$sep_MvKO_pval < 0.05)] = "One sex_Male only2"
  met1$est_dirct[(met1$pvalue_Stage1 < 0.05 & met1$pvalue_stage2 < 0.05)&(met1$sep_FvKO_pval >= 0.05&met1$sep_MvKO_pval >= 0.05 )] = "cannot classify"
  met1$est_dirct[(met1$pvalue_Stage1 >= 0.05 & met1$pvalue_stage2 < 0.05)&(met1$sep_FvKO_pval >= 0.05&met1$sep_MvKO_pval >= 0.05 )] = "cannot classify2"
  met1$est_dirct[(met1$pvalue_Stage1 < 0.05 & met1$pvalue_stage2 >= 0.05)] = "Genotype effect with no sex effect(real)"                 
  a<-sum(met1$est_dirct == "Different direction")
  a1<-sum(met1$est_dirct == "Different direction2")
  b <- sum(met1$est_dirct == "Different size")
  b1 <- sum(met1$est_dirct == "Different size2")
  c <- sum(met1$est_dirct == "One sex_Female only")
  c1 <- sum(met1$est_dirct == "One sex_Male only")
  c2 <- sum(met1$est_dirct == "One sex_Female only2")
  c3 <- sum(met1$est_dirct == "One sex_Male only2")
  d<-sum(met1$est_dirct == "cannot classify")
  d1<-sum(met1$est_dirct == "cannot classify2")
  e<-sum(met1$est_dirct == "Genotype effect with no sex effect(real)")
  n<-sum(met1$est_dirct == "Not significant")
  PieChart_Diff_Dirct<-rbind(PieChart_Diff_Dirct,c(g,a,a1,b,b1,c,c1,c2,c3,d,d1,e,n,t))
  colnames(PieChart_Diff_Dirct)<-coltitle 
}
write.csv(PieChart_Diff_Dirct,file = "5 Metabol_classification.csv")
