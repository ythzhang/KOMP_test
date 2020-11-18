#######################################################################################################
# Statistics methods are adapted from : Supporting Data of Prevalence of Sexual Dimorphism in Mammalian Phenotypic Traits (2017).
library(nlme)
library(tidyverse)
library(wcmc)
#######################################################################################################
# Load raw data file for Phenotype results
setwd("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
Raw_Phenotype_UCDavis = wcmc::read_data("Supplementary Data 2.csv")
Raw_Phenotype_UCDavis_p = Raw_Phenotype_UCDavis$p
Raw_Phenotype_UCDavis_f = Raw_Phenotype_UCDavis$f
Raw_Phenotype_UCDavis_e = Raw_Phenotype_UCDavis$e_cat_matrix

# check data structure and convert to numeric data frame 
Raw_Phenotype_UCDavis_e = apply(Raw_Phenotype_UCDavis_e,2,as.numeric)
class(Raw_Phenotype_UCDavis_e[1,1])
# ============================================================================================================================
# replace missing values by minimum values
# If some parameter is missing that's because it either couldn't be obtained or had to be removed during QC, thus the mouse with missing values are assumed to be normal and missing values were replace with minimum values
Phenotype_UCDavis_e_no_mising = Raw_Phenotype_UCDavis_e
unique_genes = unique(Raw_Phenotype_UCDavis_p$Genotype)
continuous_index = c()
index_cat = c()
n = 0
for(g in 1:length(unique_genes)){
    # read gene
    current_gene = unique_genes[g]
    # subset data for each genotype based on the gene symbol
    current_label = Raw_Phenotype_UCDavis_p$label[Raw_Phenotype_UCDavis_p$Genotype %in% c(current_gene)]
    current_e = Raw_Phenotype_UCDavis_e[, current_label]
    # subset data for each genotype female mice based on the gene symbol
    current_label_F = Raw_Phenotype_UCDavis_p$label[Raw_Phenotype_UCDavis_p$Genotype %in% c(current_gene) & Raw_Phenotype_UCDavis_p$Gender %in% "Female"]
    current_e_F = Raw_Phenotype_UCDavis_e[, current_label_F]
    # subset data for each genotype male mice based on the gene symbol
    current_label_M = Raw_Phenotype_UCDavis_p$label[Raw_Phenotype_UCDavis_p$Genotype %in% c(current_gene) & Raw_Phenotype_UCDavis_p$Gender %in% "Male"]
    current_e_M = Raw_Phenotype_UCDavis_e[,current_label_M]
    # replace missing values in WT mice
  if (current_gene == "null"){
  for(i in 1 : nrow(current_e)){
    # confirm if parameter is continuous traits or not
    continuous = sum(is.na(as.numeric(Raw_Phenotype_UCDavis_e[i,!is.na(Raw_Phenotype_UCDavis_e[i,])]))) == 0
    if(continuous){
      continuous_index[i] = TRUE
      if((sum(is.na(current_e_F[i,])) >= 14) & (sum(is.na(current_e_M[i,])) >= 14)){
        current_e_F[i,] = NA
        current_e_M[i,] = NA
      }else if((sum(is.na(current_e_F[i,])) >= 14) & (sum(is.na(current_e_M[i,]))< 14)){
        current_e_F[i,is.na(current_e_F[i,])] = runif(sum(is.na(current_e_F[i,])),min = 0.95, max = 1.05) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
        current_e_M[i,] = current_e_M[i,]
      }else if((sum(is.na(current_e_F[i,])) < 14) & (sum(is.na(current_e_M[i,])) >= 14)){
        current_e_F[i,] = current_e_F[i,]
        current_e_M[i,is.na(current_e_M[i,])] = runif(sum(is.na(current_e_M[i,])),min = 0.95,max = 1.05) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
      }else{
        current_e_F[i,] = current_e_F[i,]
        current_e_M[i,] = current_e_M[i,]
      }
    }else{
      continuous_index[i] = FALSE
      n <- n + length(index_cat)
      index_cat[n] <- i
    }
    }}else{
      # replace missing values in genotype mice
      for(i in 1 : nrow(current_e)){
        continuous = sum(is.na(as.numeric(Raw_Phenotype_UCDavis_e[i,!is.na(Raw_Phenotype_UCDavis_e[i,])]))) == 0
        # when parameter is continuous trait    
        if(continuous){
          continuous_index[i] = TRUE
        if(sum(is.na(current_e[i,])) >= 4){
          current_e_F[i,] = NA
          current_e_M[i,] = NA
        }else if(sum(is.na(current_e_F[i,])) == 3 & sum(is.na(current_e_M[i,])) == 0){
          current_e_F[i,is.na(current_e_F[i,])] = runif(sum(is.na(current_e_F[i,])),min = 0.95, max = 1.05) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
          current_e_M[i,] = current_e_M[i,]
        }else if(sum(is.na(current_e_F[i,])) == 0 & sum(is.na(current_e_M[i,])) == 3){
          current_e_F[i,] = current_e_F[i,]
          current_e_M[i,is.na(current_e_M[i,])] = runif(sum(is.na(current_e_M[i,])),min = 0.95, max = 1.05) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
        }else if(sum(is.na(current_e_F[i,])) == 2 & sum(is.na(current_e_M[i,])) == 1){
          current_e_F[i,is.na(current_e_F[i,])] = runif(sum(is.na(current_e_F[i,])), min = 0.95, max = 1.05) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
          current_e_M[i,is.na(current_e_M[i,])] = runif(sum(is.na(current_e_M[i,])),min = 0.95, max = 1.05) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
        }else if(sum(is.na(current_e_F[i,]))==1 & sum(is.na(current_e_M[i,])) == 2){
          current_e_F[i,is.na(current_e_F[i,])] = runif(sum(is.na(current_e_F[i,])),min = 0.95, max = 1.05) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
          current_e_M[i,is.na(current_e_M[i,])] = runif(sum(is.na(current_e_M[i,])),min = 0.95, max = 1.05) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
        }else if(sum(is.na(current_e[i,])) == 2){
          if(sum(is.na(current_e_F[i,])) == 2 & sum(is.na(current_e_M[i,])) == 0){
            current_e_F[i,is.na(current_e_F[i,])] = runif(sum(is.na(current_e_F[i,])),min = 0.95, max = 1.05) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
            current_e_M[i,] = current_e_M[i,]
          }else if(sum(is.na(current_e_F[i,])) == 0 & sum(is.na(current_e_M[i,])) == 2){
            current_e_F[i,] = current_e_F[i,]
            current_e_M[i,is.na(current_e_M[i,])] = runif(sum(is.na(current_e_M[i,])),min = 0.95, max = 1.05) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
          }else if(sum(is.na(current_e_F[i,])) == 1 & sum(is.na(current_e_M[i,])) == 1){
            current_e_F[i,is.na(current_e_F[i,])] = runif(sum(is.na(current_e_F[i,])),min = 0.95, max = 1.05) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
            current_e_M[i,is.na(current_e_M[i,])] = runif(sum(is.na(current_e_M[i,])),min = 0.95, max = 1.05) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
          }
        }else if(sum(is.na(current_e[i,])) <= 1){
          current_e_F[i,is.na(current_e_F[i,])] = runif(sum(is.na(current_e_F[i,])), min = 0.95, max = 1.05) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
          current_e_M[i,is.na(current_e_M[i,])] = runif(sum(is.na(current_e_M[i,])),min = 0.95, max = 1.05) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
        }
      }else{
        #when parameter is not continuous trait
        continuous_index[i] = FALSE
        n <- n + length(index_cat)
        index_cat[n] <- i
      }
    }
    }
  Phenotype_UCDavis_e_no_mising[,current_label_F] = current_e_F
  Phenotype_UCDavis_e_no_mising[,current_label_M] = current_e_M
}

# construct new data file with replaced values
Phenotype_e_no_mising = Phenotype_UCDavis_e_no_mising[continuous_index,]
nrow(Phenotype_e_no_mising)
Phenotype_e_no_mising_new <- Phenotype_e_no_mising
# class(Phenotype_e_no_mising_new[1,1])
str(Phenotype_e_no_mising_new)
Phenotype_f = Raw_Phenotype_UCDavis_f
Phenotype_p = Raw_Phenotype_UCDavis_p
# subset data for WT mice
null_label = Phenotype_p$label[Phenotype_p$Genotype %in% c("null")]
null_e = Phenotype_e_no_mising_new[, null_label]
null_p = Phenotype_p[Phenotype_p$Genotype %in% "null",]
# subset data for WT female mice
null_label = Phenotype_p$label[Phenotype_p$Genotype %in% c("null") & Phenotype_p$Gender %in% "Female"]
null_e_female = Phenotype_e_no_mising_new[, null_label]
# subset data for WT male mice
null_label = Phenotype_p$label[Phenotype_p$Genotype %in% c("null") & Phenotype_p$Gender %in% "Male"]
null_e_male = Phenotype_e_no_mising_new[, null_label]
# =========================================================================================================================================
# =========================================================================================================================================
# 5 for each genotype, check which metabolite are associated with genotype-gender, gene, gender, or not significant at all.
# --------------------------------------------------------------------------------------------------------
#load package nlme
library(nlme)
# library(car)
# Function 1 function to build model for testing of fixed effects
# Assumption: No batch effect and weight is not a independet variance 
# Goal: 
model_forFIXEDtest <- function(dataset, PHenoVAlue){
  model.formula <- as.formula(paste(PHenoVAlue, "~", paste("genotype", "sex", "genotype:sex", sep= "+")))
  model=gls(model.formula, dataset, na.action = "na.omit")
}
# -------------------------------------------------------------------------------------------------------------------------------
#Function 2: testing the fixed effects and building genotype model formula
# Goal: Genotype is always include in this model. The goal is to test genotype and sex effect and based on the output build the final genotype model formula for later testing.  
# anova function: null hypothesis, the regression coefficients are equal to zero  
#                 alternative hypothesis that the regression coefficient are not equal to zero.
# the null hypothesis was rejected when p-values < 0.05 (e.g.,accept the alternative hypothesis, the components of the model should be included in later analysis)
# Note a complexity surrounds the interaction term  - if it is significant but gender is excluded it is excluded. 
final_genotype_model <- function(dataset, PHenoVAlue){
  model_afterFIXED=model_forFIXEDtest(dataset, PHenoVAlue)
  anova_results = anova(model_afterFIXED, type="marginal")$"p-value" < 0.05
  keepSex = anova_results[3]
  keepInteraction = anova_results[4]
  if(keepSex && keepInteraction){
    return(model.formula <- as.formula(paste(PHenoVAlue, "~", paste("genotype", "sex", "genotype:sex", sep= "+"))))
  }else if(keepSex &&  !keepInteraction){
    return(model.formula <- as.formula(paste(PHenoVAlue, "~", paste("genotype", "sex", sep= "+"))))
  }else if(!keepSex &&  !keepInteraction){
    return(model.formula <- as.formula(paste(PHenoVAlue, "~", "genotype")))
  }else if(!keepSex  && keepInteraction){
    return(model.formula <- as.formula(paste(PHenoVAlue, "~", paste("genotype", "sex", "genotype:sex", sep= "+"))))
  }
}
#---------------------------------------------------------------------------------------------------------------------------
#Function 3: testing the fixed effects and building final Genotype effect null model
# Goal:  to test fixed effects of the model and based on the output build the final null model formula for later testing - as a null model it automatically excludes genotype and interaction term. 
# The anova function tests the fixed effects associated by treatment with a null hypothesis that the regression coefficients are equal to zero  and an alternative hypothesis that the regression coefficient are not equal to zero.
# If the p-values of these tests are less than 0.05 we reject the null and accept the alternative that the are significant components of the model and should be included. 
# If no terms are significant a model can be build with just an intercept element this is specified as  "model.formula <- as.formula(paste(depVariable, "~", "1"))"
null_model_genotype <- function(dataset, PHenoVAlue){ 
  model_afterFIXED = model_forFIXEDtest(dataset, PHenoVAlue)
  anova_results = anova(model_afterFIXED, type="marginal")$"p-value" < 0.05
  # anova_results = anova.lm(model_afterFIXED, type="marginal")[p-value] < 0.05
  keepSex = anova_results[3]
  keepInteraction = anova_results[4]     
  if(!keepSex && !keepInteraction){
    return(model.formula <- as.formula(paste(PHenoVAlue, "~", "1")))
  }else{
    return(model.formula <- as.formula(paste(PHenoVAlue, "~", "sex")))
  } 
}
#----------------------------------------------------------------------------------
#Function 3-b: testing the fixed effects and building final Interaction effect null model
null_model_Interaction <- function(dataset, PHenoVAlue){ 
  model_afterFIXED = model_forFIXEDtest(dataset, PHenoVAlue)
  anova_results = anova(model_afterFIXED, type="marginal")$"p-value" < 0.05
  keepSex = anova_results[3]
  keepGenotype = anova_results[2]     
  if(!keepSex && !keepGenotype){
    return(model.formula <- as.formula(paste(PHenoVAlue, "~", "1")))
  }else{
    return(model.formula <- as.formula(paste(PHenoVAlue, "~", paste("genotype","sex",sep = "+"))))
  } 
}
# ------------------------------------------------------------------------------------------------
#Function 4: testing the genotype effect
# The goal of the function is to give a pvalue of whether genotype effect is significant by compare the genotype model with the null model with the anova function
# The model_formula_null and model_formula_genotype are called to define the models for comparison.
# For each possible combination,  then an anova model is used to report the pvalue.
# For testing the genotype effect we use method= ML
testing_genotype_effect <- function(dataset,PHenoVAlue){
  #call functions to determine model.formula
  model_formula_null = null_model_genotype(dataset,PHenoVAlue)
  model_formula_genotype = final_genotype_model(dataset,PHenoVAlue)
  model_genotype = gls(model_formula_genotype,dataset,method = 'ML', na.action = "na.omit")
  model_null= gls(model_formula_null, dataset, method = "ML", na.action = "na.omit")
  return(pvalue_genotype=(anova(model_genotype,model_null)$`p-value`[2]))
}
# ------------------------------------------------------------------------------------------------
#Function 4-b: testing the interaction effect
# The goal of the function is to give a pvalue of whether genotype effect is significant by compare the genotype model with the null model with the anova function
# The model_formula_null and model_formula_genotype are called to define the models for comparison.
# For each possible combination,  then an anova model is used to report the pvalue.
# For testing the genotype effect we use method= ML
testing_Interaction_effect<-function(dataset, PHenoVAlue){
  #call functions to determine model.formula
  formula_interaction_null = null_model_Interaction(dataset,PHenoVAlue)
  model_interaction_full = gls(PHenoVAlue~genotype+sex+genotype:sex,dataset, method = 'ML',na.action = "na.omit")
  model_null_interaction = gls(formula_interaction_null, dataset,method = 'ML', na.action = "na.omit")
  return(pvalue_Stage2 = (anova(model_interaction_full, model_null_interaction)$`p-value`[2]))
}

#Function 4-d: testing the sex effect
null_model_sex <- function(dataset, PHenoVAlue){ 
  model_afterFIXED = model_forFIXEDtest(dataset, PHenoVAlue)
  anova_results = anova(model_afterFIXED, type = "marginal")$"p-value" < 0.05
  keepInteraction = anova_results[4]
  keepGenotype = anova_results[2]     
  if(!keepInteraction && !keepGenotype){
    return(model.formula <- as.formula(paste(PHenoVAlue, "~", "1")))
  }else{
    return(model.formula <- as.formula(paste(PHenoVAlue, "~", "genotype")))
  } 
}
# ------------------------------------------------------------------------------------------------
final_sex_model <- function(dataset, PHenoVAlue){
  model_afterFIXED = model_forFIXEDtest(dataset, PHenoVAlue)
  anova_results = anova(model_afterFIXED, type = "marginal")$"p-value" < 0.05
  keepGenotype = anova_results[2]
  keepInteraction = anova_results[4]
  if( keepGenotype && keepInteraction){
    return(model.formula <- as.formula(paste(PHenoVAlue, "~", paste("genotype", "sex", "genotype:sex", sep= "+"))))
  }else if(keepGenotype &&  !keepInteraction){
    return(model.formula <- as.formula(paste(PHenoVAlue, "~", paste("genotype", "sex", sep= "+"))))
  }else if(!keepGenotype &&  !keepInteraction){
    return(model.formula <- as.formula(paste(PHenoVAlue, "~", "sex" )))
  }else if(!keepGenotype  && keepInteraction){
    return(model.formula <- as.formula(paste(PHenoVAlue, "~", paste("genotype", "sex", "genotype:sex", sep= "+"))))
  }
}
testing_sex_effect <- function(dataset, PHenoVAlue){
  #call functions to determine model.formula
  formula_sex_null = null_model_sex(dataset, PHenoVAlue)
  formula_sex_full = final_sex_model(dataset,PHenoVAlue)
  model_null_sex = gls(formula_sex_null, dataset,method = "ML", na.action = "na.omit")
  model_full_sex = gls(formula_sex_full, dataset, method = "ML", na.action = "na.omit")
  return(pvalue_Stage1.5 = (anova(model_full_sex, model_null_sex)$`p-value`[2]))
}
#-----------------------------------------------------------------------------------------
#Function 5: following function returns the final model which is needed for the diagnostic plots 
#Goal of this function is to return a model output that is the final model following all the previous tests which can be used to generate diagnostics plots and output the final model details. 
# The model_formula_null and model_formula_genotype are called to define the models for comparison. 
# The keep_batch and keep_equal variance function are called this allows if but rules to build the correct model comparison ie if batch is included then we use a mixed model etc
# For each possible combination,  then the model is fitted to the data and the model reported as the output.
#this is a separate function to above as for the estimates of the fixed effects we use method=REML whilst for the genotype test we used method= ML   
#na.action = na.exclude as in this form the residue calculated from the model output have the same length as the original datafile
finalmodel <- function(dataset, PHenoVAlue){
  #call functions to determine model.formula
  model_formula_null = null_model_genotype(dataset,PHenoVAlue)
  model_genotype = final_genotype_model(dataset, PHenoVAlue)
  model_genotype = gls(model_genotype, dataset, na.action = "na.exclude")
}
# ---------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------
# Starting TwoStage significance assessment using Two-way ANOVA
Normality_before = Normality_after =Normality_F = Normality_M= list()
numb_null = numb_gene = list()
pvalue_Stage1 = pvalue_Stage1.5 = pvalue_Stage2 = list()
fdr_Stage1 = fdr_Stage1.5 = fdr_Stage2 = list()
ctrl_sex_pval = ctrl_sex_estimate = fdr_ctrl_sex = list() 
KO_sex_pval = KO_sex_estimate = fdr_KO_sex = list()
sep_allKO_pval = sep_allKO_estimate =  list()
sep_FvKO_pval = sep_FvKO_estimate = sep_MvKO_pval = sep_MvKO_estimate = list()
foldchange_all = foldchange_FvKO = foldchange_MvKO = list()
#statistical analysis for each genotype
for(g in 2:length(unique_genes)){
  #read gene
  current_gene = unique_genes[g]
  print(g)
  # subset data for each genotype
  current_label = Phenotype_p$label[Phenotype_p$Genotype %in% c(current_gene)]
  current_e = Phenotype_e_no_mising_new[, current_label]
  current_e = as.data.frame(current_e)
  current_p = Phenotype_p[Phenotype_p$label %in% current_label,]
  numb_null[[current_gene]] = numb_gene[[current_gene]] = c()
  Normality_before[[current_gene]] = Normality_after[[current_gene]] = Normality_F[[current_gene]] = Normality_M[[current_gene]] = c()
  pvalue_Stage1[[current_gene]] = pvalue_Stage1.5[[current_gene]] = pvalue_Stage2[[current_gene]] = c()
  fdr_Stage1[[current_gene]] = fdr_Stage1.5[[current_gene]] = fdr_Stage2[[current_gene]] = c()
  ctrl_sex_pval[[current_gene]] = ctrl_sex_estimate[[current_gene]] = fdr_ctrl_sex[[current_gene]] = c() 
  KO_sex_pval[[current_gene]] = KO_sex_estimate[[current_gene]] = fdr_KO_sex[[current_gene]] = c()
  sep_allKO_pval[[current_gene]] = sep_allKO_estimate[[current_gene]] = fdr_allKO = c()
  sep_FvKO_pval[[current_gene]] = sep_FvKO_estimate[[current_gene]] = sep_MvKO_pval[[current_gene]] = sep_MvKO_estimate[[current_gene]] = c()
  foldchange_all[[current_gene]] = c()
  foldchange_FvKO[[current_gene]] = c()
  foldchange_MvKO[[current_gene]] = c()
  
  # statistical test for each phenotype parameter for each genotyope
  for(i in 1:nrow(Phenotype_e_no_mising)){
    # if missing values > 70%, then not significant
    if((sum(is.na(current_e[i,])) >= (ncol(current_e)-2) )|(sum(is.na(null_e[i,])) >= (ncol(null_e)-6) )){
      numb_null[[current_gene]][i]=numb_gene[[current_gene]][i] = NA
      Normality_before[[current_gene]][i] = NA
      Normality_after[[current_gene]][i] = NA
      Normality_F[[current_gene]][i] = NA
      Normality_M[[current_gene]][i] = NA
      #p-values NA
      pvalue_Stage1[[current_gene]][i] = NA
      pvalue_Stage1.5[[current_gene]][i] = NA
      pvalue_Stage2[[current_gene]][i] = NA
      ctrl_sex_pval[[current_gene]][i] = ctrl_sex_estimate[[current_gene]][i] = fdr_ctrl_sex[[current_gene]][i] = NA
      KO_sex_pval[[current_gene]][i] = KO_sex_estimate[[current_gene]][i] = NA
      sep_allKO_pval[[current_gene]][i] = sep_allKO_estimate[[current_gene]][i] = NA
      sep_FvKO_pval[[current_gene]][i] = sep_FvKO_estimate[[current_gene]][i] = NA
      sep_MvKO_pval[[current_gene]][i] = sep_MvKO_estimate[[current_gene]][i] = NA
      #fold-change NA
      foldchange_all[[current_gene]][i] = NA
      foldchange_FvKO[[current_gene]][i] = NA
      foldchange_MvKO[[current_gene]][i] = NA
    }else{
      # if missing values < 70% do statistical analysis
      data_raw <- data.table(y =(c(current_e[i,],null_e[i,])), group = rep(c("Gene","Control"), c(ncol(current_e), ncol(null_e))), sex=c(current_p$Gender,null_p$Gender))
      colnames(data_raw) = c("PHenoVAlue_raw","genotype","sex")
      data_raw$PHenoVAlue_raw <-unlist(data_raw$PHenoVAlue_raw)
      # set as numeric 
      data_raw$PHenoVAlue_raw <- as.numeric(as.character(data_raw$PHenoVAlue_raw))
      # remove sissing values and transform data
      data_raw <- data_raw[!is.na(data_raw$PHenoVAlue_raw)]
      data_raw$PHenoVAlue <- rankNorm(data_raw$PHenoVAlue_raw)
      # check outliers
      outliers_nullF_raw <- boxplot.stats(data_raw$PHenoVAlue[data_raw$sex %in% "Female" & data_raw$genotype %in% "Control"])$out
      # remove outliers if necessary  
      if(length(outliers_nullF_raw) == 0){
        dataset2_raw <- data_raw
      }else{
        dataset2_raw <- data_raw[-which((data_raw$PHenoVAlue %in% outliers_nullF_raw) & (data_raw$sex %in% "Female") & (data_raw$genotype %in% "Control")),]
      }
      outliers_nullM_raw <- boxplot.stats(dataset2_raw$PHenoVAlue[dataset2_raw$sex %in% "Male" & dataset2_raw$genotype %in% "Control"])$out
      if(length(outliers_nullM_raw) == 0){
        dataset3_raw <- dataset2_raw
      }else{
        dataset3_raw<-dataset2_raw[-which((dataset2_raw$PHenoVAlue %in% outliers_nullM_raw) & (dataset2_raw$sex %in% "Male") & (dataset2_raw$genotype %in% "Control")),]
      }
      outliers_genoF_raw <- boxplot.stats(dataset3_raw$PHenoVAlue[dataset3_raw$sex %in% "Female" & dataset3_raw$genotype %in% "Gene"])$out
      if(length(outliers_genoF_raw) == 0){
        dataset4_raw <- dataset3_raw
      }else{
        dataset4_raw <- dataset3_raw[-which((dataset3_raw$PHenoVAlue %in% outliers_genoF_raw) & (dataset3_raw$sex %in% "Female") & (dataset3_raw$genotype %in% "Gene")),]
      }
      outliers_genoM_raw <- boxplot.stats(dataset4_raw$PHenoVAlue[dataset4_raw$sex %in% "Male" & dataset4_raw$genotype %in% "Control"])$out
      if(length(outliers_genoM_raw) == 0){
        dataset_raw <- dataset4_raw
      }else{
        dataset_raw <-dataset4_raw[-which((dataset4_raw$PHenoVAlue %in% outliers_genoM_raw) & (dataset4_raw$sex %in% "Male") & (dataset4_raw$genotype %in% "Control")),]
      }
      # remove missing values
      dataset <- dataset_raw[!is.na(dataset_raw$PHenoVAlue),]
      # test statistcal model 
      model_afterFIXED=tryCatch({
        model_forFIXEDtest(dataset, "PHenoVAlue")
      }, error = function(er){
        NA
      })
      if(length(model_afterFIXED)<=1){
        #count number for WT mice and KO mice
        numb_null[[current_gene]][i] = sum(dataset$genotype %in% "Control")
        numb_gene[[current_gene]][i] = sum(dataset$genotype %in% "Gene")
        Normality_before[[current_gene]][i] = NA
        Normality_after[[current_gene]][i] = NA
        Normality_F[[current_gene]][i] = NA
        Normality_M[[current_gene]][i] =NA
        #p-values not availabel, NA
        pvalue_Stage1[[current_gene]][i] = NA
        pvalue_Stage1.5[[current_gene]][i] = NA
        pvalue_Stage2[[current_gene]][i] = NA
        ctrl_sex_pval[[current_gene]][i] = ctrl_sex_estimate[[current_gene]][i] = fdr_ctrl_sex[[current_gene]][i] = NA 
        KO_sex_pval[[current_gene]][i] = KO_sex_estimate[[current_gene]][i] = fdr_ctrl_sex[[current_gene]][i] = NA
        sep_allKO_pval[[current_gene]][i] = sep_allKO_estimate[[current_gene]][i] = NA
        sep_FvKO_pval[[current_gene]][i] = sep_FvKO_estimate[[current_gene]][i] =NA
        sep_MvKO_pval[[current_gene]][i] = sep_MvKO_estimate [[current_gene]][i]=NA
        #calculate fold-change using WT mice as reference
        foldchange_all[[current_gene]][i]=mean(dataset$PHenoVAlue_raw[dataset$genotype %in% "Gene"], na.rm = T)/mean(dataset$PHenoVAlue_raw[dataset$genotype%in%"Control"],na.rm = T) 
        #calculate fold-change in female mice using WT mice as reference
        foldchange_FvKO[[current_gene]][i]= mean(dataset$PHenoVAlue_raw[dataset$genotype %in% "Gene" & dataset$sex %in% "Female"], na.rm = T)/ mean(dataset$PHenoVAlue_raw[dataset$genotype %in% "Control" & dataset$sex %in% "Female"], na.rm = T)
        #calculate fold-change in male mice using WT mice as reference
        foldchange_MvKO[[current_gene]][i] = mean(dataset$PHenoVAlue_raw[dataset$genotype %in% "Gene" & dataset$sex %in% "Male"], na.rm = T)/ mean(dataset$PHenoVAlue_raw[dataset$genotype %in% "Control" & dataset$sex %in% "Male"], na.rm = T)
      }else{
        # subset data by sex
        u_F <- unique(dataset$PHenoVAlue[dataset$sex %in% "Female"])
        u_M <- unique(dataset$PHenoVAlue[dataset$sex %in% "Male"])
        if(length(u_F) == 1 & length(u_M) == 1){
          Normality_before[[current_gene]][i] = NA
          Normality_after[[current_gene]][i] = NA
          Normality_F[[current_gene]][i] = NA
          Normality_M[[current_gene]][i]= NA
          # count number for WT mice and KO mice
          numb_null[[current_gene]][i]= sum(dataset$genotype %in% "Control")
          numb_gene[[current_gene]][i]= sum(dataset$genotype %in% "Gene")  
        }else if(length(u_F) > 1 & length(u_M) == 1) {
          # test normality before data transformation and after data transformation 
          Normality_before[[current_gene]][i] = shapiro.test(dataset$PHenoVAlue_raw)$p.value
          Normality_after[[current_gene]][i] = shapiro.test(dataset$PHenoVAlue)$p.value
          # test normality for male and female mice
          Normality_F[[current_gene]][i] = shapiro.test(dataset$PHenoVAlue[dataset$sex%in%"Female"] )$p.value
          Normality_M[[current_gene]][i]= NA
          # count the number of WT and KO mice
          numb_null[[current_gene]][i]= sum(dataset$genotype %in% "Control")
          numb_gene[[current_gene]][i]= sum(dataset$genotype %in% "Gene")   
        }else if(length(u_F) == 1 & length(u_M) > 1) {
          # test normality before data transformation and after data transformation 
          Normality_before[[current_gene]][i] = shapiro.test(dataset$PHenoVAlue_raw)$p.value
          Normality_after[[current_gene]][i] = shapiro.test(dataset$PHenoVAlue)$p.value
          # test normality for male and female mice
          Normality_F[[current_gene]][i] = NA
          Normality_M[[current_gene]][i] = shapiro.test(dataset$PHenoVAlue[dataset$sex %in% "Male"] )$p.value
          # count the number of WT and KO mice
          numb_null[[current_gene]][i] = sum(dataset$genotype %in% "Control")
          numb_gene[[current_gene]][i] = sum(dataset$genotype %in% "Gene")   
        }else if(length(u_F) > 1 & length(u_M) > 1) {
          # test normality before data transformation and after data transformation
          Normality_before[[current_gene]][i] = shapiro.test(dataset$PHenoVAlue_raw)$p.value
          Normality_after[[current_gene]][i] = shapiro.test(dataset$PHenoVAlue)$p.value
          # test normality for male and female mice
          Normality_F[[current_gene]][i] = shapiro.test(dataset$PHenoVAlue[dataset$sex %in% "Female"] )$p.value
          Normality_M[[current_gene]][i]= shapiro.test(dataset$PHenoVAlue[dataset$sex %in% "Male"] )$p.value
          # count the number of WT and KO mice
          numb_null[[current_gene]][i]= sum(dataset$genotype %in% "Control")
          numb_gene[[current_gene]][i]= sum(dataset$genotype %in% "Gene")   
        }
      # load null model and full model
      model_formula_null = null_model_genotype(dataset,"PHenoVAlue")
      model_formula_genotype = final_genotype_model(dataset, "PHenoVAlue")
      #test genotype effect
      geno_effect= testing_genotype_effect(dataset,"PHenoVAlue")
      #test hgenotye*sex interaction effect
      formula_interaction_null = null_model_Interaction(dataset,"PHenoVAlue")
      Interct_effect = testing_Interaction_effect(dataset, "PHenoVAlue")
      #test sex effect
      formula_sex_null = null_model_sex(dataset,"PHenoVAlue")
      sex_effect= testing_sex_effect(dataset,"PHenoVAlue")
      #p-values using gls model
      result<-gls(PHenoVAlue ~ genotype + sex + genotype*sex, data = dataset,na.action = 'na.exclude')
      summary<-nlme:::summary.gls(model_afterFIXED)$tTable
        pvalue_Stage1[[current_gene]][i] = geno_effect
        pvalue_Stage1.5[[current_gene]][i] = sex_effect
        pvalue_Stage2[[current_gene]][i] = Interct_effect
        #subset data for WT mice
        ctrl_data = dataset[dataset$genotype%in%"Control"]
        # test sex difference in WT mice using gls() model
        ctrl_test<- tryCatch({
          gls(PHenoVAlue~sex,data = ctrl_data,na.action = 'na.exclude')
        }, error = function(er){
          NA
        }) 
        if(length(ctrl_test) <= 1){
          KO_sex_pval[[current_gene]][i] <- NA
          KO_sex_estimate[[current_gene]][i] <- NA
        }else{
          ctrl_result <- nlme:::summary.gls(ctrl_test)$tTable 
          ctrl_sex_pval[[current_gene]][i] <- ctrl_result[2,4]
          ctrl_sex_estimate[[current_gene]][i] <- ctrl_result[2,1]
        }
        #subset data for KO mice
        KO_data = dataset[dataset$genotype%in%"Gene"]
        # test sex difference in WT mice using gls() model
        KO_test <- tryCatch({
          gls(PHenoVAlue~sex,data = KO_data,na.action = "na.exclude")
        }, error = function(er){
          NA
        }) 
        if(length(KO_test) <= 1){
          KO_sex_pval[[current_gene]][i] <- NA
          KO_sex_estimate[[current_gene]][i] <- NA
        }else{
          KO_result<-nlme:::summary.gls(KO_test)$tTable 
          KO_sex_pval[[current_gene]][i] <- KO_result[2,4]
          KO_sex_estimate[[current_gene]][i] <- KO_result[2,1]
        }
        # test genotyp effect without separating groups by sex
        result_all<- tryCatch({
          gls(PHenoVAlue~genotype,data = dataset,na.action = "na.exclude")
        }, error = function(er){
          NA
        }) 
        if(length(result_all) <= 1){
          sep_allKO_pval[[current_gene]][i] <- NA
          sep_allKO_estimate[[current_gene]][i] <- NA
        }else{
          sep_all <- nlme:::summary.gls(result_all)$tTable 
          sep_allKO_pval[[current_gene]][i] <- sep_all[2,4]
          sep_allKO_estimate[[current_gene]][i] <- sep_all[2,1]
        }
        
        #subset data by sex and individual comparison of genotype effect within sex
        sub_F <- dataset[dataset$sex %in% "Female",]
        sub_M <- dataset[dataset$sex %in% "Male",]
        #genotype as variable in female mice
        result_F <- tryCatch({
          gls(PHenoVAlue ~ genotype,data = sub_F, na.action = "na.exclude")
        }, error = function(er){
          NA
        })
        if(length(result_F) <= 1){
          sep_FvKO_pval[[current_gene]][i] <- NA
          sep_FvKO_estimate[[current_gene]][i] <- NA
        }else{
          sep_FvKO <- nlme:::summary.gls(result_F)$tTable 
          sep_FvKO_pval[[current_gene]][i] <- sep_FvKO[2,4]
          sep_FvKO_estimate[[current_gene]][i] <- sep_FvKO[2,1]
        }
        # genotype as variable in male mice
        result_M <- tryCatch({
          gls(PHenoVAlue ~ genotype, data = sub_M, na.action = "na.exclude")
        }, error = function(er){
          NA
        }) 
        if(length(result_M) <= 1){
          sep_MvKO_pval[[current_gene]][i] <- NA
          sep_MvKO_estimate[[current_gene]][i] <- NA
        }else{
          sep_MvKO <- nlme:::summary.gls(result_M)$tTable 
          sep_MvKO_pval[[current_gene]][i] <- sep_MvKO[2,4]
          sep_MvKO_estimate[[current_gene]][i]<- sep_MvKO[2,1]
        }
        # fold-change calculation using WT mice reference
        foldchange_all[[current_gene]][i] = mean(dataset$PHenoVAlue_raw[dataset$genotype %in% "Gene"], na.rm = T)/mean(dataset$PHenoVAlue_raw[dataset$genotype%in%"Control"],na.rm = T) 
        # fold-change calculation for female KO mice to female WT mice
        foldchange_FvKO[[current_gene]][i] = mean(dataset$PHenoVAlue_raw[dataset$genotype %in% "Gene" & dataset$sex %in% "Female"], na.rm = T)/ mean(dataset$PHenoVAlue_raw[dataset$genotype %in% "Control"& dataset$sex%in%"Female"], na.rm = T)
        # fold-change calculation for male KO mice to male WT mice
        foldchange_MvKO[[current_gene]][i] = mean(dataset$PHenoVAlue_raw[dataset$genotype %in% "Gene" & dataset$sex %in% "Male"], na.rm = T)/ mean(dataset$PHenoVAlue_raw[dataset$genotype %in% "Control" & dataset$sex %in% "Male"], na.rm = T)
              }
             }
        }
  ### adjustment method using ("BH" or its alias "fdr")
  fdr_Stage1[[current_gene]] = p.adjust(pvalue_Stage1[[current_gene]],"fdr")
  fdr_Stage1.5[[current_gene]] = p.adjust(pvalue_Stage1.5[[current_gene]],"fdr")
  fdr_Stage2[[current_gene]] = p.adjust(pvalue_Stage2[[current_gene]],"fdr")
  fdr_ctrl_sex[[current_gene]] = p.adjust(ctrl_sex_pval[[current_gene]],"fdr")
  fdr_KO_sex[[current_gene]] = p.adjust(KO_sex_pval[[current_gene]],"fdr")
  }

# calculate significance ratio in each stage
sapply(pvalue_Stage1, function(x){sum(x < 0.05, na.rm = TRUE)})
sapply(fdr_Stage1, function(x){sum(x < 0.05, na.rm = TRUE)})
sapply(fdr_Stage2, function(x){sum(x < 0.05, na.rm = TRUE)})
#Output result
setwd("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
for(i in 1:length(names(pvalue_Stage1))){
  fwrite(data.table(label = Phenotype_f$label, Phenotype = Phenotype_f$`Parameter name`, Normality_before = Normality_before[[i]], Normality_after = Normality_after[[i]], Normality_F =Normality_F[[i]], Normality_M =Normality_M[[i]], numb_null=numb_null[[i]],numb_gene=numb_gene[[i]], pvalue_Stage1=pvalue_Stage1[[i]],pvalue_stage1.5=pvalue_Stage1.5[[i]], pvalue_stage2=pvalue_Stage2[[i]], fdr_Stage1=fdr_Stage1[[i]], fdr_Stage1.5=fdr_Stage1.5[[i]], fdr_Stage2=fdr_Stage2[[i]], sep_allKO_pval = sep_allKO_pval[[i]], sep_allKO_estimate= sep_allKO_estimate[[i]], sep_FvKO_pval=sep_FvKO_pval[[i]], sep_FvKO_estimate= sep_FvKO_estimate[[i]], sep_MvKO_pval=sep_MvKO_pval[[i]], sep_MvKO_estimate= sep_MvKO_estimate[[i]],sep_FvKO_fdr=p.adjust(sep_FvKO_pval[[i]],'fdr'), sep_MvKO_fdr=p.adjust(sep_MvKO_pval[[i]],'fdr'), foldchange_all= foldchange_all[[i]], foldchange_FvKO = foldchange_FvKO[[i]], foldchange_MvKO = foldchange_MvKO[[i]], ctrl_sex_pval= ctrl_sex_pval[[i]], ctrl_sex_estimate = ctrl_sex_estimate[[i]], fdr_ctrl_sex = fdr_ctrl_sex[[i]], KO_sex_pval= KO_sex_pval[[i]], KO_sex_estimate = KO_sex_estimate[[i]], fdr_KO_sex = fdr_KO_sex[[i]]), paste0("5 TwoStage PHENOtype (",names(pvalue_Stage1)[i],"),genoSD.csv"))
}


#6 information extraction for pie chart: Phenotype / genptype*sex interaction effect and direction
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
  pheno<-read.csv(files_metmerge[i])  
  g<-i
  print(i)
  pheno1<-met[!is.na(met$pvalue_Stage1),]
  pheno1$est_dirct<-"Not significant"
  t<-length(pheno1$pvalue_Stage1)
  pheno1$est_dirct[(pheno1$pvalue_Stage1 < 0.05&pheno1$pvalue_stage2 < 0.05)&(pheno1$sep_FvKO_pval < 0.05&pheno1$sep_MvKO_pval < 0.05 )& ((pheno1$sep_FvKO_estimate>0& pheno1$sep_MvKO_estimate<0) | (pheno1$sep_FvKO_estimate<0 & pheno1$sep_MvKO_estimate>0))] = "Different direction"
  pheno1$est_dirct[(pheno1$pvalue_Stage1 < 0.05&pheno1$pvalue_stage2 < 0.05)&(pheno1$sep_FvKO_pval < 0.05 & pheno1$sep_MvKO_pval < 0.05) & ((pheno1$sep_FvKO_estimate>0 & pheno1$sep_MvKO_estimate>0)  | (pheno1$sep_FvKO_estimate<0 & pheno1$sep_MvKO_estimate<0))] = "Different size"
  pheno1$est_dirct[(pheno1$pvalue_Stage1 < 0.05 & pheno1$pvalue_stage2 < 0.05)&(pheno1$sep_FvKO_pval < 0.05&pheno1$sep_MvKO_pval >= 0.05)] = "One sex_Female only"
  pheno1$est_dirct[(pheno1$pvalue_Stage1 < 0.05 & pheno1$pvalue_stage2 < 0.05)&(pheno1$sep_FvKO_pval >= 0.05&pheno1$sep_MvKO_pval < 0.05)] = "One sex_Male only"
  pheno1$est_dirct[(pheno1$pvalue_Stage1 >= 0.05 & pheno1$pvalue_stage2 < 0.05)&(pheno1$sep_FvKO_pval < 0.05&pheno1$sep_MvKO_pval < 0.05 )& ((pheno1$sep_FvKO_estimate>0& pheno1$sep_MvKO_estimate<0) | (pheno1$sep_FvKO_estimate<0 & pheno1$sep_MvKO_estimate>0))] = "Different direction2"
  pheno1$est_dirct[(pheno1$pvalue_Stage1 >= 0.05 & pheno1$pvalue_stage2 < 0.05)&(pheno1$sep_FvKO_pval < 0.05 & pheno1$sep_MvKO_pval < 0.05) & ((pheno1$sep_FvKO_estimate>0 & pheno1$sep_MvKO_estimate>0)  | (pheno1$sep_FvKO_estimate<0 & pheno1$sep_MvKO_estimate<0))] = "Different size2"
  pheno1$est_dirct[(pheno1$pvalue_Stage1 >= 0.05 & pheno1$pvalue_stage2 < 0.05)&(pheno1$sep_FvKO_pval < 0.05&pheno1$sep_MvKO_pval >= 0.05)] = "One sex_Female only2"
  pheno1$est_dirct[(pheno1$pvalue_Stage1 >= 0.05 & pheno1$pvalue_stage2 < 0.05)&(pheno1$sep_FvKO_pval >= 0.05&pheno1$sep_MvKO_pval < 0.05)] = "One sex_Male only2"
  pheno1$est_dirct[(pheno1$pvalue_Stage1 < 0.05 & pheno1$pvalue_stage2 < 0.05)&(pheno1$sep_FvKO_pval >= 0.05&pheno1$sep_MvKO_pval >= 0.05 )] = "cannot classify"
  pheno1$est_dirct[(pheno1$pvalue_Stage1 >= 0.05 & pheno1$pvalue_stage2 < 0.05)&(pheno1$sep_FvKO_pval >= 0.05&pheno1$sep_MvKO_pval >= 0.05 )] = "cannot classify2"
  pheno1$est_dirct[(pheno1$pvalue_Stage1 < 0.05 & pheno1$pvalue_stage2 >= 0.05)] = "Genotype effect with no sex effect(real)"                 
  
  a<-sum(pheno1$est_dirct == "Different direction")
  a1<-sum(pheno1$est_dirct == "Different direction2")
  b <- sum(pheno1$est_dirct == "Different size")
  b1 <- sum(pheno1$est_dirct == "Different size2")
  c <- sum(pheno1$est_dirct == "One sex_Female only")
  c1 <- sum(pheno1$est_dirct == "One sex_Male only")
  c2 <- sum(pheno1$est_dirct == "One sex_Female only2")
  c3 <- sum(pheno1$est_dirct == "One sex_Male only2")
  d<-sum(pheno1$est_dirct=="cannot classify")
  d1<-sum(pheno1$est_dirct=="cannot classify2")
  e<-sum(pheno1$est_dirct=="Genotype effect with no sex effect(real)")
  n<-sum(pheno1$est_dirct=="Not significant")
  PieChart_Diff_Dirct<-rbind(PieChart_Diff_Dirct,c(g,a,a1,b,b1,c,c1,c2,c3,d,d1,e,n,t))
  colnames(PieChart_Diff_Dirct)<-coltitle 
}
write.csv(PieChart_Diff_Dirct,file = "6 PHENO_classification.csv")

