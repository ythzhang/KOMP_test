library(wcmc)
# setwd("C:\\Users\\yzhan\\Documents\\Fiehn\\3_KOMP\\7_Sexual D_20190605\\Jennly-Zhang-MetaboliteSD_ver03-master")

Raw_Phenotype_UCDavis = wcmc::read_data("C:\\Users\\yzhan\\Documents\\Fiehn\\3_KOMP\\7_Sexual D_20190605\\Jennly Zhang KOMP SD_ver10\\Input\\Raw_Phenotype_UCDavis.xlsx")
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

nrow(Raw_Phenotype_UCDavis_e_merge)

library(tidyverse)
duplicated(Raw_Phenotype_UCDavis_f$label_phenotype)
Phenotype_f_merge<-Raw_Phenotype_UCDavis_f[!duplicated(Raw_Phenotype_UCDavis_f$label_phenotype),]

colnames(Raw_Phenotype_UCDavis_e_merge) = colnames(Raw_Phenotype_UCDavis_e)
# rownames(Raw_Phenotype_UCDavis_e_merge) = Phenotype_f_merge$label_phenotype
rownames(Raw_Phenotype_UCDavis_e_merge) = Phenotype_f_merge$label

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


# ============================================================================================================================

# ========================================================================================================================

Phenotype_UCDavis_e_no_mising = Raw_Phenotype_UCDavis_e_merge
unique_genes = unique(Raw_Phenotype_UCDavis_p$Genotype)
continuous_index = c()
index_cat=c()
n=0

for(g in 1:length(unique_genes)){
    current_gene = unique_genes[g]

    current_label = Raw_Phenotype_UCDavis_p$label[Raw_Phenotype_UCDavis_p$Genotype %in% c(current_gene)]
    current_e = Raw_Phenotype_UCDavis_e_merge[, current_label]

    current_label_F = Raw_Phenotype_UCDavis_p$label[Raw_Phenotype_UCDavis_p$Genotype %in% c(current_gene)&Raw_Phenotype_UCDavis_p$Gender%in%'Female']
    current_e_F = Raw_Phenotype_UCDavis_e_merge[, current_label_F]

    current_label_M = Raw_Phenotype_UCDavis_p$label[Raw_Phenotype_UCDavis_p$Genotype %in% c(current_gene)&Raw_Phenotype_UCDavis_p$Gender%in%'Male']
    current_e_M = Raw_Phenotype_UCDavis_e_merge[,current_label_M]
    
  if (current_gene=='null'){
  for(i in 1:nrow(current_e)){

    continuous = sum(is.na(as.numeric(Raw_Phenotype_UCDavis_e_merge[i,!is.na(Raw_Phenotype_UCDavis_e_merge[i,])]))) == 0


    if(continuous){
      continuous_index[i] = TRUE
      if((sum(is.na(current_e_F[i,]))>=14) & (sum(is.na(current_e_M[i,]))>=14)){
        current_e_F[i,] = NA
        current_e_M[i,] = NA
      }else if((sum(is.na(current_e_F[i,]))>=14) & (sum(is.na(current_e_M[i,]))< 14)){
        current_e_F[i,is.na(current_e_F[i,])] = runif(sum(is.na(current_e_F[i,])),min=0.49,max=0.51) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
        current_e_M[i,] = current_e_M[i,]
        
      }else if((sum(is.na(current_e_F[i,])) < 14) & (sum(is.na(current_e_M[i,])) >= 14)){
        current_e_F[i,] = current_e_F[i,]
        current_e_M[i,is.na(current_e_M[i,])] = runif(sum(is.na(current_e_M[i,])),min=0.49,max=0.51) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
      }else{
        current_e_F[i,] = current_e_F[i,]
        current_e_M[i,] = current_e_M[i,]
      }
    }else{
      continuous_index[i] = FALSE
      n<-n+length(index_cat)
      index_cat[n]<-i
    }
    }}else{
      
      for(i in 1:nrow(current_e)){
        continuous = sum(is.na(as.numeric(Raw_Phenotype_UCDavis_e_merge[i,!is.na(Raw_Phenotype_UCDavis_e_merge[i,])]))) == 0
        
        if(continuous){
          continuous_index[i] = TRUE
          
        if(sum(is.na(current_e[i,]))>=4){
          current_e_F[i,] = NA
          current_e_M[i,] = NA
        }else if(sum(is.na(current_e_F[i,]))==3 & sum(is.na(current_e_M[i,]))==0){
          current_e_F[i,is.na(current_e_F[i,])] = runif(sum(is.na(current_e_F[i,])),min=0.49,max=0.51) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
          current_e_M[i,] = current_e_M[i,]
        }else if(sum(is.na(current_e_F[i,])) ==0 & sum(is.na(current_e_M[i,]))==3){
          current_e_F[i,] = current_e_F[i,]
          current_e_M[i,is.na(current_e_M[i,])] = runif(sum(is.na(current_e_M[i,])),min=0.49,max=0.51) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
        }else if(sum(is.na(current_e_F[i,]))==2 & sum(is.na(current_e_M[i,]))==1){
          current_e_F[i,is.na(current_e_F[i,])] = runif(sum(is.na(current_e_F[i,])),min=0.49,max=0.51) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
          current_e_M[i,is.na(current_e_M[i,])] = runif(sum(is.na(current_e_M[i,])),min=0.98,max=1.02) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
        }else if(sum(is.na(current_e_F[i,]))==1 & sum(is.na(current_e_M[i,]))==2){
          current_e_F[i,is.na(current_e_F[i,])] = runif(sum(is.na(current_e_F[i,])),min=0.98,max=1.02) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
          current_e_M[i,is.na(current_e_M[i,])] = runif(sum(is.na(current_e_M[i,])),min=0.49,max=0.51) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
        }else if(sum(is.na(current_e[i,]))==2){
          if(sum(is.na(current_e_F[i,]))==2 & sum(is.na(current_e_M[i,]))==0){
            current_e_F[i,is.na(current_e_F[i,])] = runif(sum(is.na(current_e_F[i,])),min=0.98,max=1.02) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
            current_e_M[i,] = current_e_M[i,]
          }else if(sum(is.na(current_e_F[i,]))==0 & sum(is.na(current_e_M[i,]))==2){
            current_e_F[i,] = current_e_F[i,]
            current_e_M[i,is.na(current_e_M[i,])] = runif(sum(is.na(current_e_M[i,])),min=0.98,max=1.02) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
          }else if(sum(is.na(current_e_F[i,]))==1 & sum(is.na(current_e_M[i,]))==1){
            current_e_F[i,is.na(current_e_F[i,])] = runif(sum(is.na(current_e_F[i,])),min=0.98,max=1.02) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
            current_e_M[i,is.na(current_e_M[i,])] = runif(sum(is.na(current_e_M[i,])),min=0.98,max=1.02) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
          }
        }else if(sum(is.na(current_e[i,]))<=1){
          current_e_F[i,is.na(current_e_F[i,])] = runif(sum(is.na(current_e_F[i,])),min=0.98,max=1.02) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
          current_e_M[i,is.na(current_e_M[i,])] = runif(sum(is.na(current_e_M[i,])),min=0.98,max=1.02) * min(as.numeric(current_e[i,!is.na(current_e[i,])]))
        }
      }else{
        continuous_index[i] = FALSE
        n<-n+length(index_cat)
        index_cat[n]<-i
      }
    }
    }
  Phenotype_UCDavis_e_no_mising[,current_label_F] = current_e_F
  Phenotype_UCDavis_e_no_mising[,current_label_M] = current_e_M
}


Phenotype_e_no_mising = Phenotype_UCDavis_e_no_mising[continuous_index,]
Phenotype_e_no_mising = as.data.frame(unlist(Phenotype_e_no_mising))
# str(Phenotype_e_no_mising)
# Phenotype_e_no_mising_new <- matrix(Phenotype_e_no_mising)
# summary(Phenotype_e_no_mising_new)
# str(Phenotype_e_no_mising_new)
# Phenotype_e_no_mising_new[2,1] <- as.numeric (as.character( as.factor(unlist(Phenotype_e_no_mising[2,1])))) 
# Phenotype_e_no_mising <- as.data.frame(sapply(Phenotype_e_no_mising, as.numeric))
# # str(Phenotype_e_no_mising)
# str(Phenotype_e_no_mising[1,1])

Phenotype_e_no_mising_new <- matrix(nrow= nrow(Phenotype_e_no_mising), ncol=ncol(Phenotype_e_no_mising))
colnames(Phenotype_e_no_mising_new) <- colnames(Phenotype_e_no_mising)
rownames(Phenotype_e_no_mising_new) <- rownames(Phenotype_e_no_mising)



for (i in 1: nrow(Phenotype_e_no_mising_new)){
  
  for (c in 1: length(Phenotype_e_no_mising_new[i,])){
    # Phenotype_e_no_mising_new[i,c]  <- NA
    # row_value[c] <- as.numeric(unlist(as.numeric(Phenotype_e_no_mising_new[i,c])) )
    Phenotype_e_no_mising_new[i,c] <- as.numeric (as.character( as.factor(unlist(Phenotype_e_no_mising[i,c])))) 
    # Phenotype_e_no_mising[i,c] <- as.numeric (as.character(Phenotype_e_no_mising[i,c]))
  }
  
}

# Phenotype_e_no_mising[i,] <- as.numeric(unlist (Phenotype_e_no_mising[i,]))

str(Phenotype_e_no_mising_new)
# Phenotype_e_no_mising = as.numeric(as.character(Phenotype_e_no_mising))
Phenotype_f = Phenotype_f_merge[continuous_index,]
Phenotype_p = Raw_Phenotype_UCDavis_p


null_label = Phenotype_p$label[Phenotype_p$Genotype %in% c('null')]
null_e = Phenotype_e_no_mising_new[, null_label]


null_label = Phenotype_p$label[Phenotype_p$Genotype %in% c('null') & Phenotype_p$Gender %in% "Female"]
null_e_female = Phenotype_e_no_mising_new[, null_label]
null_label = Phenotype_p$label[Phenotype_p$Genotype %in% c('null') & Phenotype_p$Gender %in% "Male"]
null_e_male = Phenotype_e_no_mising_new[, null_label]

# =========================================================================================================================================
# =========================================================================================================================================
# 5 for each genotype, check which metabolite are associated with gene * gender, gene, gender, or not significant at all.



# --------------------------------------------------------------------------------------------------------


# https://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0052410.s004
# https://www.researchgate.net/profile/Stano_Pekar/publication/304001371_Marginal_Models_Via_GLS_A_Convenient_Yet_Neglected_Tool_for_the_Analysis_of_Correlated_Data_in_the_Behavioural_Sciences/links/5779816e08aeb9427e2c0017/Marginal-Models-Via-GLS-A-Convenient-Yet-Neglected-Tool-for-the-Analysis-of-Correlated-Data-in-the-Behavioural-Sciences.pdf
# https://www.zoology.ubc.ca/~schluter/R/fit-model/
####### https://www.youtube.com/watch?v=vN5cNN2-HWE

#load package nlme
library(nlme)
# library(car)
# Function 1 function to build model for testing of fixed effects
# Assumption: No batch effect and weight is not a independet variance 
# Goal: 

model_forFIXEDtest<-function(dataset, PHenoVAlue){
  
  model.formula <- as.formula(paste(PHenoVAlue, "~", paste("genotype", "sex", "genotype*sex", sep= "+")))
  model=gls(model.formula, dataset, na.action="na.omit")
  # model=lm(model.formula, dataset, na.action="na.omit")
}



# -------------------------------------------------------------------------------------------------------------------------------
#Function 2: testing the fixed effects and building genotype model formula
# Goal: Genotype is always include in this model. The goal is to test genotype and sex effect and based on the output build the final genotype model formula for later testing.  
# anova function: null hypothesis, the regression coefficients are equal to zero  
#                 alternative hypothesis that the regression coefficient are not equal to zero.
# the null hypothesis was rejected when p-values < 0.05 (e.g.,accept the alternative hypothesis, the components of the model should be included in later analysis)
# Note a complexity surrounds the interaction term  - if it is significant but gender is excluded it is excluded. 


final_genotype_model<-function(dataset, PHenoVAlue){
  model_afterFIXED=model_forFIXEDtest(dataset, PHenoVAlue)
  anova_results = anova(model_afterFIXED, type="marginal")$"p-value" < 0.05
  keepSex = anova_results[3]
  keepInteraction = anova_results[4]
  
  if(keepSex && keepInteraction){
    return(model.formula <- as.formula(paste(PHenoVAlue, "~", paste("genotype", "sex", "genotype*sex", sep= "+"))))
  }else if(keepSex &&  !keepInteraction){
    return(model.formula <- as.formula(paste(PHenoVAlue, "~", paste("genotype", "sex", sep= "+"))))
  }else if(!keepSex &&  !keepInteraction){
    return(model.formula <- as.formula(paste(PHenoVAlue, "~", "genotype")))
  }else if(!keepSex  && keepInteraction){
    return(model.formula <- as.formula(paste(PHenoVAlue, "~", paste("genotype", "sex", "genotype*sex", sep= "+"))))
  }
}


#---------------------------------------------------------------------------------------------------------------------------
#Function 3: testing the fixed effects and building final Genotype effect null model
# Goal:  to test fixed effects of the model and based on the output build the final null model formula for later testing - as a null model it automatically excludes genotype and interaction term. 
# The anova function tests the fixed effects associated by treatment with a null hypothesis that the regression coefficients are equal to zero  and an alternative hypothesis that the regression coefficient are not equal to zero.
# If the p-values of these tests are less than 0.05 we reject the null and accept the alternative that the are significant components of the model and should be included. 
# If no terms are significant a model can be build with just an intercept element this is specified as  "model.formula <- as.formula(paste(depVariable, "~", "1"))"


null_model_genotype<-function(dataset, PHenoVAlue){ 
  model_afterFIXED=model_forFIXEDtest(dataset, PHenoVAlue)
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
null_model_Interaction<-function(dataset, PHenoVAlue){ 
  model_afterFIXED=model_forFIXEDtest(dataset, PHenoVAlue)
  anova_results = anova(model_afterFIXED, type="marginal")$"p-value" < 0.05
  keepSex = anova_results[3]
  keepGenotype = anova_results[2]     
  
  if(!keepSex && !keepGenotype){
    return(model.formula <- as.formula(paste(PHenoVAlue, "~", "1")))
  }else{
    return(model.formula <- as.formula(paste(PHenoVAlue, "~", paste("genotype","sex",sep="+"))))
  } 
}



# ------------------------------------------------------------------------------------------------

#Function 4: testing the genotype effect
# The goal of the function is to give a pvalue of whether genotype effect is significant by compare the genotype model with the null model with the anova function
# The model_formula_null and model_formula_genotype are called to define the models for comparison.
# For each possible combination,  then an anova model is used to report the pvalue.
# For testing the genotype effect we use method= ML


testing_genotype_effect<-function(dataset,PHenoVAlue){
  #call functions to determine model.formula
  model_formula_null = null_model_genotype(dataset,PHenoVAlue)
  model_formula_genotype = final_genotype_model(dataset,PHenoVAlue)
  model_genotype=gls(model_formula_genotype,dataset,method='ML', na.action="na.omit")
  # model_genotype=lm(model_formula_genotype,dataset,na.action="na.omit")
  
  model_null=gls(model_formula_null, dataset,method='ML',na.action="na.omit")
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
  model_interaction_full=gls(PHenoVAlue~genotype+sex+genotype*sex,dataset, method='ML',na.action="na.omit")
  model_null_interaction =gls(formula_interaction_null, dataset,method='ML', na.action="na.omit")
  return(pvalue_Stage2=(anova(model_interaction_full, model_null_interaction)$`p-value`[2]))
}


#Function 4-d: testing the sex effect
null_model_sex<-function(dataset, PHenoVAlue){ 
  model_afterFIXED=model_forFIXEDtest(dataset, PHenoVAlue)
  anova_results = anova(model_afterFIXED, type="marginal")$"p-value" < 0.05
  keepSex = anova_results[3]
  keepGenotype = anova_results[2]     
  
  if(!keepSex && !keepGenotype){
    return(model.formula <- as.formula(paste(PHenoVAlue, "~", "1")))
  }else{
    return(model.formula <- as.formula(paste(PHenoVAlue, "~", "genotype")))
  } 
}

# ------------------------------------------------------------------------------------------------

testing_sex_effect<-function(dataset, PHenoVAlue){
  #call functions to determine model.formula
  formula_sex_null = null_model_sex(dataset,"PHenoVAlue")
  model_sex_full=gls(PHenoVAlue~genotype+sex,dataset, method='ML',na.action="na.omit")
  model_null_sex =gls(formula_sex_null, dataset,method='ML', na.action="na.omit")
  return(pvalue_Stage1.5=(anova(model_sex_full, model_null_sex)$`p-value`[2]))
}





#-----------------------------------------------------------------------------------------
#Function 5: following function returns the final model which is needed for the diagnostic plots 
#Goal of this function is to return a model output that is the final model following all the previous tests which can be used to generate diagnostics plots and output the final model details. 
# The model_formula_null and model_formula_genotype are called to define the models for comparison. 
# The keep_batch and keep_equal variance function are called this allows if but rules to build the correct model comparison ie if batch is included then we use a mixed model etc
# For each possible combination,  then the model is fitted to the data and the model reported as the output.
#this is a separate function to above as for the estimates of the fixed effects we use method=REML whilst for the genotype test we used method= ML   
#na.action = na.exclude as in this form the residue calculated from the model output have the same length as the original datafile
#see http://stackoverflow.com/questions/6882709/how-do-i-deal-with-nas-in-residuals-in-a-regression-in-r


finalmodel<-function(dataset, PHenoVAlue){
  #call functions to determine model.formula
  model_formula_null = null_model_genotype(dataset,PHenoVAlue)
  model_genotype = final_genotype_model(dataset, PHenoVAlue)
  
  model_genotype=gls(model_genotype,dataset, na.action="na.exclude")
}


# ---------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------
#Starting TwoStage significance calculation

Normality_before = Normality_after =Normality_F = Normality_M= list()

numb_null = numb_gene = list()

pvalue_Stage1 = pvalue_Stage1.5 = pvalue_Stage2 = list()
fdr_Stage1 = fdr_Stage1.5 = fdr_Stage2 = list()

ctrl_sex_pval = ctrl_sex_estimate = fdr_ctrl_sex = list() 
KO_sex_pval = KO_sex_estimate = fdr_KO_sex = list()

sep_allKO_pval = sep_allKO_estimate =  list()
sep_FvKO_pval = sep_FvKO_estimate = sep_MvKO_pval = sep_MvKO_estimate = list()
foldchange_all = foldchange_FvKO = foldchange_MvKO = list()



for(g in 2:length(unique_genes)){
  current_gene = unique_genes[g]
  print(g)
  current_label = Phenotype_p$label[Phenotype_p$Genotype %in% c(current_gene)]
  current_e = Phenotype_e_no_mising_new[, current_label]
  current_e=as.data.frame(current_e)
  current_p = Phenotype_p[Phenotype_p$label %in% current_label,]
  
  null_p = Phenotype_p[Phenotype_p$Genotype %in% 'null',]
  
  # sex_index[[current_gene]]=c()
  numb_null[[current_gene]]=numb_gene[[current_gene]]=c()
  
  Normality_before[[current_gene]] = Normality_after[[current_gene]] =Normality_F[[current_gene]] = Normality_M[[current_gene]]= c()
  
  pvalue_Stage1[[current_gene]]=pvalue_Stage1.5[[current_gene]]=pvalue_Stage2[[current_gene]]=c()
  
  fdr_Stage1[[current_gene]]= fdr_Stage1.5[[current_gene]] =fdr_Stage2[[current_gene]] =c()
  
  ctrl_sex_pval[[current_gene]] = ctrl_sex_estimate[[current_gene]] = fdr_ctrl_sex[[current_gene]] = c() 
  KO_sex_pval[[current_gene]] = KO_sex_estimate[[current_gene]] = fdr_KO_sex[[current_gene]] = c()
  
  sep_allKO_pval[[current_gene]] = sep_allKO_estimate[[current_gene]] = fdr_allKO = c()
  
  
  sep_FvKO_pval[[current_gene]] = sep_FvKO_estimate[[current_gene]] = sep_MvKO_pval[[current_gene]] = sep_MvKO_estimate[[current_gene]] = c()
  
  foldchange_all[[current_gene]] = c()
  foldchange_FvKO[[current_gene]] = c()
  foldchange_MvKO[[current_gene]] = c()
  
  
  for(i in 1:nrow(Phenotype_e_no_mising)){
    
    if((sum(is.na(current_e[i,])) >= (ncol(current_e)-2) )|(sum(is.na(null_e[i,])) >= (ncol(null_e)-6) )){
      numb_null[[current_gene]][i]=numb_gene[[current_gene]][i]=NA
      Normality_before[[current_gene]][i] =NA
      Normality_after[[current_gene]][i] =NA
      Normality_F[[current_gene]][i] =NA
      Normality_M[[current_gene]][i] =NA
      
      
      pvalue_Stage1[[current_gene]][i]=NA
      pvalue_Stage1.5[[current_gene]][i]=NA
      pvalue_Stage2[[current_gene]][i]=NA
      
      ctrl_sex_pval[[current_gene]][i] = ctrl_sex_estimate[[current_gene]][i] = fdr_ctrl_sex[[current_gene]][i] = NA
      KO_sex_pval[[current_gene]][i] = KO_sex_estimate[[current_gene]][i] = NA
      
      sep_allKO_pval[[current_gene]][i] = sep_allKO_estimate[[current_gene]][i] = NA
      sep_FvKO_pval[[current_gene]][i]= sep_FvKO_estimate[[current_gene]][i] =NA
      sep_MvKO_pval[[current_gene]][i]= sep_MvKO_estimate[[current_gene]][i] =NA
      
      
      foldchange_all[[current_gene]][i]= NA
      foldchange_FvKO[[current_gene]][i]= NA
      foldchange_MvKO[[current_gene]][i]= NA
      
      
    }else{
      
      
      data_raw<- data.table(y =(c(current_e[i,],null_e[i,])), group = rep(c("Gene","Control"), c(ncol(current_e), ncol(null_e))), sex=c(current_p$Gender,null_p$Gender))
      
      colnames(data_raw) = c("PHenoVAlue_raw","genotype",'sex')
      data_raw$PHenoVAlue_raw <-unlist(data_raw$PHenoVAlue_raw)
      # data_raw<-(data_raw$PHenoVAlue)
      data_raw$PHenoVAlue_raw<-as.numeric(as.character(data_raw$PHenoVAlue_raw))
      data_raw<-data_raw[!is.na(data_raw$PHenoVAlue_raw)]
      
      
      data_raw$PHenoVAlue<-rankNorm(data_raw$PHenoVAlue_raw)
      
      # dataset<-data_raw[!is.na(data_raw$PHenoVAlue),]
      
      # mean(data_raw$PHenoVAlue)
      
      outliers_nullF_raw<-boxplot.stats(data_raw$PHenoVAlue[data_raw$sex%in%'Female'&data_raw$genotype%in%'Control'])$out

      if(length(outliers_nullF_raw)==0){
        dataset2_raw<-data_raw
      }else{
        dataset2_raw<-data_raw[-which((data_raw$PHenoVAlue%in% outliers_nullF_raw)&(data_raw$sex%in%'Female')&(data_raw$genotype%in%'Control')),]
      }

      outliers_nullM_raw<-boxplot.stats(dataset2_raw$PHenoVAlue[dataset2_raw$sex%in%'Male'&dataset2_raw$genotype%in%'Control'])$out
      if(length(outliers_nullM_raw)==0){
        dataset3_raw<-dataset2_raw
      }else{
        dataset3_raw<-dataset2_raw[-which((dataset2_raw$PHenoVAlue%in% outliers_nullM_raw)&(dataset2_raw$sex%in%'Male')&(dataset2_raw$genotype%in%'Control')),]
      }


      outliers_genoF_raw<-boxplot.stats(dataset3_raw$PHenoVAlue[dataset3_raw$sex%in%'Female'&dataset3_raw$genotype%in%'Gene'])$out
      if(length(outliers_genoF_raw)==0){
        dataset4_raw<-dataset3_raw
      }else{
        dataset4_raw<-dataset3_raw[-which((dataset3_raw$PHenoVAlue%in% outliers_genoF_raw)&(dataset3_raw$sex%in%'Female')&(dataset3_raw$genotype%in%'Gene')),]
      }

      outliers_genoM_raw<-boxplot.stats(dataset4_raw$PHenoVAlue[dataset4_raw$sex%in%'Male'&dataset4_raw$genotype%in%'Control'])$out
      if(length(outliers_genoM_raw)==0){
        dataset_raw<-dataset4_raw
      }else{
        dataset_raw<-dataset4_raw[-which((dataset4_raw$PHenoVAlue%in% outliers_genoM_raw)&(dataset4_raw$sex%in%'Male')&(dataset4_raw$genotype%in%'Control')),]
      }
      
      dataset<-dataset_raw[!is.na(dataset_raw$PHenoVAlue),]
      
      
      model_afterFIXED=tryCatch({
        model_forFIXEDtest(dataset, "PHenoVAlue")
      }, error = function(er){
        NA
      })
      if(length(model_afterFIXED)<=1){
        
        numb_null[[current_gene]][i]= sum(dataset$genotype%in%"Control")
        numb_gene[[current_gene]][i]= sum(dataset$genotype%in%"Gene")
        
        Normality_before[[current_gene]][i] =NA
        Normality_after[[current_gene]][i] =NA
        Normality_F[[current_gene]][i] =NA
        Normality_M[[current_gene]][i] =NA
        
        pvalue_Stage1[[current_gene]][i] = NA
        pvalue_Stage1.5[[current_gene]][i] = NA
        pvalue_Stage2[[current_gene]][i] = NA
        
        ctrl_sex_pval[[current_gene]][i] = ctrl_sex_estimate[[current_gene]][i] = fdr_ctrl_sex[[current_gene]][i] = NA 
        KO_sex_pval[[current_gene]][i] = KO_sex_estimate[[current_gene]][i] = fdr_ctrl_sex[[current_gene]][i] = NA
        
        sep_allKO_pval[[current_gene]][i] = sep_allKO_estimate[[current_gene]][i] = NA
        sep_FvKO_pval[[current_gene]][i] = sep_FvKO_estimate[[current_gene]][i] =NA
        sep_MvKO_pval[[current_gene]][i] = sep_MvKO_estimate [[current_gene]][i]=NA
        
        
        foldchange_all[[current_gene]][i]=mean(dataset$PHenoVAlue_raw[dataset$genotype%in%"Gene"],na.rm = T)/mean(dataset$PHenoVAlue_raw[dataset$genotype%in%"Control"],na.rm = T) 
        foldchange_FvKO[[current_gene]][i]= mean(dataset$PHenoVAlue_raw[dataset$genotype%in%"Gene"&dataset$sex%in%"Female"],na.rm = T)/ mean(dataset$PHenoVAlue_raw[dataset$genotype%in%"Control"&dataset$sex%in%"Female"],na.rm = T)
        foldchange_MvKO[[current_gene]][i]= mean(dataset$PHenoVAlue_raw[dataset$genotype%in%"Gene"&dataset$sex%in%"Male"],na.rm = T)/ mean(dataset$PHenoVAlue_raw[dataset$genotype%in%"Control"&dataset$sex%in%"Male"],na.rm = T)
        
        
        
      }else{
        
        
        u_F<-unique(dataset$PHenoVAlue[dataset$sex%in%"Female"])
        u_M<-unique(dataset$PHenoVAlue[dataset$sex%in%"Male"])
        
        if(length(u_F)==1 &length(u_M)==1){
          Normality_before[[current_gene]][i] = NA
          Normality_after[[current_gene]][i] = NA
          Normality_F[[current_gene]][i] = NA
          Normality_M[[current_gene]][i]= NA
          
          numb_null[[current_gene]][i]= sum(dataset$genotype %in%"Control")
          numb_gene[[current_gene]][i]= sum(dataset$genotype %in%"Gene")  
          
        }else if(length(u_F)>1&length(u_M)==1) {
          Normality_before[[current_gene]][i] = shapiro.test(dataset$PHenoVAlue_raw)$p.value
          Normality_after[[current_gene]][i] = shapiro.test(dataset$PHenoVAlue)$p.value
          Normality_F[[current_gene]][i] = shapiro.test(dataset$PHenoVAlue[dataset$sex%in%"Female"] )$p.value
          Normality_M[[current_gene]][i]= NA
          
          numb_null[[current_gene]][i]= sum(dataset$genotype %in%"Control")
          numb_gene[[current_gene]][i]= sum(dataset$genotype %in%"Gene")   
        }else if(length(u_F)==1&length(u_M)>1) {
          Normality_before[[current_gene]][i] = shapiro.test(dataset$PHenoVAlue_raw)$p.value
          Normality_after[[current_gene]][i] = shapiro.test(dataset$PHenoVAlue)$p.value
          Normality_F[[current_gene]][i] = NA
          Normality_M[[current_gene]][i]= shapiro.test(dataset$PHenoVAlue[dataset$sex%in%"Male"] )$p.value
          
          numb_null[[current_gene]][i]= sum(dataset$genotype %in%"Control")
          numb_gene[[current_gene]][i]= sum(dataset$genotype %in%"Gene")   
        }else if(length(u_F)>1&length(u_M)>1) {
          Normality_before[[current_gene]][i] = shapiro.test(dataset$PHenoVAlue_raw)$p.value
          Normality_after[[current_gene]][i] = shapiro.test(dataset$PHenoVAlue)$p.value
          Normality_F[[current_gene]][i] = shapiro.test(dataset$PHenoVAlue[dataset$sex%in%"Female"] )$p.value
          Normality_M[[current_gene]][i]= shapiro.test(dataset$PHenoVAlue[dataset$sex%in%"Male"] )$p.value
          
          numb_null[[current_gene]][i]= sum(dataset$genotype %in%"Control")
          numb_gene[[current_gene]][i]= sum(dataset$genotype %in%"Gene")   
        }
        
        


      
      model_formula_null=null_model_genotype(dataset,"PHenoVAlue")
      model_formula_genotype=final_genotype_model(dataset, "PHenoVAlue")
      
      geno_effect= testing_genotype_effect(dataset,"PHenoVAlue")
      
      #pvalue_genotype effect
      formula_interaction_null=null_model_Interaction(dataset,"PHenoVAlue")
      Interct_effect= testing_Interaction_effect(dataset,"PHenoVAlue")
      
      formula_sex_null = null_model_sex(dataset,"PHenoVAlue")
      sex_effect= testing_sex_effect(dataset,"PHenoVAlue")
      #pvalue_SD
      result<-gls(PHenoVAlue~genotype+sex+genotype*sex,data = dataset,na.action = 'na.exclude')
      summary<-nlme:::summary.gls(model_afterFIXED)$tTable
      
      # if(nrow(summary) == 4){
       
        pvalue_Stage1[[current_gene]][i]=geno_effect
        pvalue_Stage1.5[[current_gene]][i]=sex_effect
        pvalue_Stage2[[current_gene]][i]=Interct_effect
        
        
        
        ctrl_data = dataset[dataset$genotype%in%"Control"]
        
        ctrl_test<- tryCatch({
          gls(PHenoVAlue~sex,data = ctrl_data,na.action = 'na.exclude')
        }, error = function(er){
          NA
        }) 
        if(length(ctrl_test)<=1){
          KO_sex_pval[[current_gene]][i]<- NA
          KO_sex_estimate[[current_gene]][i]<- NA
        }else{
          ctrl_result<-nlme:::summary.gls(ctrl_test)$tTable 
          ctrl_sex_pval[[current_gene]][i]<- ctrl_result[2,4]
          ctrl_sex_estimate[[current_gene]][i]<- ctrl_result[2,1]
        }
        
        
        KO_data = dataset[dataset$genotype%in%"Gene"]
        
        KO_test<- tryCatch({
          gls(PHenoVAlue~sex,data = KO_data,na.action = 'na.exclude')
        }, error = function(er){
          NA
        }) 
        if(length(KO_test)<=1){
          KO_sex_pval[[current_gene]][i]<- NA
          KO_sex_estimate[[current_gene]][i]<- NA
        }else{
          KO_result<-nlme:::summary.gls(KO_test)$tTable 
          KO_sex_pval[[current_gene]][i]<- KO_result[2,4]
          KO_sex_estimate[[current_gene]][i]<- KO_result[2,1]
        }
        
        
        
        result_all<- tryCatch({
          gls(PHenoVAlue~genotype,data = dataset,na.action = 'na.exclude')
        }, error = function(er){
          NA
        }) 
        
        if(length(result_all)<=1){
          sep_allKO_pval[[current_gene]][i]<- NA
          sep_allKO_estimate[[current_gene]][i]<- NA
        }else{
          sep_all<-nlme:::summary.gls(result_all)$tTable 
          sep_allKO_pval[[current_gene]][i] <- sep_all[2,4]
          sep_allKO_estimate[[current_gene]][i]<- sep_all[2,1]
        }
        
        
        
        
        sub_F<- dataset[dataset$sex%in%"Female",]
        sub_M<-dataset[dataset$sex%in%"Male",]
        
        result_F<-tryCatch({
          gls(PHenoVAlue~genotype,data = sub_F,na.action = 'na.exclude')
        }, error = function(er){
          NA
        })
        if(length(result_F)<=1){
          sep_FvKO_pval[[current_gene]][i]<- NA
          sep_FvKO_estimate[[current_gene]][i]<- NA
        }else{
          sep_FvKO<-nlme:::summary.gls(result_F)$tTable 
          sep_FvKO_pval[[current_gene]][i]<- sep_FvKO[2,4]
          sep_FvKO_estimate[[current_gene]][i]<- sep_FvKO[2,1]
        }
        
        result_M<- tryCatch({
          gls(PHenoVAlue~genotype,data = sub_M,na.action = 'na.exclude')
        }, error = function(er){
          NA
        }) 
        if(length(result_M)<=1){
          sep_MvKO_pval[[current_gene]][i]<- NA
          sep_MvKO_estimate[[current_gene]][i]<- NA
        }else{
          sep_MvKO<-nlme:::summary.gls(result_M)$tTable 
          sep_MvKO_pval[[current_gene]][i]<- sep_MvKO[2,4]
          sep_MvKO_estimate[[current_gene]][i]<- sep_MvKO[2,1]
        }
        
        foldchange_all[[current_gene]][i]=mean(dataset$PHenoVAlue_raw[dataset$genotype%in%"Gene"],na.rm = T)/mean(dataset$PHenoVAlue_raw[dataset$genotype%in%"Control"],na.rm = T) 
        foldchange_FvKO[[current_gene]][i]= mean(dataset$PHenoVAlue_raw[dataset$genotype%in%"Gene"&dataset$sex%in%"Female"],na.rm = T)/ mean(dataset$PHenoVAlue_raw[dataset$genotype%in%"Control"&dataset$sex%in%"Female"],na.rm = T)
        foldchange_MvKO[[current_gene]][i]= mean(dataset$PHenoVAlue_raw[dataset$genotype%in%"Gene"&dataset$sex%in%"Male"],na.rm = T)/ mean(dataset$PHenoVAlue_raw[dataset$genotype%in%"Control"&dataset$sex%in%"Male"],na.rm = T)
        
              }
             }
        }
  
  fdr_Stage1[[current_gene]] = p.adjust(pvalue_Stage1[[current_gene]],'fdr')
  fdr_Stage1.5[[current_gene]] = p.adjust(pvalue_Stage1.5[[current_gene]],'fdr')
  fdr_Stage2[[current_gene]] = p.adjust(pvalue_Stage2[[current_gene]],'fdr')
  
  fdr_ctrl_sex[[current_gene]] = p.adjust(ctrl_sex_pval[[current_gene]],'fdr')
  fdr_KO_sex[[current_gene]] = p.adjust(KO_sex_pval[[current_gene]],'fdr')
  }


sapply(pvalue_Stage1,function(x){sum(x<0.05,na.rm = TRUE)})
sapply(fdr_Stage1,function(x){sum(x<0.05,na.rm = TRUE)})
sapply(fdr_Stage1.5,function(x){sum(x<0.05,na.rm = TRUE)})
sapply(fdr_Stage2,function(x){sum(x<0.05,na.rm = TRUE)})



setwd("C:\\Users\\yzhan\\Documents\\Fiehn\\3_KOMP\\7_Sexual D_20190605\\Jennly Zhang KOMP SD_ver10\\Output\\20200410_PHENOtype_genoSD")
for(i in 1:length(names(pvalue_Stage1))){
  fwrite(data.table(label = Phenotype_f$label_phenotype, Phenotype=Phenotype_f$Parameter,Normality_before = Normality_before[[i]], Normality_after = Normality_after[[i]], Normality_F =Normality_F[[i]],  Normality_M =Normality_M[[i]] , numb_null=numb_null[[i]],numb_gene=numb_gene[[i]], pvalue_Stage1=pvalue_Stage1[[i]],pvalue_stage1.5=pvalue_Stage1.5[[i]], pvalue_stage2=pvalue_Stage2[[i]], fdr_Stage1=fdr_Stage1[[i]], fdr_Stage1.5=fdr_Stage1.5[[i]], fdr_Stage2=fdr_Stage2[[i]], sep_allKO_pval = sep_allKO_pval[[i]], sep_allKO_estimate= sep_allKO_estimate[[i]], sep_FvKO_pval=sep_FvKO_pval[[i]], sep_FvKO_estimate= sep_FvKO_estimate[[i]], sep_MvKO_pval=sep_MvKO_pval[[i]], sep_MvKO_estimate= sep_MvKO_estimate[[i]],sep_FvKO_fdr=p.adjust(sep_FvKO_pval[[i]],'fdr'), sep_MvKO_fdr=p.adjust(sep_MvKO_pval[[i]],'fdr'), foldchange_all= foldchange_all[[i]], foldchange_FvKO = foldchange_FvKO[[i]], foldchange_MvKO = foldchange_MvKO[[i]], ctrl_sex_pval= ctrl_sex_pval[[i]], ctrl_sex_estimate = ctrl_sex_estimate[[i]], fdr_ctrl_sex = fdr_ctrl_sex[[i]], KO_sex_pval= KO_sex_pval[[i]], KO_sex_estimate = KO_sex_estimate[[i]], fdr_KO_sex = fdr_KO_sex[[i]]), paste0("5 TwoStage PHENOtype (",names(pvalue_Stage1)[i],"),genoSD_0410.csv"))
}

  


# =======================================================================================================
# =======================================================================================================
# =======================================================================================================
# =======================================================================================================

#5d Figure 4 information extraction pie chart-metabolite~gene&gender&direction

# ==============================================================================================================================
# End of run using 
# RAW pvalue
# ==============================================================================================================================
# ==============================================================================================================================
# ==============================================================================================================================
# ==============================================================================================================================
# ==============================================================================================================================

# sTART of run using 0.02 RAW P value as threshold

##### calculate the metabolites in each significant group  pie chart gene*sex interaction

setwd("C:\\Users\\yzhan\\Documents\\Fiehn\\3_KOMP\\7_Sexual D_20190605\\Jennly Zhang KOMP SD_ver10\\Output\\20200410_PHENOtype_genoSD")
# files_cmb<-list.files(pattern="^5 COMBINED DATA for each genotype (.*)check the direction in male and female",full.names = T)
files_metmerge<-list.files(pattern="^5 TwoStage PHENOtype (.*),genoSD_0410",full.names = T)

GeneList<-read.csv("C:\\Users\\yzhan\\Documents\\Fiehn\\3_KOMP\\7_Sexual D_20190605\\Jennly Zhang KOMP SD_ver10\\Input\\GeneList_pheno.csv",header = TRUE)
GeneList<-(GeneList$Genotype)
Genelist<-as.character(GeneList)


# GeneList<-as.matrix(GeneList)
PieChart_Diff_Dirct<-data.frame(matrix(ncol = 10, nrow = 0))

coltitle<-c('Genotype','Different direction','Different size','One sex', "Cannot Classify", "Genotype effect with no sex effect(real)", "Genotype effect with no sex effect(cannot classify)","interaction effect with no genotype effect",'Not Significant',"Phenotype count")
colnames(PieChart_Diff_Dirct)<-coltitle




for (i in 1:length(Genelist)){
  met<-read.csv(files_metmerge[i])  #met<-read.csv(files_metmerge[1])
  g<-i
  print(i)
  met1<-met[!is.na(met$pvalue_Stage1),]
  met1$est_dirct<-"Not significant"
  t<-length(met1$pvalue_Stage1)
  # print(t)
  # for (k in 1:length(met1$pvalue_Stage1)){
  met1$est_dirct[(met1$pvalue_Stage1<=0.05&met1$pvalue_stage2<=0.05)&(met1$sep_FvKO_pval<=0.05&met1$sep_MvKO_pval<=0.05)& ((met1$sep_FvKO_estimate>0& met1$sep_MvKO_estimate<0) | (met1$sep_FvKO_estimate<0 &met1$sep_MvKO_estimate>0))] = "Different direction"
  
  met1$est_dirct[(met1$pvalue_Stage1<=0.05&met1$pvalue_stage2<=0.05)&(met1$sep_FvKO_pval<=0.05&met1$sep_MvKO_pval<=0.05) & ((met1$sep_FvKO_estimate>0 & met1$sep_MvKO_estimate>0)  | (met1$sep_FvKO_estimate<0 & met1$sep_MvKO_estimate<0))] = "Different size"
  
  
  met1$est_dirct[(met1$pvalue_Stage1<=0.05&met1$pvalue_stage2<=0.05)&((met1$sep_FvKO_pval<=0.05&met1$sep_MvKO_pval>0.05)|(met1$sep_FvKO_pval>0.05&met1$sep_MvKO_pval<=0.05))] = "One sex"
  
  met1$est_dirct[(met1$pvalue_Stage1<=0.05&met1$pvalue_stage2<=0.05)&(met1$sep_FvKO_pval>0.05&met1$sep_MvKO_pval>0.05 & met1$sep_allKO_pval > 0.05 )] = "Cannot Classify"
  
  # met1$est_dirct[(met1$pvalue_Stage1<=0.05&met1$pvalue_stage2<=0.05)&(met1$sep_FvKO_pval>0.05&met1$sep_MvKO_pval>0.05 & met1$sep_allKO_pval <= 0.05 )] = "Interaction_lack power"
  # met1$est_dirct[(met1$pvalue_Stage1<=0.05&met1$pvalue_stage2<=0.05)&(met1$sep_FvKO_pval>0.05&met1$sep_MvKO_pval>0.05 & met1$sep_allKO_pval <= 0.05 )] = "Genotype effect with no sex effect(real)"
  
  met1$est_dirct[(met1$pvalue_Stage1<=0.05&met1$pvalue_stage2>0.05)&(met1$sep_FvKO_pval <= 0.05 | met1$sep_MvKO_pval <= 0.05 | met1$sep_allKO_pval <= 0.05)] = "Genotype effect with no sex effect(real)"
  met1$est_dirct[(met1$pvalue_Stage1<=0.05&met1$pvalue_stage2>0.05)&(met1$sep_FvKO_pval > 0.05 & met1$sep_MvKO_pval > 0.05 & met1$sep_allKO_pval > 0.05)] = "Genotype effect with no sex effect(cannot classify)"
  
  
  met1$est_dirct[(met1$pvalue_Stage1>0.05 & met1$pvalue_stage2<=0.05)] = "interaction effect with no genotype effect"
  
  
  a<-sum(met1$est_dirct %in% "Different direction")
  b<-sum(met1$est_dirct %in% "Different size")
  c<-sum(met1$est_dirct %in%  "One sex")
  d<-sum(met1$est_dirct == "Cannot Classify")
  e<-sum(met1$est_dirct == "Genotype effect with no sex effect(real)")
  e1<-sum(met1$est_dirct == "Genotype effect with no sex effect(cannot classify)")
  f<-sum(met1$est_dirct == "interaction effect with no genotype effect")
  n<-sum(met1$est_dirct == "Not significant")
  
  
  PieChart_Diff_Dirct<-rbind(PieChart_Diff_Dirct,c(g,a,b,c,d,e,e1,f,n,t))
  colnames(PieChart_Diff_Dirct)<-coltitle 
}
write.csv(PieChart_Diff_Dirct,file = "5 Graph Piechart_raw p5%_PHENOtype_genoSD_0412.csv")

# ----------------------------------------------------------------------------------------------------------------------------
