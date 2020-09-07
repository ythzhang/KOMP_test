# =============================================================================================================================
# 5 for each genotype, check which metabolite are associated with gene * gender, gene, gender, or not significant at all.
# --------------------------------------------------------------------------------------------------------
# https://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0052410.s004
# https://www.researchgate.net/profile/Stano_Pekar/publication/304001371_Marginal_Models_Via_GLS_A_Convenient_Yet_Neglected_Tool_for_the_Analysis_of_Correlated_Data_in_the_Behavioural_Sciences/links/5779816e08aeb9427e2c0017/Marginal-Models-Via-GLS-A-Convenient-Yet-Neglected-Tool-for-the-Analysis-of-Correlated-Data-in-the-Behavioural-Sciences.pdf
# https://www.zoology.ubc.ca/~schluter/R/fit-model/
####### https://www.youtube.com/watch?v=vN5cNN2-HWE

#load package nlme
library(nlme)
library(RNOmni)
# library(car)
# Function 1 function to build model for testing of fixed effects
# Assumption: No batch effect and weight is not a independet variance 
# Goal: 

model_forFIXEDtest<-function(dataset, Intensity){
  
  model.formula <- as.formula(paste(Intensity, "~", paste("genotype", "sex", "genotype*sex", sep= "+")))
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


final_genotype_model<-function(dataset, Intensity){
  model_afterFIXED=model_forFIXEDtest(dataset, Intensity)
  anova_results = anova(model_afterFIXED, type="marginal")$"p-value" < 0.05
  keepSex = anova_results[3]
  keepInteraction = anova_results[4]
  
  if(keepSex && keepInteraction){
    return(model.formula <- as.formula(paste(Intensity, "~", paste("genotype", "sex", "genotype*sex", sep= "+"))))
  }else if(keepSex &&  !keepInteraction){
    return(model.formula <- as.formula(paste(Intensity, "~", paste("genotype", "sex", sep= "+"))))
  }else if(!keepSex &&  !keepInteraction){
    return(model.formula <- as.formula(paste(Intensity, "~", "genotype")))
  }else if(!keepSex  && keepInteraction){
    return(model.formula <- as.formula(paste(Intensity, "~", paste("genotype", "sex", "genotype*sex", sep= "+"))))
  }
}


#---------------------------------------------------------------------------------------------------------------------------
#Function 3: testing the fixed effects and building final Genotype effect null model
# Goal:  to test fixed effects of the model and based on the output build the final null model formula for later testing - as a null model it automatically excludes genotype and interaction term. 
# The anova function tests the fixed effects associated by treatment with a null hypothesis that the regression coefficients are equal to zero  and an alternative hypothesis that the regression coefficient are not equal to zero.
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


testing_Interaction_effect<-function(dataset, Intensity){
  #call functions to determine model.formula
  formula_interaction_null = null_model_Interaction(dataset,Intensity)
  model_interaction_full=gls(Intensity~genotype+sex+genotype*sex,dataset, method='ML',na.action="na.omit")
  model_null_interaction =gls(formula_interaction_null, dataset,method='ML', na.action="na.omit")
  return(pvalue_Stage2=(anova(model_interaction_full, model_null_interaction)$`p-value`[2]))
}


# ------------------------------------------------------------------------------------------------
#Function 4-d: testing the sex effect
null_model_sex<-function(dataset, Intensity){ 
  model_afterFIXED=model_forFIXEDtest(dataset, Intensity)
  anova_results = anova(model_afterFIXED, type="marginal")$"p-value" < 0.05
  keepSex = anova_results[3]
  keepGenotype = anova_results[2]     
  
  if(!keepSex && !keepGenotype){
    return(model.formula <- as.formula(paste(Intensity, "~", "1")))
  }else{
    return(model.formula <- as.formula(paste(Intensity, "~", "genotype")))
  } 
}

# ------------------------------------------------------------------------------------------------

testing_sex_effect<-function(dataset, Intensity){
  #call functions to determine model.formula
  formula_sex_null = null_model_sex(dataset,"Intensity")
  model_sex_full=gls(Intensity~genotype+sex,dataset, method='ML',na.action="na.omit")
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


finalmodel<-function(dataset, Intensity){
  #call functions to determine model.formula
  model_formula_null = null_model_genotype(dataset,Intensity)
  model_genotype = final_genotype_model(dataset, Intensity)
  
  model_genotype=gls(model_genotype,dataset, na.action="na.exclude")
}


# ---------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------


library(wcmc)
# setwd("C:\\Users\\yzhan\\Documents\\Fiehn\\3_KOMP\\7_Sexual D_20190605\\Jennly Zhang KOMP SD_ver05\\Input")
All_Metabolomics = wcmc::read_data("C:\\Users\\yzhan\\Documents\\Fiehn\\3_KOMP\\7_Sexual D_20190605\\Jennly Zhang KOMP SD_ver10\\Input\\ONLINE_data\\All Metabolomics_online.xlsx")
# All_Metabolomics = wcmc::read_data("All Metabolomics_norm.xlsx")
# All_Metabolomics = read.csv("C:\\Users\\yzhan\\Documents\\Fiehn\\3_KOMP\\7_Sexual D_20190605\\Jennly-Zhang-MetaboliteSD_ver03-master\\Jennly-Zhang-MetaboliteSD_ver03-master\\All_Metabolomics.csv",header=T)

All_Metabolomics_p = All_Metabolomics$p ##(first row plus label row)
All_Metabolomics_f = All_Metabolomics$f  ##(first column plus label column)
All_Metabolomics_e = All_Metabolomics$e_matrix   ###(just matrix values)


#All_Metabolomics_p = All_Metabolomics[1:3,]
#All_Metabolomics_f = All_Metabolomics[,1:2]
unique_genes = unique(All_Metabolomics_p$Genotype)

All_Metabolomics_e_no_mising = All_Metabolomics_e

All_Metabolomics_e[!is.finite(All_Metabolomics_e)] = NA
All_Metabolomics_e[All_Metabolomics_e==0]<-NA
All_Metabolomics_e[All_Metabolomics_e=="NA"]<-NA
All_Metabolomics_e[All_Metabolomics_e=="NaN"]<-NA
All_Metabolomics_e[All_Metabolomics_e=="Inf"]<-NA
sum(All_Metabolomics_e=="Inf")


# All_Metabolomics_e[All_Metabolomics_e==0] = NA
#All_Metabolomics_e[!is.finite(All_Metabolomics_e)] =0

# --------------------------------------------------------------------------------------------------------------------------------------------



for(g in 1:length(unique_genes)){
  
  current_gene = unique_genes[g]
  
  current_label = All_Metabolomics_p$label[All_Metabolomics_p$Genotype %in% c(current_gene)]
  current_e = All_Metabolomics_e[, current_label]
  
  current_label_F = All_Metabolomics_p$label[All_Metabolomics_p$Genotype %in% c(current_gene)&All_Metabolomics_p$Gender%in%'Female']
  current_e_F = All_Metabolomics_e[, current_label_F]
  
  current_label_M = All_Metabolomics_p$label[All_Metabolomics_p$Genotype %in% c(current_gene)&All_Metabolomics_p$Gender%in%'Male']
  current_e_M =All_Metabolomics_e[,current_label_M]
  
  if (current_gene=='null'){
    for(i in 1:nrow(current_e)){
      #Replace female NA
      if((sum(is.na(current_e_F[i,]))>=14) & (sum(is.na(current_e_M[i,]))>=14)){
        current_e_F[i,] = NA
        current_e_M[i,] = NA
      }else if((sum(is.na(current_e_F[i,]))>=14) & (sum(is.na(current_e_M[i,]))< 14)){
        current_e_F[i,is.na(current_e_F[i,])] = runif(sum(is.na(current_e_F[i,])),min=0.49,max=0.51) * min(current_e[i,!is.na(current_e[i,])])
        current_e_M[i,] = current_e_M[i,]
        
      }else if((sum(is.na(current_e_F[i,])) < 14) & (sum(is.na(current_e_M[i,])) >= 14)){
        current_e_F[i,] = current_e_F[i,]
        current_e_M[i,is.na(current_e_M[i,])] = runif(sum(is.na(current_e_M[i,])),min=0.49,max=0.51) * min(current_e[i,!is.na(current_e[i,])])
      }else{
        current_e_F[i,] = current_e_F[i,]
        current_e_M[i,] = current_e_M[i,]
      }
    }
  }else{
    
    for(i in 1:nrow(current_e)){
      if(sum(is.na(current_e[i,]))>=4){
        current_e_F[i,] = NA
        current_e_M[i,] = NA
      }else if(sum(is.na(current_e_F[i,]))==3 & sum(is.na(current_e_M[i,]))==0){
        current_e_F[i,is.na(current_e_F[i,])] = runif(sum(is.na(current_e_F[i,])),min=0.49,max=0.51) * min(current_e[i,!is.na(current_e[i,])])
        current_e_M[i,] = current_e_M[i,]
      }else if(sum(is.na(current_e_F[i,])) ==0 & sum(is.na(current_e_M[i,]))==3){
        current_e_F[i,] = current_e_F[i,]
        current_e_M[i,is.na(current_e_M[i,])] = runif(sum(is.na(current_e_M[i,])),min=0.49,max=0.51) * min(current_e[i,!is.na(current_e[i,])])
      }else if(sum(is.na(current_e_F[i,]))==2 & sum(is.na(current_e_M[i,]))==1){
        current_e_F[i,is.na(current_e_F[i,])] = runif(sum(is.na(current_e_F[i,])),min=0.49,max=0.51) * min(current_e[i,!is.na(current_e[i,])])
        current_e_M[i,is.na(current_e_M[i,])] = runif(sum(is.na(current_e_M[i,])),min=0.98,max=1.02) * min(current_e[i,!is.na(current_e[i,])])
      }else if(sum(is.na(current_e_F[i,]))==1 & sum(is.na(current_e_M[i,]))==2){
        current_e_F[i,is.na(current_e_F[i,])] = runif(sum(is.na(current_e_F[i,])),min=0.98,max=1.02) * min(current_e[i,!is.na(current_e[i,])])
        current_e_M[i,is.na(current_e_M[i,])] = runif(sum(is.na(current_e_M[i,])),min=0.49,max=0.51) * min(current_e[i,!is.na(current_e[i,])])
      }else if(sum(is.na(current_e[i,]))==2){
        if(sum(is.na(current_e_F[i,]))==2 & sum(is.na(current_e_M[i,]))==0){
          current_e_F[i,is.na(current_e_F[i,])] = runif(sum(is.na(current_e_F[i,])),min=0.98,max=1.02) * min(current_e[i,!is.na(current_e[i,])])
          current_e_M[i,] = current_e_M[i,]
        }else if(sum(is.na(current_e_F[i,]))==0 & sum(is.na(current_e_M[i,]))==2){
          current_e_F[i,] = current_e_F[i,]
          current_e_M[i,is.na(current_e_M[i,])] = runif(sum(is.na(current_e_M[i,])),min=0.98,max=1.02) * min(current_e[i,!is.na(current_e[i,])])
        }else if(sum(is.na(current_e_F[i,]))==1 & sum(is.na(current_e_M[i,]))==1){
          current_e_F[i,is.na(current_e_F[i,])] = runif(sum(is.na(current_e_F[i,])),min=0.98,max=1.02) * min(current_e[i,!is.na(current_e[i,])])
          current_e_M[i,is.na(current_e_M[i,])] = runif(sum(is.na(current_e_M[i,])),min=0.98,max=1.02) * min(current_e[i,!is.na(current_e[i,])])
        }
      }else if(sum(is.na(current_e[i,]))<=1){
        current_e_F[i,is.na(current_e_F[i,])] = runif(sum(is.na(current_e_F[i,])),min=0.98,max=1.02) * min(current_e[i,!is.na(current_e[i,])])
        current_e_M[i,is.na(current_e_M[i,])] = runif(sum(is.na(current_e_M[i,])),min=0.98,max=1.02) * min(current_e[i,!is.na(current_e[i,])])
      }
    }
  }
  All_Metabolomics_e_no_mising[,current_label_F] = current_e_F
  All_Metabolomics_e_no_mising[,current_label_M] = current_e_M
}



# =================================================    
#i = 759

for(i in 1: nrow(All_Metabolomics_e_no_mising)){
  
  All_Metabolomics_f$CompoundName[i]

# ---------------------------------------------------------
  # FEMALE
  
  outliers_F_raw_null <- boxplot.stats(All_Metabolomics_e_no_mising[i, All_Metabolomics_p$Gender %in% 'Female' & All_Metabolomics_p$Genotype %in% "null"])$out
  
  outliers_F_label_null <- All_Metabolomics_p$label[which(All_Metabolomics_p$Gender == "Female" & All_Metabolomics_p$Genotype %in% "null" & All_Metabolomics_e_no_mising[i,] %in% outliers_F_raw_null)]
  outliers_F_genotype_null <- All_Metabolomics_p$Genotype[which(All_Metabolomics_p$Gender == "Female" & All_Metabolomics_p$Genotype %in% "null" & All_Metabolomics_e_no_mising[i,] %in% outliers_F_raw_null)]
  
  if (length(outliers_F_raw_null) == 0){
    
    All_Metabolomics_e_no_mising[i,All_Metabolomics_p$Gender %in% 'Female'  & All_Metabolomics_p$Genotype %in% "null"] <- All_Metabolomics_e_no_mising[i,All_Metabolomics_p$Gender %in% 'Female' & All_Metabolomics_p$Genotype %in% "null"]
    
  }else{
    if(length(outliers_F_raw_null) > (sum(!is.na(All_Metabolomics_e_no_mising[i, All_Metabolomics_p$Gender %in% 'Female' & All_Metabolomics_p$Genotype %in% "null"])) /6)){
      
      outliers_F_raw_null <- NA
      outliers_F_genotype_null <- NA
      
    }else{
      
      outliers_F_raw_null <- outliers_F_raw_null
      outliers_F_label_null <- outliers_F_label_null
      outliers_F_genotype_null <- outliers_F_genotype_null
      
    }
    All_Metabolomics_e_no_mising[i, All_Metabolomics_p$label %in% outliers_F_label_null] <- NA
  }

  
  outliers_F_raw_KO <- boxplot.stats(All_Metabolomics_e_no_mising[i, All_Metabolomics_p$Gender %in% 'Female'])$out
  
  outliers_F_label_KO <- All_Metabolomics_p$label[which(All_Metabolomics_p$Gender == "Female" & All_Metabolomics_e_no_mising[i,] %in% outliers_F_raw_KO)]
  
  outliers_F_genotype_KO <- All_Metabolomics_p$Genotype[which(All_Metabolomics_p$label %in% outliers_F_label_KO)]
  
  ref_F_gene <- unique(outliers_F_genotype_KO)
  
  
  if(length(outliers_F_raw_KO)==0){
    All_Metabolomics_e_no_mising[i,All_Metabolomics_p$Gender %in% 'Female'] <- All_Metabolomics_e_no_mising[i,All_Metabolomics_p$Gender %in% 'Female']
    
  }else{
    
    for (l in 1: length(ref_F_gene)){
      
      if (ref_F_gene[l] == "null"){
        
        
        null_kp_F_lable <- which(outliers_F_genotype_KO == "null")
        outliers_F_raw_KO <- outliers_F_raw_KO[-null_kp_F_lable]
        outliers_F_genotype_KO <- outliers_F_genotype_KO[-null_kp_F_lable]
        outliers_F_label_KO <- outliers_F_label_KO[-null_kp_F_lable]
        
        # outliers_F_raw <- outliers_F_raw[-which(outliers_F_genotype== "null")]
        # outliers_F_genotype <- outliers_F_genotype[-which(outliers_F_genotype== "null")]
        # outliers_F_label <- outliers_F_label[-which(outliers_F_genotype == "null")]
        
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
    # All_Metabolomics_e_no_mising[i, which((All_Metabolomics_e_no_mising[i,] %in% outliers_F_raw) & (All_Metabolomics_p$Gender %in% 'Female')) ] <- NA
  }
  
  
  # ---------------------------------------------------------------------------------------------------
  # MALE
  
  outliers_M_raw_null <- boxplot.stats(All_Metabolomics_e_no_mising[i, All_Metabolomics_p$Gender %in% 'Male' & All_Metabolomics_p$Genotype %in% "null"])$out
  
  outliers_M_label_null <- All_Metabolomics_p$label[which(All_Metabolomics_p$Gender == "Male" & All_Metabolomics_p$Genotype %in% "null" & All_Metabolomics_e_no_mising[i,] %in% outliers_M_raw_null)]
  outliers_M_genotype_null <- All_Metabolomics_p$Genotype[which(All_Metabolomics_p$Gender == "Male" & All_Metabolomics_p$Genotype %in% "null" & All_Metabolomics_e_no_mising[i,] %in% outliers_M_raw_null)]
  
  if (length(outliers_M_raw_null) == 0){
    
    All_Metabolomics_e_no_mising[i,All_Metabolomics_p$Gender %in% 'Male' & All_Metabolomics_p$Genotype %in% "null"] <- All_Metabolomics_e_no_mising[i,All_Metabolomics_p$Gender %in% 'Male' & All_Metabolomics_p$Genotype %in% "null"]
    
  }else{
    if(length(outliers_M_raw_null) > (sum(!is.na(All_Metabolomics_e_no_mising[i, All_Metabolomics_p$Gender %in% 'Male' & All_Metabolomics_p$Genotype %in% "null"])) /6)){
      
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
  
  outliers_M_raw_KO <- boxplot.stats(All_Metabolomics_e_no_mising[i, All_Metabolomics_p$Gender %in% 'Male'])$out
  
  outliers_M_label_KO <- All_Metabolomics_p$label[which(All_Metabolomics_p$Gender == "Male" & All_Metabolomics_e_no_mising[i,] %in% outliers_M_raw_KO)]
  
  outliers_M_genotype_KO <- All_Metabolomics_p$Genotype[which(All_Metabolomics_p$label %in% outliers_M_label_KO)]
  
  ref_M_gene <- unique(outliers_M_genotype_KO)
  
  
  if(length(outliers_M_raw_KO)==0){
    All_Metabolomics_e_no_mising[i, All_Metabolomics_p$Gender %in% 'Male'] <- All_Metabolomics_e_no_mising[i, All_Metabolomics_p$Gender %in% 'Male']
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
    # All_Metabolomics_e_no_mising[i, which((All_Metabolomics_e_no_mising[i,] %in% outliers_M_raw) & (All_Metabolomics_p$Gender %in% 'Male')) ] <- NA 
  }
  
}


# -------------------------------------------------------------------------------


null_label = All_Metabolomics_p$label[All_Metabolomics_p$Genotype %in% c('null')]
null_e = All_Metabolomics_e_no_mising[, null_label]


null_label = All_Metabolomics_p$label[All_Metabolomics_p$Genotype %in% c('null') & All_Metabolomics_p$Gender %in% "Female"]
null_e_female = All_Metabolomics_e_no_mising[, null_label]
null_label = All_Metabolomics_p$label[All_Metabolomics_p$Genotype %in% c('null') & All_Metabolomics_p$Gender %in% "Male"]
null_e_male = All_Metabolomics_e_no_mising[, null_label]

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

# which(All_Metabolomics_f$label %in% "ONLI_1095")
# a<-All_Metabolomics_e_no_mising[759,]
# b <- All_Metabolomics_e[759,]



for(g in 2:length(unique_genes)){
  current_gene = unique_genes[g]
  print(g)
  current_label = All_Metabolomics_p$label[All_Metabolomics_p$Genotype %in% c(current_gene)]
  
  current_e = All_Metabolomics_e_no_mising[, current_label]
  
  current_p = All_Metabolomics_p[All_Metabolomics_p$label %in% current_label,]
  null_p = All_Metabolomics_p[All_Metabolomics_p$Genotype %in% 'null',]
  
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
  
  
  for(i in 1:nrow(All_Metabolomics_e_no_mising)){
    
    if((sum(is.na(current_e[i,])) >= (ncol(current_e)-2) )|(sum(is.na(null_e[i,])) >= (ncol(null_e)-6) )){
      # sex_index[[current_gene]][i]=NA
      
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
      # ==================================================================================
      data_raw<- data.table(y =c(current_e[i,], null_e[i,]), group = rep(c("Gene","Control"), c(ncol(current_e), ncol(null_e))), sex=c(current_p$Gender,null_p$Gender))
      colnames(data_raw) = c("Intensity_raw","genotype",'sex')
      data_raw<-data_raw[!is.na(data_raw$Intensity_raw), ]
      
      
      All_Metabolomics_f$CompoundName[i]
      
      Normality_before[[current_gene]][i]<-shapiro.test(data_raw$Intensity_raw)$p.value
      
      data_raw$Intensity<-rankNorm(data_raw$Intensity_raw)
      Normality_after[[current_gene]][i]<-shapiro.test(data_raw$Intensity)$p.value
      dataset<-data_raw[!is.na(data_raw$Intensity),]

      # -------------------------------------------------
      
      # outliers_nullF_raw<-boxplot.stats(data_raw$Intensity[data_raw$sex%in%'Female'&data_raw$genotype%in%'Control'])$out
      # if(length(outliers_nullF_raw)==0){
      #   dataset2_raw<-data_raw
      # }else{
      #   dataset2_raw<-data_raw[-which((data_raw$Intensity%in% outliers_nullF_raw)&(data_raw$sex%in%'Female')&(data_raw$genotype%in%'Control')),]
      # }
      # 
      # outliers_nullM_raw<-boxplot.stats(dataset2_raw$Intensity[dataset2_raw$sex%in%'Male'&dataset2_raw$genotype%in%'Control'])$out
      # if(length(outliers_nullM_raw)==0){
      #   dataset3_raw<-dataset2_raw
      # }else{
      #   dataset3_raw<-dataset2_raw[-which((dataset2_raw$Intensity%in% outliers_nullM_raw)&(dataset2_raw$sex%in%'Male')&(dataset2_raw$genotype%in%'Control')),]
      # }

      # dataset<-dataset3_raw[!is.na(dataset3_raw$Intensity),]
      
      # # -------------------------------------------------
      # # Outliers for genotype
      # outliers_genoF_raw<-boxplot.stats(dataset3_raw$Intensity[dataset3_raw$sex%in%'Female'&dataset3_raw$genotype%in%'Gene'])$out
      # 
      # 
      # if(length(outliers_genoF_raw)==0){
      #   dataset4_raw<-dataset3_raw
      # }else{
      #   dataset4_raw<-dataset3_raw[-which((dataset3_raw$Intensity%in% outliers_genoF_raw)&(dataset3_raw$sex%in%'Female')&(dataset3_raw$genotype%in%'Gene')),]
      # }
      # 
      # outliers_genoM_raw<-boxplot.stats(dataset4_raw$Intensity[dataset4_raw$sex%in%'Male'&dataset4_raw$genotype%in%'Gene'])$out
      # if(length(outliers_genoM_raw)==0){
      #   dataset_raw<-dataset4_raw
      # }else{
      #   dataset_raw<-dataset4_raw[-which((dataset4_raw$Intensity%in% outliers_genoM_raw)&(dataset4_raw$sex%in%'Male')&(dataset4_raw$genotype%in%'Gene')),]
      # }
      
      # dataset<-dataset_raw[!is.na(dataset_raw$Intensity),]
      
      
      # =========================================================================
      
      model_afterFIXED=tryCatch({
        model_forFIXEDtest(dataset, "Intensity")
      }, error = function(er){
        NA
      })
      if(length(model_afterFIXED)<=1){
        # sex_index[[current_gene]][i]=NA
        
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
        
        
        
        foldchange_all[[current_gene]][i]=mean(dataset$Intensity[dataset$genotype%in%"Gene"],na.rm = T)/mean(dataset$Intensity[dataset$genotype%in%"Control"],na.rm = T) 
        foldchange_FvKO[[current_gene]][i]= mean(dataset$Intensity[dataset$genotype%in%"Gene"&dataset$sex%in%"Female"],na.rm = T)/ mean(dataset$Intensity[dataset$genotype%in%"Control"&dataset$sex%in%"Female"],na.rm = T)
        foldchange_MvKO[[current_gene]][i]= mean(dataset$Intensity[dataset$genotype%in%"Gene"&dataset$sex%in%"Male"],na.rm = T)/ mean(dataset$Intensity[dataset$genotype%in%"Control"&dataset$sex%in%"Male"],na.rm = T)
        
        
        
      }else{
        
        model_formula_null=null_model_genotype(dataset,"Intensity")
        
        model_formula_genotype=final_genotype_model(dataset, "Intensity")
        
        geno_effect= testing_genotype_effect(dataset,"Intensity")
        
        #pvalue_genotype effect
        formula_interaction_null=null_model_Interaction(dataset,"Intensity")
        Interct_effect= testing_Interaction_effect(dataset,"Intensity")
        
        #pvalue_genotype effect
        formula_sex_null = null_model_sex(dataset,"Intensity")
        sex_effect= testing_sex_effect(dataset,"Intensity")
        
        
        #pvalue_SD
        result<-gls(Intensity~genotype+sex+genotype*sex,data = dataset,na.action = 'na.exclude')
        summary<-nlme:::summary.gls(result)$tTable 
        # summary = summary(lm(Intensity~genotype+sex+genotype*sex,data = dataset,na.action = 'na.omit'))$coef
        # cs <- as.data.frame(summary$tTable)
        
        # if(nrow(summary) == 4){
        
        
        Normality_before[[current_gene]][i]<-shapiro.test(data_raw$Intensity_raw)$p.value
        # data_raw$Intensity<-rankNorm(data_raw$Intensity_raw)
        Normality_after[[current_gene]][i]<-shapiro.test(data_raw$Intensity)$p.value
        Normality_F[[current_gene]][i]<-shapiro.test(data_raw$Intensity[data_raw$sex%in%"Female"])$p.value
        Normality_M[[current_gene]][i]<-shapiro.test(data_raw$Intensity[data_raw$sex%in%"Male"])$p.value
        
        numb_null[[current_gene]][i]= sum(dataset$genotype%in%"Control")
        
        numb_gene[[current_gene]][i]= sum(dataset$genotype%in%"Gene")
        
        pvalue_Stage1[[current_gene]][i]=geno_effect
        pvalue_Stage1.5[[current_gene]][i]=sex_effect
        pvalue_Stage2[[current_gene]][i]=Interct_effect
        
        # subset data for individual gls model in female and male    
        
        ctrl_data = dataset[dataset$genotype%in%"Control"]
        
        ctrl_test<- tryCatch({
          gls(Intensity~sex,data = ctrl_data,na.action = 'na.exclude')
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
          gls(Intensity~sex,data = KO_data,na.action = 'na.exclude')
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
          gls(Intensity~genotype,data = dataset,na.action = 'na.exclude')
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
          gls(Intensity~genotype,data = sub_F,na.action = 'na.exclude')
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
          gls(Intensity~genotype,data = sub_M,na.action = 'na.exclude')
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
        
        
        
        foldchange_all[[current_gene]][i]=mean(dataset$Intensity_raw[dataset$genotype%in%"Gene"],na.rm = T)/mean(dataset$Intensity_raw[dataset$genotype%in%"Control"],na.rm = T) 
        foldchange_FvKO[[current_gene]][i]= mean(dataset$Intensity_raw[dataset$genotype%in%"Gene"&dataset$sex%in%"Female"],na.rm = T)/ mean(dataset$Intensity_raw[dataset$genotype%in%"Control"&dataset$sex%in%"Female"],na.rm = T)
        foldchange_MvKO[[current_gene]][i]= mean(dataset$Intensity_raw[dataset$genotype%in%"Gene"&dataset$sex%in%"Male"],na.rm = T)/ mean(dataset$Intensity_raw[dataset$genotype%in%"Control"&dataset$sex%in%"Male"],na.rm = T)
        
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
sapply(pvalue_Stage2,function(x){sum(x<0.05,na.rm = TRUE)})
sapply(fdr_Stage1,function(x){sum(x<0.05,na.rm = TRUE)})
sapply(fdr_Stage1.5,function(x){sum(x<0.05,na.rm = TRUE)})
sapply(fdr_Stage2,function(x){sum(x<0.05,na.rm = TRUE)})

 
for(i in 1:length(names(pvalue_Stage1))){
  fwrite(data.table(label = All_Metabolomics_f$label, assay = All_Metabolomics_f$Assay, Metabolite=All_Metabolomics_f$CompoundName, Normality_before = Normality_before[[i]], Normality_after = Normality_after[[i]], Normality_F =Normality_F[[i]],  Normality_M =Normality_M[[i]] , numb_null=numb_null[[i]],numb_gene=numb_gene[[i]], pvalue_Stage1=pvalue_Stage1[[i]],pvalue_stage1.5=pvalue_Stage1.5[[i]], pvalue_stage2=pvalue_Stage2[[i]], fdr_Stage1=fdr_Stage1[[i]], fdr_Stage1.5=fdr_Stage1.5[[i]], fdr_Stage2=fdr_Stage2[[i]], sep_allKO_pval = sep_allKO_pval[[i]], sep_allKO_estimate= sep_allKO_estimate[[i]], sep_FvKO_pval=sep_FvKO_pval[[i]], sep_FvKO_estimate= sep_FvKO_estimate[[i]], sep_MvKO_pval=sep_MvKO_pval[[i]], sep_MvKO_estimate= sep_MvKO_estimate[[i]],sep_FvKO_fdr=p.adjust(sep_FvKO_pval[[i]],'fdr'), sep_MvKO_fdr=p.adjust(sep_MvKO_pval[[i]],'fdr'), foldchange_all= foldchange_all[[i]], foldchange_FvKO = foldchange_FvKO[[i]], foldchange_MvKO = foldchange_MvKO[[i]], ctrl_sex_pval= ctrl_sex_pval[[i]], ctrl_sex_estimate = ctrl_sex_estimate[[i]], fdr_ctrl_sex = fdr_ctrl_sex[[i]], KO_sex_pval= KO_sex_pval[[i]], KO_sex_estimate = KO_sex_estimate[[i]], fdr_KO_sex = fdr_KO_sex[[i]] ), paste0("C:\\Users\\yzhan\\Documents\\Fiehn\\3_KOMP\\7_Sexual D_20190605\\Jennly Zhang KOMP SD_ver10\\Output\\20200528_METAbo_genoSD_koSEX_rm first\\5 TwoStage METAbol_koSEX (",names(pvalue_Stage1)[i],"),genoSD_0528.csv"))
}




# =======================================================================================================

#5d Figure 4 information extraction pie chart-metabolite~gene&gender&direction
# combine dataframe of female and male
# ==============================================================================================================================
# End of run using 
# RAW pvalue
# ==============================================================================================================================
# ==============================================================================================================================
# ==============================================================================================================================
# ==============================================================================================================================
# ==============================================================================================================================
# sTART of run using



library(wcmc)
setwd("C:\\Users\\yzhan\\Documents\\Fiehn\\3_KOMP\\7_Sexual D_20190605\\Jennly Zhang KOMP SD_ver10\\Output\\20200528_METAbo_genoSD_koSEX_rm first")
GeneList<-read.csv("C:\\Users\\yzhan\\Documents\\Fiehn\\3_KOMP\\7_Sexual D_20190605\\Jennly Zhang KOMP SD_ver10\\Input\\GeneList_Metabol.csv",header = TRUE)

files_metmerge<-list.files(pattern = "^5 TwoStage METAbol_koSEX (.*),genoSD_0528.csv",full.names = T)


# GeneList<-as.matrix(GeneList)
PieChart_Diff_Dirct<-data.frame(matrix(ncol = 16, nrow = 0))
coltitle<-c('Genotype','Different direction',"Different direction(maybe)",'Different size',"Different size(maybe)","One sex_Female only", "One sex_Male only","One sex_Female only (maybe)","One sex_Male only (maybe)", "Genotype/sex interaction(maybe)", "Genotype effect with no sex effect(real)", "Genotype effect with no sex effect(maybe)", "Genotype effect with no sex effect(cannot classify)","interaction effect with no genotype effect",'Not Significant',"Metabolites count")
colnames(PieChart_Diff_Dirct)<-coltitle


for (i in 1:length(files_metmerge)){
  met<-read.csv(files_metmerge[i])  #met<-read.csv(files_metmerge[1])
  g<-i
  print(i)
  met1<-met[!is.na(met$pvalue_Stage1),]
  met1$est_dirct<-"Not significant"
  t<-length(met1$pvalue_Stage1)
  # print(t)
  # for (k in 1:length(met1$pvalue_Stage1)){
  met1$est_dirct[(met1$pvalue_Stage1<=0.05&met1$pvalue_stage2<=0.05)&(met1$sep_FvKO_pval<=0.05&met1$sep_MvKO_pval<=0.05 )& ((met1$sep_FvKO_estimate>0& met1$sep_MvKO_estimate<0) | (met1$sep_FvKO_estimate<0 & met1$sep_MvKO_estimate>0))] = "Different direction"
  
  met1$est_dirct[(met1$pvalue_Stage1<=0.05&met1$pvalue_stage2<=0.05)&(met1$sep_FvKO_pval<=0.05&met1$sep_MvKO_pval<=0.05) & ((met1$sep_FvKO_estimate>0 & met1$sep_MvKO_estimate>0)  | (met1$sep_FvKO_estimate<0 & met1$sep_MvKO_estimate<0))] = "Different size"
  
  met1$est_dirct[(met1$pvalue_Stage1<=0.05&met1$pvalue_stage2<=0.05)&(met1$sep_FvKO_pval<=0.05&met1$sep_MvKO_pval>0.05)] = "One sex_Female only"
  met1$est_dirct[(met1$pvalue_Stage1<=0.05&met1$pvalue_stage2<=0.05)&(met1$sep_FvKO_pval>0.05&met1$sep_MvKO_pval<=0.05)] = "One sex_Male only"
  
  met1$est_dirct[(met1$pvalue_Stage1<=0.05&met1$pvalue_stage2<=0.05)&(met1$sep_FvKO_pval>0.05&met1$sep_MvKO_pval>0.05 & met1$sep_allKO_pval > 0.05 )] = "Genotype/sex interaction(maybe)"
  
  # met1$est_dirct[(met1$pvalue_Stage1<=0.05&met1$pvalue_stage2<=0.05)&(met1$sep_FvKO_pval>0.05&met1$sep_MvKO_pval>0.05 & met1$sep_allKO_pval <= 0.05 )] = "Genotype effect with no sex effect(real)"
  
  
  met1$est_dirct[(met1$pvalue_Stage1<=0.05&met1$pvalue_stage2>0.05)&(met1$sep_FvKO_pval<=0.05&met1$sep_MvKO_pval<=0.05 )& ((met1$sep_FvKO_estimate>0& met1$sep_MvKO_estimate<0) | (met1$sep_FvKO_estimate<0 & met1$sep_MvKO_estimate>0))] = "Different direction(maybe)"
  met1$est_dirct[(met1$pvalue_Stage1<=0.05&met1$pvalue_stage2>0.05)&(met1$sep_FvKO_pval<=0.05&met1$sep_MvKO_pval<=0.05) & ((met1$sep_FvKO_estimate>0 & met1$sep_MvKO_estimate>0)  | (met1$sep_FvKO_estimate<0 & met1$sep_MvKO_estimate<0))] = "Different size(maybe)"
  
  met1$est_dirct[(met1$pvalue_Stage1<=0.05&met1$pvalue_stage2>0.05)&(met1$sep_allKO_pval > 0.05 & (met1$sep_FvKO_pval<=0.05&met1$sep_MvKO_pval>0.05))] = "One sex_Female only (maybe)"
  met1$est_dirct[(met1$pvalue_Stage1<=0.05&met1$pvalue_stage2>0.05)&(met1$sep_allKO_pval > 0.05 & (met1$sep_FvKO_pval>0.05&met1$sep_MvKO_pval<=0.05))] = "One sex_Male only (maybe)"
  met1$est_dirct[(met1$pvalue_Stage1<=0.05&met1$pvalue_stage2>0.05)&( (met1$sep_FvKO_pval>0.05 | met1$sep_MvKO_pval>0.05) & met1$sep_allKO_pval <= 0.05 )] = "Genotype effect with no sex effect(maybe)"                
  met1$est_dirct[(met1$pvalue_Stage1<=0.05&met1$pvalue_stage2>0.05)&(met1$sep_FvKO_pval<=0.05&met1$sep_MvKO_pval<=0.05)] = "Genotype effect with no sex effect(real)"                 
  
  met1$est_dirct[(met1$pvalue_Stage1<=0.05&met1$pvalue_stage2>0.05)&(met1$sep_FvKO_pval > 0.05 & met1$sep_MvKO_pval > 0.05 & met1$sep_allKO_pval > 0.05)] = "Genotype effect with no sex effect(cannot classify)"
  
  
  met1$est_dirct[(met1$pvalue_Stage1 > 0.05 & met1$pvalue_stage2 <= 0.05)] = "interaction effect with no genotype effect"
  
  
  a<-sum(met1$est_dirct == "Different direction")
  a1<-sum(met1$est_dirct == "Different direction(maybe)")
  b <- sum(met1$est_dirct == "Different size")
  b1 <- sum(met1$est_dirct == "Different size(maybe)")
  c<-sum(met1$est_dirct == "One sex_Female only")
  c1 <- sum(met1$est_dirct == "One sex_Male only")
  c2 <- sum(met1$est_dirct == "One sex_Female only (maybe)")
  c3 <- sum(met1$est_dirct == "One sex_Male only (maybe)")
  d<-sum(met1$est_dirct=="Genotype/sex interaction(maybe)")
  e<-sum(met1$est_dirct=="Genotype effect with no sex effect(real)")
  e1<-sum(met1$est_dirct=="Genotype effect with no sex effect(maybe)")
  e2<-sum(met1$est_dirct=="Genotype effect with no sex effect(cannot classify)")
  f<-sum(met1$est_dirct=="interaction effect with no genotype effect")
  n<-sum(met1$est_dirct=="Not significant")
  
  PieChart_Diff_Dirct<-rbind(PieChart_Diff_Dirct,c(g,a,a1,b,b1,c,c1,c2,c3,d,e,e1,e2,f,n,t))
  colnames(PieChart_Diff_Dirct)<-coltitle 
}

write.csv(PieChart_Diff_Dirct,file = "5 Graph raw p5%_For Piechart_Metabolomics_more included_0528.csv")



# ----------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------
# ==============================================================================================================================
# ==============================================================================================================================
# sTART of run using


library(wcmc)
setwd("C:\\Users\\yzhan\\Documents\\Fiehn\\3_KOMP\\7_Sexual D_20190605\\Jennly Zhang KOMP SD_ver10\\Output\\20200412_METAbo_genoSD_koSEX_rm first")
GeneList<-read.csv("C:\\Users\\yzhan\\Documents\\Fiehn\\3_KOMP\\7_Sexual D_20190605\\Jennly Zhang KOMP SD_ver10\\Input\\GeneList_Metabol.csv",header = TRUE)

files_metmerge<-list.files(pattern = "^5 TwoStage METAbol_koSEX (.*),genoSD_0412.csv",full.names = T)


# GeneList<-as.matrix(GeneList)
PieChart_Diff_Dirct<-data.frame(matrix(ncol = 10, nrow = 0))
coltitle<-c('Genotype','Different direction','Different size','One sex', "Cannot Classify", "Genotype effect with no sex effect(real)", "Genotype effect with no sex effect(cannot classify)","interaction effect with no genotype effect",'Not Significant',"Metabolites count")
colnames(PieChart_Diff_Dirct)<-coltitle


for (i in 1:length(files_metmerge)){
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

write.csv(PieChart_Diff_Dirct,file = "5 Graph raw p5%_For Piechart_Metabolomics_less included_0412.csv")



# ----------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------
