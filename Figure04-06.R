# ---------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------
library(patchwork)
library(wcmc)
library(ggplot2)
library(dplyr)
library(RNOmni)
library(forcats)
setwd("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
All_Metabolomics = wcmc::read_data("Supplementary Data1.xlsx")
All_Metabolomics_p = All_Metabolomics$p ##(first row plus label row)
All_Metabolomics_f = All_Metabolomics$f  ##(first column plus label column)
All_Metabolomics_e = All_Metabolomics$e_matrix   ###(just matrix values)
GeneList <- read.csv(file="GeneList_Metabol.csv",header = TRUE)
files_met <- list.files(pattern = "^5 TwoStage METAbol_koSEX (.*),genoSD.csv",full.names = T)
Speci_Metabol <- data.frame(matrix(ncol = 32, nrow = 0))
coltitle <- c("label", "assay",	"Metabolite",	"Normality_before",	"Normality_after",	"Normality_F",	"Normality_M",	"numb_null",	"numb_gene",	"pvalue_Stage1",	"pvalue_stage1.5",	"pvalue_stage2",	"fdr_Stage1",	"fdr_Stage1.5",	"fdr_Stage2",	"sep_allKO_pval",	"sep_allKO_estimate",	"sep_FvKO_pval",	"sep_FvKO_estimate",	"sep_MvKO_pval",	"sep_MvKO_estimate",	"sep_FvKO_fdr",	"sep_MvKO_fdr",	"foldchange_all",	"foldchange_FvKO",	"foldchange_MvKO",	"ctrl_sex_pval",	"ctrl_sex_estimate",	"fdr_ctrl_sex",	"KO_sex_pval",	"KO_sex_estimate",	"fdr_KO_sex")
colnames(Speci_Metabol) <- coltitle
GeneList <- (GeneList$Genotype)
Genelist <- as.character(GeneList)
for (g in 1:length(files_met)){
  gene <- Genelist[g]
  met <- read.csv(files_met[g]) 
  met1 <- na.omit(met)
  met1$Genotype <- gene
  Speci_Metabol <- rbind(Speci_Metabol,met1)
}
pvalFrame_e = All_Metabolomics_e
pvalFrame_p = All_Metabolomics_p
pvalFrame_f = All_Metabolomics_f
# pvalFrame_e = NA
ColorFrame_e = All_Metabolomics_e
ColorFrame_p = All_Metabolomics_p
ColorFrame_f = All_Metabolomics_f
# ONLI_0896,703
which(Speci_Metabol$label == "ONLI_0703")
which(All_Metabolomics_f$label == "ONLI_0703")
All_Metabolomics_f$CompoundName[7]
List <- unique(pvalFrame_p$Genotype)
for (i in 1 : nrow(pvalFrame_e)){
  data_A<-Speci_Metabol[Speci_Metabol$label %in% pvalFrame_f$label[i],]
  # data_A$label  length(data_A$label)
  pvalFrame_e[i, pvalFrame_p$Genotype %in% "null"] = "Not significant"
  ColorFrame_e[i, ColorFrame_p$Genotype %in% "null"] = "Not significant"
  for (g in 2 : length(List)){
    Gene=List[g]
    if(Gene %in% data_A$Genotype){
      pvalFrame_e[i, pvalFrame_p$Genotype %in% Gene & pvalFrame_p$Gender == "Female"] = data_A$sep_FvKO_pval[data_A$Genotype %in% Gene]
      pvalFrame_e[i, pvalFrame_p$Genotype %in% Gene & pvalFrame_p$Gender == "Male"] = data_A$sep_MvKO_pval[data_A$Genotype %in% Gene]
    }else{
      pvalFrame_e[i, pvalFrame_p$Genotype %in% Gene & pvalFrame_p$Gender == "Female"] = "Not significant"
      pvalFrame_e[i, pvalFrame_p$Genotype %in% Gene & pvalFrame_p$Gender == "Male"] = "Not significant"
    }
  }
  ColorFrame_e[i, pvalFrame_e[i,] <= 0.05] = "Significant"
  ColorFrame_e[i, pvalFrame_e[i,] > 0.05] = "Not significant"
}
# unique(ColorFrame_e[7,])
# ColorFrame_f$label[7]
Classi <- c()
for(i in 1 : nrow(ColorFrame_e)){
  u <- unique(ColorFrame_e[i,])
  Classi <- c(Classi,u)
}
unique(Classi)
# ------------------------------------------------------------------------------------------------------------------------------------------
#Starting boxplot
which(Speci_Metabol$label== "ONLI_0896")
which(All_Metabolomics_f$label == "ONLI_1095")
for(i in 1:nrow(All_Metabolomics_e_no_mising)){
  print(i)
  Metabolite_Label <- All_Metabolomics_f$label[i]
  Metabolite_Name <- All_Metabolomics_f$CompoundName[i]
  data_raw <- data.table(All_Metabolomics_e_no_mising[i, ],ColorFrame_e[i,], All_Metabolomics_p)
  colnames(data_raw)[1:2] <- c("MetaValues","Classify")
  unique(data_raw$Classify)
  data_raw <- data_raw[!is.na(data_raw$MetaValues),]
  data_F <- data_raw[data_raw$Gender %in% "Female"]
  data_M <- data_raw[data_raw$Gender %in% "Male"]
  # -------------------------------------------------
  # calculations
  data_M$mean_M[data_M$Genotype %in% "null"]<-mean(data_M$MetaValues[data_M$Genotype %in% "null"], na.rm=TRUE)
  data_M$foldchange[data_M$Genotype %in%  "null"]<-   data_M$MetaValues[data_M$Genotype %in% "null"] /unique(data_M$mean_M[data_M$Genotype %in% "null"])
  data_M$fc_mean[data_M$Genotype %in% "null"] <- mean(data_M$foldchange[data_M$Genotype %in% "null"], na.rm = T)
  data_M$fc_sd[data_M$Genotype %in% "null"] <- sd(data_M$foldchange[data_M$Genotype %in% "null"], na.rm = T)
  data_M$fc_z[data_M$Genotype %in% "null"] <- (data_M$foldchange[data_M$Genotype %in% "null"] - unique(data_M$fc_mean[data_M$Genotype %in% "null"]))/unique(data_M$fc_sd[data_M$Genotype %in% "null"])
  data_M$fc_z_mean[data_M$Genotype %in% "null"] <- mean(data_M$fc_z[data_M$Genotype %in% "null"], na.rm = T)
  data_M$fc_z_stdev[data_M$Genotype %in% "null"] <- sd(data_M$fc_z[data_M$Genotype %in% "null"], na.rm = T)
  # head(data_M)
  for (s in 2:length(unique(data_M$Genotype))){
    Gene_sd<-unique(data_M$Genotype)[s]
    data_M$mean_M[data_M$Genotype %in% Gene_sd] <- mean(data_M$MetaValues[data_M$Genotype %in% Gene_sd],na.rm=TRUE)
    data_M$foldchange[data_M$Genotype %in% Gene_sd]<- data_M$MetaValues[data_M$Genotype %in% Gene_sd]/unique(data_M$mean_M[data_M$Genotype %in% "null"])
    data_M$fc_mean[data_M$Genotype %in% Gene_sd] <- mean(data_M$foldchange[data_M$Genotype %in% Gene_sd], na.rm = T)
    data_M$fc_sd[data_M$Genotype %in% Gene_sd] <- sd(data_M$foldchange[data_M$Genotype%in% Gene_sd], na.rm = T)
    data_M$fc_z[data_M$Genotype %in% Gene_sd] <- (data_M$foldchange[data_M$Genotype %in% Gene_sd] - unique(data_M$fc_mean[data_M$Genotype %in% "null"])) /  unique(data_M$fc_sd[data_M$Genotype %in% "null"])
    data_M$fc_z_mean[data_M$Genotype %in% Gene_sd] <- mean(data_M$fc_z[data_M$Genotype %in% Gene_sd], na.rm = T)
    data_M$fc_z_stdev[data_M$Genotype %in% Gene_sd] <- sd(data_M$fc_z[data_M$Genotype %in% Gene_sd], na.rm = T)
  }
  data_F$mean_F[data_F$Genotype %in% "null"] <- mean(data_F$MetaValues[data_F$Genotype %in% "null"],na.rm=TRUE)
  data_F$foldchange[data_F$Genotype %in% "null"] <- data_F$MetaValues[data_F$Genotype %in% "null"] /unique(data_F$mean_F[data_F$Genotype %in% "null"])   
  data_F$fc_mean[data_F$Genotype %in% "null"] <- mean(data_F$foldchange[data_F$Genotype %in% "null"],na.rm=TRUE)
  data_F$fc_sd[data_F$Genotype %in% "null"] <- sd(data_F$foldchange[data_F$Genotype %in% "null"],na.rm=TRUE)
  data_F$fc_z[data_F$Genotype %in% "null"] <- (data_F$foldchange[data_F$Genotype %in% "null"] - unique(data_F$fc_mean[data_F$Genotype %in% "null"]))/ unique(data_F$fc_sd[data_F$Genotype %in% "null"])
  data_F$fc_z_mean[data_F$Genotype %in% "null"] <- mean(data_F$fc_z[data_F$Genotype %in% "null"],na.rm=TRUE)
  data_F$fc_z_stdev[data_F$Genotype %in% "null"] <- sd(data_F$fc_z[data_F$Genotype %in% "null"],na.rm=TRUE)
  for (t in 2:length(unique(data_F$Genotype))){
    Gene_sd <- unique(data_F$Genotype)[t]
    data_F$mean_F[data_F$Genotype %in% Gene_sd] <- mean(data_F$MetaValues[data_F$Genotype %in% Gene_sd],na.rm=TRUE)
    data_F$foldchange[data_F$Genotype %in% Gene_sd] <- data_F$MetaValues[data_F$Genotype %in% Gene_sd]/unique(data_F$mean_F[data_F$Genotype %in% "null"])   
    data_F$fc_mean[data_F$Genotype %in% Gene_sd] <- mean(data_F$foldchange[data_F$Genotype %in% Gene_sd],na.rm=TRUE)
    data_F$fc_sd[data_F$Genotype %in% Gene_sd] <- sd(data_F$foldchange[data_F$Genotype %in% Gene_sd],na.rm=TRUE)
    data_F$fc_z[data_F$Genotype %in% Gene_sd] <- (data_F$foldchange[data_F$Genotype %in% Gene_sd] - unique(data_F$fc_mean[data_F$Genotype %in% "null"])) /  unique(data_F$fc_sd[data_F$Genotype %in% "null"])
    data_F$fc_z_mean[data_F$Genotype %in% Gene_sd] <- mean(data_F$fc_z[data_F$Genotype %in% Gene_sd],na.rm=TRUE)
    data_F$fc_z_stdev[data_F$Genotype %in% Gene_sd] <- sd(data_F$fc_z[data_F$Genotype %in% Gene_sd],na.rm=TRUE)
  }
  # head(data_F)
  data_M_new <- data_M[!data_M$Genotype %in% "null"& !is.na(data_M$MetaValues),]
  data_F_new <- data_F[!data_F$Genotype %in%"null"& !is.na(data_F$MetaValues),]
  for(l in 1: length(unique(data_M_new$Genotype))){
    Gene_length <-unique(data_M_new$Genotype)[l]
    if(sum(data_M_new$Genotype%in% Gene_length) == sum(data_M_new$Genotype%in% Gene_length)){
      print(l)
      data_M_new <- data_M_new
      data_F_new <- data_F_new
    }else{
      data_M_new <- data_M_new[!data_M$Genotype %in% Gene_length,]
      data_F_new <- data_F_new[!data_F$Genotype %in% Gene_length,]
    }
  }
    data_M_new[data_M_new$Genotype == "Mfap4",]
    p_M <- data_M_new %>%
      ggplot(aes(x = fct_reorder(Genotype,desc(Order)), y = fc_z,  color = Classify)) +
      geom_errorbar(aes(ymax = fc_z_mean + fc_z_stdev,  ymin = fc_z_mean - fc_z_stdev),position = "dodge", width = 0.6, size= 1) +
      stat_summary(fun.y="mean" , geom = "point", position=position_dodge(width=0.01), size=2) +
      coord_flip() +
      scale_colour_manual(values =c("black", "red"), aesthetics = c("colour", "fill")) +
      theme(axis.text.x = element_text(angle=0,hjust=0.5)) + 
      geom_hline(yintercept = 0,color="dimgray") +
      scale_y_continuous(name = "Effect size (Male)", limits=c(-10,10),breaks = seq(-10,10, by = 5)) +
      theme_bw() +
      theme(legend.position =  "bottom",  
            legend.title =  element_text(size = 10,color = "black"),
            axis.text.x = element_text(colour = "black",hjust=0.5,vjust=0.5,margin = margin(t = 0, r = 100, b = 0, l = 100)),
            axis.text.y = element_text(colour = "black",hjust=0.5,vjust=0.5))+ #legend.position = c(0.88, 0.88),
      xlab("Genotype")
    p_F <- data_F_new %>%
      ggplot(aes(x= fct_reorder(Genotype,desc(Order)), y=fc_z, color = Classify)) +
      geom_errorbar(aes(ymax = fc_z_mean + fc_z_stdev,  ymin = fc_z_mean - fc_z_stdev),position = "dodge", width = 0.6, size= 1) +
      stat_summary(fun.y="mean", geom = "point", position=position_dodge(width=0.01), size=2)+
      coord_flip() +
      scale_colour_manual(values =c("black", "red"), aesthetics = c("colour", "fill")) +   #"darkgreen",
      theme(axis.text.x = element_text(angle=0,hjust=0.5)) +
      geom_hline(yintercept = 0,color="dimgray") +
      scale_y_continuous(name = "Effect size (Female)", limits=c(-10,10),breaks = seq(-10,10, by = 5)) +
      theme_bw() +
      theme(legend.position =  "bottom",  
            legend.title =  element_text(size = 10,color = "black"),
            axis.text.x = element_text(colour = "black",hjust=0.5,vjust=0.5,margin = margin(t = 0, r = 100, b = 0, l = 100)),
            axis.text.y = element_text(colour = "black",hjust=0.5,vjust=0.5))+ #legend.position = c(0.88, 0.88),
      xlab("Genotype")
     x11(width=10,height=10)
     p_F|p_M
     savePlot(file=paste("Boxplot_genoSD,",i,",", Metabolite_Name,".emf"))
     dev.off()
}
End of Figure 4
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
#