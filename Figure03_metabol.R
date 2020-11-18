##Figure 3a pie chart metabolomics: genotype-gender interaction
library(ggplot2)  # Data visualization
library(dplyr)    # Data manipulation
Mydata<-read.csv("Figure 4a_Metabol.csv")
df <- data.frame(Mydata)
df
# mycols <- c("#CD534CFF","#0073C2FF","#868686FF")
# colors()
library(colorspace)
hcl_palettes(plot = T)
col1 <- diverge_hcl(9,"Blue-Red2")
barplot(1:9, col = col1)
dd <- col1[7]
ds <- col1[9]
os <- col1[1]
cc <- col1[3]
gns <- col1[5]
sng <- col1[6]
col1 <- c(cc, dd, ds, os, gns, sng)
blank_theme <- theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
pie <- ggplot(df, aes(x = "", y = Percent, fill = Group)) + blank_theme +
  theme(axis.text.x=element_blank()) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = pi/2,direction = 1) +
  geom_text(aes(label = paste0(df$Percent, "%")), position = position_stack(vjust = 0.5), size=5) +
  scale_fill_manual(values = col1)
setwd("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
x11(width=8, height=8)
pie
savePlot("Figure 4a_Metabol Pie_genoSD.emf",type="emf")
dev.off()
# End of Figure 3a
# ----------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------
# Figure 3b Pie plot_Phenotype
library(ggplot2)  # Data visualization
library(dplyr)    # Data manipulation
Mydata<-read.csv("Figure 3b_Pheno_Pie.csv")
df <- data.frame(Mydata)
df
library(colorspace)
hcl_palettes(plot = T)
col1 <- diverge_hcl(10,'Blue-Red2')
barplot(1:10,col=col1)
dd <- col1[7]
ds <- col1[9]
fo <- col1[1]
mo <- col1[4]
cc <- col1[2]
gns <- col1[5]
col1 <- c(dd,ds,fo,cc,mo,gns)
blank_theme <- theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
pie <- ggplot(df, aes(x = "", y = Percent, fill = Group)) + blank_theme +
  theme(axis.text.x=element_blank())+
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = pi/2.6,direction = 1)+
  geom_text(aes(label = paste0(df$Percent, "%")), position = position_stack(vjust = 0.5), size = 5) +
  scale_fill_manual(values = col1)
setwd("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
x11(width=8,height=8)
pie
savePlot("Figure3b Pheno Pie_genoSD.emf",type="emf")
dev.off()
# End of Figure 3b
# ----------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------
# Figure 3c bar plot_Metabolomics sexual dimorphism for each gene knockout strain
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(ggthemes)
library(stringr)
display.brewer.all()
df<-read.csv("Figure 4b_Metabol_Bar.csv")
head(df)
Order_df <-
  df %>% 
  filter(Group == "Percent_Inter") %>%
  arrange((Percent_gene))
the_order <- Order_df$Genotype
head(Order_df)
# pal <- brewer.pal(5,"RdBu")
# plot(c(1,2,3,4,5), col = pal) 
pal<-c("antiquewhite", "sienna1", "gray70", "E69F00")
legend.pal <- pal
df$col <- rep(pal, each = 30)
interac_effect <-df[(length(df[,1])/4) + 1: (length(df[,1])/4),]
geno_effect <-df[1: (length(df[,1])/4),]
aNot_sig <- df[(length(df[,1])*0.5) + 1 : (length(df[,1])*0.25) ,]
SD_lipids <- df[(length(df[,1])*0.75) + 1 : length(df[,1]),][1:30,]
Levels_group <- c("geno_effect", "interact_effect", "SD_lipids", "Not_sig")
trans <- function(x){pmin(x,30) + 1*pmax(x-30,0)}
yticks <- c(0, 5, 10, 15, 20, 25, 30, 100)
#Transform the data onto the display scale
aNot_sig$Percent_gene_t <- trans(aNot_sig$Percent_gene)
geno_effect$Percent_gene_t <- trans(geno_effect$Percent_gene)
interac_effect$Percent_gene_t <- trans(interac_effect$Percent_gene) 
SD_lipids$Percent_gene_t <- trans(SD_lipids$Percent_gene) 
p<- 
  df %>% 
  ggplot() +
  geom_bar(data= aNot_sig, aes(x = Genotype, y=Percent_gene, fill="gray70"), color="black",position="stack", stat="identity") +
  geom_bar(data= geno_effect, aes(x = Genotype, y=Percent_gene, fill="aquamarine4"), color="black",position="stack", stat="identity") +
  geom_bar(data= interac_effect, aes(x = Genotype, y=Percent_gene, fill= "sienna"),color="black",position="stack", stat="identity") +
  geom_bar(data= SD_lipids, aes(x = Genotype, y=Percent_gene, fill= "#E69F00"),color="black",position="stack", stat="identity") +
  geom_hline(yintercept = 0, color =c("white")) +
  scale_fill_identity("Group", labels = Levels_group, breaks=legend.pal, guide="legend") + 
  # theme_fivethirtyeight() + 
  # coord_flip() +
  scale_x_discrete(limits = the_order) +
  labs(title="overlap", y="",x="") +
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(plot.title = element_text(size=14, hjust=0.5)) +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black"),
        axis.text.y = element_text(hjust=0.5,vjust = 0.5,colour = "black")) +
  # scale_y_continuous(breaks=trans(yticks), labels=yticks, limits=c(0,35),sec.axis = sec_axis(~.*1,breaks =trans(yticks)))
  scale_y_continuous(breaks=seq(0,100,5), limits=c(0,100),sec.axis = sec_axis(~.*1,breaks = seq(0,100,5)))
x11(width=14,height=10)
p
savePlot("Figure 4b_Metabol Bar.emf",type="emf")
dev.off()  
# End of Figure 3c
# -----------------------------------------------------------------------------------------------------------------------
# Figure 3d-f PLSDA for each KO
library(ggplot2)
library(devtools)
library(devEMF)
library(ropls)
library(tidyverse)
library(wcmc)
setwd("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
All_Metabolomics = wcmc::read_data("Dataset for Figure3d_Metabolomics.xlsx")
All_Metabolomics_p = All_Metabolomics$p ##(first row plus label row)
All_Metabolomics_f = All_Metabolomics$f  ##(first column plus label column)
All_Metabolomics_e = All_Metabolomics$e_matrix   ###(just matrix values)
unique_genes = unique(All_Metabolomics_p$Genotype)
for(g in 2:length(unique_genes)){
  current_gene = unique_genes[g]
  print(g)
  current_label = All_Metabolomics_p$label[All_Metabolomics_p$Genotype %in% c("null", current_gene)]
  current_e = All_Metabolomics_e_no_mising[, current_label]
  current_p = All_Metabolomics_p[All_Metabolomics_p$label %in% current_label,]
  current_p$Gender <-  as.factor(current_p$Gender)
  current_p$group <- paste0(current_p$Genotype, current_p$Gender)
  Metabol_genderFc <- current_p$group  
  # nrow(current_e_plsda)
  current_e_plsda <- current_e
  rownames(current_e_plsda) <- All_Metabolomics_f$label
  colnames(current_e_plsda) <- current_label
  row_index <- c()
  for(i in 1 : nrow(current_e_plsda)){
    if(sum(is.na(current_e_plsda[i,])) > 0.7* ncol(current_e_plsda) ){
      row_index <- c(row_index,i)
    }else{
      current_e_plsda[i,is.na(current_e_plsda[i,])] = runif(sum(is.na(current_e_plsda[i,])),min=0.49,max=0.51) * min(current_e_plsda[i,!is.na(current_e_plsda[i,])])
    }
  }
  current_e_plsda <- current_e_plsda[-row_index,]
  opls(t(current_e_plsda), Metabol_genderFc)
  Metabo.splsda <- splsda(t(current_e_plsda), Metabol_genderFc)
  # plotIndiv(Metabo.splsda)    
  p <- plotIndiv(Metabo.splsda, group = Metabol_genderFc, ind.names = FALSE,
                 ellipse = TRUE, legend = TRUE, legend.position = "right", title = paste("PLSDA of", current_gene))
  x11(width=10,height=8)
  print(p)
  savePlot(file=paste("Figure3d",paste(current_gene,"PLSDA"),".emf"))
  dev.off()
}  
#End of Figure 3d