##Figure 1 pie chart Metabolomics_sex difference  adjusted by BW

setwd("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
met_pie<-read.csv("Supplementary Data 3.csv")
nrow(met_pie)
names(met_pie)
#Data for Figure 1a
met_pie$Classify <- ifelse(((met_pie$Coefficient < 0 & met_pie$fold_change < 1)|(met_pie$Coefficient > 0 & met_pie$fold_change>1)), "CORRECT","Cannot classify" )
m_higher <- sum(met_pie$Classify == "CORRECT" & met_pie$p_value < 0.05 & met_pie$fold_change > 1)
f_higher <- sum(met_pie$Classify == "CORRECT" & met_pie$p_value < 0.05 & met_pie$fold_change < 1)
cannot_classify<-sum(met_pie$Classify == "Cannot classify" & met_pie$p_value < 0.05)
non_sig <- nrow(met_pie) - m_higher - f_higher - cannot_classify
piedta <- matrix(c(m_higher, f_higher,non_sig) , ncol = 1, byrow = TRUE, dimnames = list(c("m_higher", "f_higher", "non_sig"), c("NumberOf Metabolites")))
met_pie$Metbaolites[met_pie$Classify == "Cannot classify" & met_pie$p_value < 0.05]
met_pie$label[met_pie$Classify == "Cannot classify" & met_pie$p_value < 0.05]
# ----------------------------------------------------------------------------------------------------------------
##Data for Figure 1b Platform
###criteria p-value 0.05 threshold, using female as reference:
setwd("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
met_pie <- read.csv("Supplementary Data 3.csv")
nrow(met_pie)
GC_TOF <- sum(met_pie$p_value < 0.05 & met_pie$Platform %in% "gctof")
Lipidomics_Pos <- sum(met_pie$p_value < 0.05 & met_pie$Platform %in% "CSHPOS")
Lipidomics_Neg <- sum(met_pie$p_value < 0.05 & met_pie$Platform %in% "CSHNEG")
HILIC_Neg <- sum(met_pie$p_value < 0.05 & met_pie$Platform %in% "HILICPOS")
HILIC_Pos <- sum(met_pie$p_value < 0.05 & met_pie$Platform %in% "HILICNEG")
Oxylipins <- sum(met_pie$p_value < 0.05 & met_pie$Platform %in% "Oxylipins")
BA_Steroids <- sum(met_pie$p_value < 0.05 & met_pie$Platform %in% "BileAcid_Steroids")
Fig1b <- matrix(c(GC_TOF, Lipidomics_Pos, Lipidomics_Neg, HILIC_Neg, HILIC_Pos, Oxylipins, BA_Steroids), ncol = 1, byrow = TRUE, dimnames = list(c("GC_TOF","Lipidomics_Pos", "Lipidomics_Neg", "HILIC_Neg", "HILIC_Pos", "Oxylipins", "BA_Steroids"), c("NumberOf Metabolites")))

# ----------------------------------------------------------------------------------------------------------------
##Data for Figure 1c
###criteria p-value 0.05 for threshold, using female as reference:
met_pie <- met_pie[!is.na(met_pie$p_value),]
GC_TOF_sum <- sum(met_pie$Platform %in% "gctof")
GC_TOF_M <- sum(met_pie$p_value < 0.05 & met_pie$fold_change > 1 & met_pie$Platform %in% "gctof")
GC_TOF_F <- sum(met_pie$p_value < 0.05 & met_pie$fold_change < 1 & met_pie$Platform %in% "gctof")
Lipidomics_Pos_sum <- sum(met_pie$Platform %in% "CSHPOS")
Lipidomics_Pos_M <- sum(met_pie$p_value < 0.05 & met_pie$fold_change > 1 & met_pie$Platform %in% "CSHPOS")
Lipidomics_Pos_F <- sum(met_pie$p_value < 0.05 & met_pie$fold_change < 1 & met_pie$Platform %in% "CSHPOS")
Lipidomics_Neg_sum <- sum(met_pie$Platform %in% "CSHNEG")
Lipidomics_Neg_M <- sum(met_pie$p_value < 0.05 & met_pie$fold_change > 1 & met_pie$Platform %in% "CSHNEG")
Lipidomics_Neg_F <- sum(met_pie$p_value < 0.05 & met_pie$fold_change < 1 & met_pie$Platform %in% "CSHNEG")
HILIC_Neg_sum <- sum(met_pie$Platform %in% "HILICNEG")
HILIC_Neg_M <- sum(met_pie$p_value < 0.05 & met_pie$fold_change > 1 & met_pie$Platform %in% "HILICNEG")
HILIC_Neg_F <- sum(met_pie$p_value < 0.05 & met_pie$fold_change < 1 & met_pie$Platform %in% "HILICNEG")
HILIC_Pos_sum <- sum(met_pie$Platform %in% "HILICPOS")
HILIC_Pos_M <- sum(met_pie$p_value < 0.05 & met_pie$fold_change > 1 & met_pie$Platform %in% "HILICPOS")
HILIC_Pos_F <- sum(met_pie$p_value < 0.05 & met_pie$fold_change < 1 & met_pie$Platform%in% "HILICPOS")
Oxylipins_sum <- sum(met_pie$Platform %in% "Oxylipins")
Oxylipins_M <- sum(met_pie$p_value < 0.05 & met_pie$fold_change > 1 & met_pie$Platform %in% "Oxylipins")
Oxylipins_F <- sum(met_pie$p_value < 0.05 & met_pie$fold_change < 1 & met_pie$Platform %in% "Oxylipins")
BA_Steroids_sum <- sum(met_pie$Platform %in% "BileAcid_Steroids")
BA_Steroids_M <- sum(met_pie$p_value < 0.05 & met_pie$fold_change > 1 & met_pie$Platform %in% "BileAcid_Steroids")
BA_Steroids_F <- sum(met_pie$p_value < 0.05 & met_pie$fold_change < 1 & met_pie$Platform %in% "BileAcid_Steroids")
Fig1b <- matrix(c(GC_TOF_sum, GC_TOF_M, GC_TOF_F, Lipidomics_Pos_sum, Lipidomics_Pos_M, Lipidomics_Pos_F, Lipidomics_Neg_sum, Lipidomics_Neg_M, Lipidomics_Neg_F, HILIC_Neg_sum, HILIC_Neg_M, HILIC_Neg_F, HILIC_Pos_sum, HILIC_Pos_M, HILIC_Pos_F, Oxylipins_sum, Oxylipins_M, Oxylipins_F, BA_Steroids_sum, BA_Steroids_M, BA_Steroids_F), ncol = 1, byrow = TRUE, dimnames = list(c("GC_TOF_sum", "GC_TOF_M", "GC_TOF_F", "Lipidomics_Pos_sum", "Lipidomics_Pos_M", "Lipidomics_Pos_F", "Lipidomics_Neg_sum", "Lipidomics_Neg_M", "Lipidomics_Neg_F", "HILIC_Neg_sum", "HILIC_Neg_M", "HILIC_Neg_F", "HILIC_Pos_sum", "HILIC_Pos_M", "HILIC_Pos_F", "Oxylipins_sum", "Oxylipins_M", "Oxylipins_F", "BA_Steroids_sum", "BA_Steroids_M", "BA_Steroids_F")))
# =============================================================================================================================================
library(ggplot2)  # Data visualization
library(dplyr)    # Data manipulation
Mydata<-read.csv("Figure 1a.csv")
df <- data.frame(Mydata)
df
library(colorspace)
hcl_palettes(plot = T)
col1 <- diverge_hcl(3, "Cyan-Mage")
col1
fh <- col1[2]
ns <- col1[3]
mh <- col1[1]
col1 <- c(ns, mh, fh)
# set figure theme
blank_theme <- theme_minimal() + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size = 14, face = "bold")
  )

pie <- ggplot(df, aes(x = "", y = Percent005, fill = Group)) + blank_theme +
  theme(axis.text.x=element_blank())+
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = pi/2)+
  geom_text(aes(label = paste0(df$Percent005, "%")), position = position_stack(vjust = 0.5), size=5) +
  scale_fill_manual(values = col1)
pie
setwd("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
x11(w=8, h=8)
pie
savePlot("Figure 1a.emf", type="emf")
dev.off()

# Figure 1a
------------------------------------------------------------------------------------------------------------------------------
####Figure 1b pie chart_Platforms
library(colorspace)
library(RColorBrewer)
plot(1 : 10, rep(1, 10), col = rainbow(10, alpha = 0.5), cex = 5, pch = 15)
display.brewer.all()
col2 <- brewer.pal(8,"Set1")
plot(1 : 8, rep(1, 8), col = col2, pch = 16, cex = 2)
Mydata <- read.csv("Figure 1b_Platforms.csv")
df <- data.frame(Mydata)
df
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size = 14, face = "bold")
  )
pie <- ggplot(df, aes(x = " ", y = Percent005, fill = Platform)) + blank_theme +
  theme(axis.text.x=element_blank())+
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = pi/7.5, direction = 1)+
  geom_text(aes(label = paste0(Percent005, "%")), position = position_stack(vjust = 0.5), size=5) +
  scale_fill_manual(values = col2)
pie
x11(w = 8, h = 8)
pie
savePlot("Figure1b_platform.emf", type="emf")
dev.off()
####Figure 1b pie chart_Platforms
# ------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------
####Figure 1c
library(dplyr)   # For data manipulation
library(ggplot2) # For data visualization
Mydata<-read.csv("Figure 1c.csv",stringsAsFactors = F)
colnames(Mydata)
df<-as.data.frame(Mydata)
library(RColorBrewer)
display.brewer.all()
# Make a stacked barplot--> it will be in %!
library(colorspace)
hcl_palettes(plot=T)
col1<-diverge_hcl(3,'Cyan-Mage')
barplot(1:3,col=col1)
col1
fh<-col1[3]
ns<-col1[2]
mh<-col1[1]
col1<-c(ns,mh,fh)
bar <- ggplot(df, aes(y=Percent005, x=Platform, fill=Group)) +  
  geom_bar(stat = "identity", position = "stack", width = 0.9, colour = "black") +   
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        plot.margin = unit(c(1,1,1,0.3), "cm"))+
  scale_fill_manual("Group", values = c("Higher in female" = fh, "Higher in male" = mh, "Not significant" = ns))
x11(w=9, h=6)
bar
savePlot("Figure 1c_WT stack bar.emf",type="emf")
dev.off()
#### End of Figure 1c
# ------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------