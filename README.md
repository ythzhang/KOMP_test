# KOMP_test


### The content of this folder

The materials listed in this folder are provided to give statistical analysis for manuscript "Strong sexual dimorphism of plasma metabolites in 30 KO strains by multi-platform mass spectrometry analysis"

Note: 
1. Folders are neither linear nor independent.  Some figures/tables require multiple preceding steps which due to their size and their goal are housed in separate folders.

2. Within each folder is a README file which will tell you the inputs needed to achieve the goal of that folder and point you to other folders if necessary. 



```
Summary output files
1.	Supplementary Data 1.xlsx is the raw metabolomics data file for 220 mice (wildtype & KOs)
2.  Supplementary Data 2.xlsx is the raw phenotype data file extracted from IMPC database for 220 mice (wildtype & KOs)

This folder has the summary files obtained from the analysis collected together in one place.  Below is a quick summary of how to interpret the data and information on how the file was generated is provided.

```

Focusing on the role of sex in wildtype data-Mtabolomics
3.	Supplementary Data 3.xlsx
3a.	Spreadsheet 1:

•	This file is generated from the R.code: (1_Wildtype)
•	Column “p_value” gives the p value for the role of sex without adjustment
•	Column “Coefficient” or "fold_change" gives the comparison with male vs. female mice	
•	Column “adjusted_p_value” gives the p value for the role after adjusting using Benjamini & Hochberg method. For a 5% FDR look at rows with value <0.05. 

3b.	Spreadsheet 2:
•	This file is generated from the R.code: (1_Wildtype), and data was adjusted by bodyweight
Similar to 3a,
•	Column "p_value" gives the p value for the role of sex without adjustment
•	Column "Coefficient" or "fold_change" gives the comparison with male vs. female mice	
•	Column "adjusted_p_value" gives the p value for the role after adjusting using Benjamini & Hochberg method. For a 5% FDR look at rows with value <0.05. 


Focusing on the role of sex in wildtype data - Phenotype Data
4.	Supplementary Data 4.xlsx
4a.	Spreadsheet 1:
•	This file is generated from the R.code: (1_Wildtype)
•	Column "p_value" gives the p value for the role of sex without adjustment
•	Column "Coefficient" or "fold_change" gives the comparison with male vs. female mice for phenotype
•	Column "adjusted_p_value" gives the p value for the role after adjusting using Benjamini & Hochberg method. For a 5% FDR look at rows with value < 0.05. 

4b.	Spreadsheet 2:
•	This file is generated from the R.code: (1_Wildtype), and data was adjusted by bodyweight
Similar to 4a,
•	Column "p_value" gives the p value for the role of sex without adjustment
•	Column "Coefficient" or "fold_change" gives the comparison with male vs. female mice for phenotype	
•	Column "adjusted_p_value" gives the p value for the role after adjusting using Benjamini & Hochberg method. For a 5% FDR look at rows with value < 0.05. 


Focusing on the role of sex in the genotype effect - Metabolomics Data
5.	Supplementary Data 5.xlsx
5a.	Spreadsheet 1:
•	This file is generated from the the R.code: (5_Metabol_genoSD) and is for classifying the genotype-sex interaction effect with pvalue_Stage1 < 0.05
•	 The data is used to create bar plot in figure 3c
•	Column "pvalue_Stage1" gives the p value for the role of genotype without adjustment
•	Column "pvalue_Stage2" gives the  p value for the role of genotype-sex interaction without adjustment
•	Column "sep_FvKO_pval" and "sep_MvKO_pval" give the individual p values for comparing KO female vs. WT female or comparing KO male vs. WT male
•	Column "ctrl_sex_pval" and "KO_sex_pval" give the individual p values for comparing WT male vs. WT female comparing KO male vs. KO female, which is used for plotting boxplot in Figure 4, 5, and 6.

5b.	Spreadsheet 2:
•	This Spreadsheet is generated from the the R.code: (5_Metabol_genoSD) and is a summary of 30 KO starins to generate pie chart in Figure 3a
•	Column "Metabolites count" gives the total number of metabolites tested for each KO strain 
•	Column "Significant" gives the sum of all metabolites showed to be affected by genotype 
•	Other columns give detailed summary of metabolites under each classification


Focusing on the role of sex in the genotype effect - Phenotype Data
6.	Supplementary Data 6.xlsx
6a.	Spreadsheet 1:
•	This spreadsheet is generated from the the R.code: (5_PHENO_genoSD) and is for classifying the genotype-sex interaction effect with pvalue_Stage1 < 0.05
•	Column "pvalue_Stage1" gives the p value for the role of genotype without adjustment
•	Column "pvalue_Stage2" gives the  p value for the role of genotype-sex interaction without adjustment
•	Column "sep_FvKO_pval" and "sep_MvKO_pval" give the individual p values for comparing KO female vs. WT female or comparing KO male vs. WT male
•	Column "ctrl_sex_pval" and "KO_sex_pval" give the individual p values for comparing WT male vs. WT female comparing KO male vs. KO female


6b.	Spreadsheet 2:
•	This Spreadsheet is generated from the the R.code: (5_PHENO_genoSD) and is a summary of 30 KO starins to generate pie chart in Figure 3b
•	Column "Phenotype count" gives the total number of phenotypes tested for each KO strain 
•	Column "Significant" gives the sum of all phenotypes showed to be affected by genotype 
•	Other columns give detailed summary of phenotypes under each classification


