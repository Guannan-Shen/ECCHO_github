# ECCHO_Guannan

ECCHO grant (chemicals and methylation project of the Healthy Start cohort).

## 9/27/2018
linear regression model for 9 outcomes on obesity and CpG  
$ y = CpG + race + CellTypes (Male)$  
Outcomes: "birth_weight", "ipv3_pp_fm_pct", "Chol_IPV3", "FFA_IPV3", "Gluc_IPV3", "HDL_IPV3", "Insu_IPV3", "Trig_IPV3", "Leptin_actual__ng_ml_".  
Although we discussed include the variable “infant sex” as a covariate, we actually did the chemicals-methylation analysis stratified by sex. So the DMRs were only detected among male offspring. Therefore we should probably restrict our analysis to individuals with infant_sex = 2 (males), at least to begin with.  
There will be missing data for some outcomes, I think it is OK if the sample size differs for each model.  

### Data set
1. HS_450K_CB_Mval_normbatch_StarlingSubset_10-01-18: a csv file which is the subset of the M-values of 600 participants from the full dataset based on the top 300 CpG names.  
2. healthy_start_cordblood_cellcounts: here is the cell type dataset from Weiming. Columns are cell type and rows are subjectID.  
3. pfas_methyl_di (PFAS-methylation outcomes list Sep 2018): limited dataset containing the outcomes and covariates for the analysis we discussed, and a data dictionary explaining the variable names.  

## 01/19/2019
Redo the DMP, DMR and Correlation analysis.  
### All analysis are done for male and female, separately.  
**DMP**: Differentially methylated Probes.  
$ Mvalue of individual CpG = ln(Chemical concentration) + race (Binary) + CellType1 +  CellType2 ... + CellType7 $  
**DMR**: Differentially methylated regions.  
DMR was done by DMRcate package in R, pcutoff=0.05 was set as the thresholding. Initial analysis used fdr = 0.20 to have more than 300 significant DMRs.  
**Obesity Parameters**  
$ log10(Obesity outcomes) = Mvalue of individual CpG + race (Binary) + CellType1 +  CellType2 ... + CellType7 $  

