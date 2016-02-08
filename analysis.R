##########################################################

# DE analysis for Bat Virus Array

##########################################################

#setwd("/Users/choonoo/DE_array_analysis")

## Load libraries and functions for array processing and QA/QC plots
source('functions.R')

# save.image("~/de_test_v3.RData")

# load("./norm_exp.RData")

##########################################################

# Read in raw expression data

##########################################################

# save file where raw expression is saved

file = "./norm_exp.RData"

load(file)

##########################################################

# DE Analysis

##########################################################

# save normalized data
output_normalized_expression(norm.exprs=norm.exprs)

# feature filter
norm.exprs.filter = filter_features(norm.exprs=norm.exprs)

# sample filter
norm.exprs.filter.category = sample_filter(norm.exprs.filter=norm.exprs.filter)

# Compute DE analysis and save table of results
de_table = de_analysis_table(norm.exprs.filter.category=norm.exprs.filter.category, category='Category')

head(de_table)

##########################################################

# Pathway Analysis

##########################################################

# GO stats enrichment analysis of DE genes
path_results = pathway_analysis(norm.exprs.filter.category=norm.exprs.filter.category, de_table=de_table, pvalue_cutoff=0.05)

head(path_results)

# plots
parameter='OddsRatio'
pathway_plots(path_results=path_results,parameter=parameter,pvalue='raw')
# pathway_plots(path_results=path_results,parameter="OddsRatio",pvalue='adjusted')
# pathway_plots(path_results=path_results,parameter="generatio",pvalue='raw')
# pathway_plots(path_results=path_results,parameter="generatio",pvalue='adjusted')
