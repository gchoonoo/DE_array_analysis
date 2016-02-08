##########################################################

# DE analysis for Bat Virus Array

##########################################################

## Load libraries and functions for array processing and QA/QC plots
source('DE_pathway_functions.R')

# save.image("~/de_test_v3.RData")

# load("~/de_test_v3.RData")

##########################################################

# Read in raw expression data

##########################################################

# save file where raw expression is saved

file = "~/de_test_v2.RData"

load(file)

##########################################################

# DE Analysis

##########################################################

# normalize data and save
norm.exprs = normalize_expression(raw.exprs=raw.exprs)

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
pathway_plots(path_results=path_results,parameter="OddsRatio",pvalue='raw')
# pathway_plots(path_results=path_results,parameter="OddsRatio",pvalue='adjusted')
# pathway_plots(path_results=path_results,parameter="generatio",pvalue='raw')
# pathway_plots(path_results=path_results,parameter="generatio",pvalue='adjusted')

