
## Load libraries and functions for DE and pathway analysis
source('functions.R')

# save.image("~/DE_array_analysis/de_analysis_v2.RData")

# load("~/DE_array_analysis/de_analysis_v2.RData")

## Set data directory and load raw expression data
## You will have to change the file name and paths
data_dir = '.'
raw_exprs_file = 'de_test_v2.RData'
## This will load the raw.exprs object 
load(file.path(data_dir, raw_exprs_file))

## In this case we are selecting samples based on extreme values of the day 4 weight percentage
## We will create a 'Category' variable to indicate the two sample groups to compare
pData(raw.exprs)$Category = NA
pData(raw.exprs)[which(pData(raw.exprs)[,'D4_percent'] < 0.85),'Category'] <- 'Sensitive'
pData(raw.exprs)[which(pData(raw.exprs)[,'D4_percent'] > 0.98),'Category'] <- 'Resistant'
## Subset the ExpressionSet
raw.exprs.filter <- raw.exprs[, !is.na(raw.exprs$Category)]

dim(raw.exprs)

dim(raw.exprs.filter)

## Normalize the data
norm.exprs.filter = rma(raw.exprs.filter, normalize=TRUE, target="core")

dim(norm.exprs.filter)

# Save normalized expression as .csv file (optional)
output_normalized_expression(norm.exprs=norm.exprs.filter, save.dir=data_dir)

## This function removes transcript cluster IDs that do not map to an Entrez gene
## Also, if multiple transcript cluster IDs map to a gene, the one with the highest IQR is kept
norm.exprs.filter = filter_features(norm.exprs=norm.exprs.filter, save.dir=getwd())

dim(norm.exprs.filter)

# Compute DE analysis and save table of results
de_table = de_analysis_table(norm.exprs=norm.exprs.filter, category='Category', probe_mapping_file='probe_mapping.txt')

#write.table(file="./de_table.txt", x=de_table, sep="\t",quote=F)

# GO stats enrichment analysis of DE genes
path_results = pathway_analysis(norm.exprs=norm.exprs.filter, de_table=de_table, path_pvalue=0.05, DE_p="adjusted", DE_pvalue=0.05, DE_fc=NULL, top_genes=10)

head(path_results)

# Set parameter equal to either Odds Ratio (OddsRatio) or Gene Ratio (generatio)
pathway_plots(path_results=path_results, N=10, parameter='OddsRatio', pvalue='adjusted')
