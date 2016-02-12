##########################################################

# DE and Pathway Analysis Functions

##########################################################

# setwd("/Users/choonoo/DE_array_analysis")

# source("http://bioconductor.org/biocLite.R")

#biocLite('BiocInstaller')
require(gdata)
require(plyr)
require(oligo)
require(pd.mogene.2.1.st)
require(Heatplus)
require(reshape2)
require(hwriter)
require(ggplot2)
require(limma)
require(clusterProfiler)
require(XML)
require(mogene21sttranscriptcluster.db)
require(GOstats)
require(genefilter)
require(AnnotationDbi)

##########################################################

# Optional save norm expression
# Input: 
  # norm.exprs: normalized expression
  # save.dir: path file to save data
# Output: Saved normalized expression file
output_normalized_expression = function(norm.exprs=norm.exprs, save.dir='.'){
  
  # save normalized expression to file
  norm.exprs_file = exprs(norm.exprs)
  
  # insert folder location to save normalized expression
  #reports.dir = "./"
  
  # save file
  write.csv(file=file.path(save.dir,"normalized_expression.csv"), x=norm.exprs_file, quote=F)
  
  writeLines(paste0("Saving normalized expression to file: \n", file.path(save.dir,'normalized_expression.csv')))
}

# Filter expression IQRs, used within the filter_features function
#from genefilter:::rowIQRs
# Input: 
  # eSet: Expression set
# Output: Difference between 75th and 25th quantile for each probe
rowIQRs <- function (eSet) 
{
  numSamp <- ncol(eSet)
  lowQ <- rowQ(eSet, floor(0.25 * numSamp))
  upQ <- rowQ(eSet, ceiling(0.75 * numSamp))
  upQ - lowQ
}

# Filter probe features
# Input:
  # norm.exprs: normalized expression set
  # save.dir: file path directory where to save probe mapping file
# Output: filtered normalized expression set
filter_features = function(norm.exprs=norm.exprs, save.dir='.'){
  
  # Map probe IDs to Entrez ID
  suppressWarnings(probe.ents <- select(mogene21sttranscriptcluster.db, keys=featureNames(norm.exprs), columns=c("ENTREZID"), keytype="PROBEID"))
  
  # Filter NA
  use.probe.ents <- probe.ents[!is.na(probe.ents$ENTREZID),]
  
  # Remove duplicated mappings
  sub.eset <- norm.exprs[unique(use.probe.ents$PROBEID),]
  
  # Multiple probes per gene so this function finds all replicates and selects ones with largest value of the test statistic
  unique.probes <- findLargest(featureNames(sub.eset), rowIQRs(sub.eset), 'mogene21sttranscriptcluster.db')
  
  # Code chunk to output mapping of entrez to probe used in de analysis if necessary
  use.probe.ents[as.vector(unlist(sapply(unique.probes, function(x)which(x==use.probe.ents[,1])))),] -> probe_mapping
  
  
  write.table(file=file.path(save.dir, "probe_mapping.txt"), probe_mapping, row.names=F, quote=F, sep="\t")
  writeLines(paste0("Saving probe mapping to file: \n", file.path(save.dir,'probe_mapping.txt')))
  
  # Filter expression set to include unique probes
  norm.exprs.filter <- sub.eset[unique.probes,]
  
  # return filered data
  return(norm.exprs.filter)
}

#Computes differential expression between 2 specified categories
#Input:
#       use.exprs: normalized expression
#       contrast.correction: How to correct for multiple testing should be a valid value to provide to the limma::decideTests function
#       category: column of pData to use for de analysis between categories        
#Output:
#       A list containing the following elements:
#           -results: A list containing the following elements named by Mating:
#               1) coefs: The contrast coefficients
#               2) infecteds: The matrix of 0,1 and -1's formed from running limma::decideTests with method=contrast.correction
#           -fit: The MArrayLM object from limma.
#           -cont.counts: The number of samples for the two parts of each formed contrast or NA if the samples were not available

# add timepoint contrasts within this function

de.analysis <- function(use.exprs, contrast.correction="separate", category="Category")
{
  
  exprs.dta <- pData(use.exprs)
  
  category = pData(use.exprs)[,names(pData(use.exprs)) == category]
  
  contrasts <- with(exprs.dta, paste(category, sep="."))
  
  cont.counts = mapply(function(x,y) return(c(x,y)), as.numeric(table(contrasts)[names(summary(as.factor(contrasts)))[2]]), as.numeric(table(contrasts)[names(summary(as.factor(contrasts)))[1]]), SIMPLIFY=F)
  
  mouse.conts = paste(names(summary(as.factor(contrasts)))[1], names(summary(as.factor(contrasts)))[2],sep="-")
  #mouse.conts = c("Resistant-Sensitive")
  
  names(cont.counts) <- mouse.conts
  
  mod <- model.matrix(~0+contrasts)
  
  colnames(mod) <- sub("contrasts", "", colnames(mod)) 
  
  fit.1 <- lmFit(use.exprs, mod)
  
  conts <- makeContrasts(contrasts=mouse.conts, levels=mod)
  
  fit.2 <- contrasts.fit(fit.1, conts)
  
  fit.2 <- eBayes(fit.2)
  
  unique.matings <- unique(contrasts)
  
  de_results <- list(infecteds=decideTests(fit.2, method=contrast.correction), coefs=fit.2$coefficients)
  
  de_results <- list(de_results)
  
  names(de_results) <- colnames(fit.2$coefficients)
  
  #null.list <- sapply(de_results, is.null)
  
  return(list(results=de_results, fit=fit.2, cont.counts=cont.counts))
}

#Summarizes DE results
#Input:
#       de_results: Results from the de.analysis function
#       by: summarize with gene symbols or probeset IDs.
#Output:
#       A list containing the following elements:
#           -de.res: data.frame containing the probeset, gene symbol, coeffecients (log2 scale) and whether or not the comparison was significant named as before but with .Signif instead of .logFC.
#           -summary.dta: A data.frame containing either symbols or probeset IDs followed by a column or 1s and 0s for each Mating, indicating whether at least one contrast was significant over the timecourse after application of any supplied filters.
de.summary <- function(de_results, by=c("symbol", "probe"), ...)
{
  
  by <- match.arg(by)
  
  filt.list <- list(...)
  
  annot.temp <- select(mogene21sttranscriptcluster.db, keys=rownames(de_results[[1]][[1]]), columns=c("SYMBOL"), keytype="PROBEID")
  
  stopifnot(anyDuplicated(annot.temp$SYMBOL) == 0)
  rownames(annot.temp) <- annot.temp$PROBEID
  
  de_results_v2 <- lapply(de_results, function(x)
  {
    stopifnot(all(rownames(x$coefs) == rownames(x$infecteds)))
    
    coef.mat <- x$coefs
    sig.mat <- x$infecteds
    
    #any.sig <- apply(sig.mat, 1, function(x) any(x != 0))
    
    #any.sig <- apply(sig.mat, 1, function(x) any(x != 0))
    
    colnames(coef.mat) <- paste(gsub("M\\d+[xX]\\d+\\.", "", colnames(coef.mat)), "logFC", sep=" ")
    colnames(sig.mat) <- paste(gsub("M\\d+[xX]\\d+\\.", "", colnames(sig.mat)), "Signif", sep=" ")
    
    temp.mat <- data.frame(ProbeId=rownames(coef.mat), Symbol=annot.temp[rownames(coef.mat),"SYMBOL"], coef.mat, sig.mat, stringsAsFactors=F)
    
    for(i in names(filt.list))
    {
      temp.mat <- cbind(temp.mat, as.integer(filt.list[[i]][rownames(temp.mat)]))
      names(temp.mat)[ncol(temp.mat)] <- i
    }
    
    return(temp.mat)
  })
  
  #see if all columns are the same, if not add in the necessary columns with NAs
  
  cols <- lapply(de_results_v2, colnames)
  col.lens <- sapply(cols, length)
  col.max <- max(col.lens)
  
  exp.cols <- unique(unlist(cols))
  
  if (any(col.lens < col.max))
  {
    max.cols <- cols[[which.max(col.lens)]]
    
    for(i in which(col.lens < col.max))
    {
      diff.cols <- setdiff(max.cols, cols[[i]])
      fill.mat <- matrix(NA, ncol=length(diff.cols), nrow=nrow(de_results_v2[[i]]), dimnames=list(NULL, diff.cols))
      
      de_results_v2[[i]] <- cbind(de_results_v2[[i]], fill.mat)[,max.cols]
    }
  }
  
  de_results_v2 <- de_results_v2[sapply(de_results_v2, nrow) > 0]
  
  cont.summary.dta <- data.frame(do.call("rbind", lapply(1:length(de_results_v2), function(x) cbind(de_results_v2[[x]], Mating=names(de_results_v2)[x]))))
  
  if (length(filt.list) > 0)
  {
    filt.mat <- apply(do.call("cbind", filt.list), 1, function(x) all(x == F))
    cont.summary.dta <- cont.summary.dta[cont.summary.dta$ProbeId %in% names(filt.mat)[filt.mat],]
  }
  
  use.form <- formula(paste0(switch(by, probe="ProbeId", symbol="Symbol"), "~Mating"))
  
  cont.summary <- dcast(use.form, data=cont.summary.dta, fill=0, value.var="Mating", fun.aggregate=length)
  
  return(list(de.res=de_results_v2, summary.dta=cont.summary))
}

# Differental expression analysis
# Input:
  # norm.exprs: normalized expression set
  # category: Column of phenotype data of expression to that distinguishes the groups to compute differential expression
  # probe_mapping_file: probe mappings to gene IDs that was saved in the filter features function
#Output:
  # DE analysis table: ProbeId (Probe ID from expression set), Symbol (Gene Symbol), category.logFC (Fold change of mean expression between 2 categories), category.Signif (-1=downregulated,1=upregulated,0=not significant), p.value (raw pvalue), category1_mean (mean expression for category 1), category2_mean (mean expression for category 2), p.value.BY (adjusted pvalue), ENTREZID (Entrez ID)
de_analysis_table = function(norm.exprs, category, probe_mapping_file){
  
  # Compute DE analysis
  de_analysis <- de.analysis(use.exprs=norm.exprs, contrast.correction="separate", category=category)
  
  # Summarize DE results
  de_results <- de.summary(de_analysis$results)
  
  # Save DE table of results with pvalue, adjusted pvalue, and expression means for each group
  
  de_results$de.res[[1]] -> de_matrix
  
  # Add Pvalue
  rows = sapply(row.names(de_matrix),function(x)which(x==row.names(de_analysis$fit$p.value)))
  
  pvalue = de_analysis$fit$p.value
  
  pvalue_v2 = as.matrix(pvalue[as.vector(unlist(rows)),])
  
  colnames(pvalue_v2) = "p.value"
  
  sum(row.names(de_matrix) == row.names(pvalue_v2))
  
  de_table = cbind(de_matrix, pvalue_v2)
  
  # Order by Pvalue
  
  de_table[order(de_table[,"p.value"]),] -> de_table_v2
  
  # Annotate category
  
  sum(colnames(norm.exprs) == row.names(pData(norm.exprs))) == length(colnames(norm.exprs))
  
  colnames(norm.exprs) <- paste(colnames(norm.exprs),pData(norm.exprs)[,category],sep="_")
  
  norm_table = exprs(norm.exprs)
  
  # Add category means  
  group1 = sapply(strsplit(names(de_table_v2)[3],"\\."),"[",1)
  group1_mean = rowMeans(norm_table[,grep(group1, colnames(norm_table))])
  
  group2 = sapply(strsplit(names(de_table_v2)[3],"\\."),"[",2)
  group2_mean = rowMeans(norm_table[,grep(group2, colnames(norm_table))])
  
  de_means = data.frame(group1_mean, group2_mean)
  colnames(de_means) = c(paste0(group1, '_mean'), paste0(group2, '_mean'))
  
  de_means$ProbeId = row.names(de_means)
  
  merge(de_table_v2, de_means, by="ProbeId",all.x=T) -> de_table_v3
  
  # Add adjusted pvalue
  de_table_v3$p.value.BY = p.adjust(de_table_v3[,"p.value"],method="BY")
  
  # read in probe mapping
  read.table(file=probe_mapping_file, header=T, sep="\t") -> probe_mapping
  
  # Add entrez symbol
  names(probe_mapping)[1] <- "ProbeId"
  
  merge(de_table_v3, probe_mapping, by="ProbeId",all.x=T) -> de_table_v4
  
  # order by pvalue
  de_table_v4[order(de_table_v4[,"p.value.BY"]),] -> de_table_v5
  
  # Write table to file
  #print("Saving DE Table to file...")
  #write.table(file="./de_table.txt", x=de_table_v5, sep="\t",quote=F)
  
  # return DE table
  return(de_table_v5)
}

# GO enrichment
# Input: 
  #norm.exprs: normalized expression set
  # de_table: differentially expressed gene analysis table
  # path_pvalue: pvalue cutoff value for significantly enriched pathways
  # DE_p: consider adjusted or raw pvalue for DE genes to use in pathway analysis
  # DE_pvalue: DE pvalue cutoff value for significant DE genes
  # DE_fc: Fold change cutoff value for DE genes (subsets absolute FC changes greater than this value)
  # top_genes: The number of top genes to analyze path enrichment for (i.e. 10)
#Output:
  # Pathway analysis data frame: GOBPID (GO Biological Process ID), Pvalue (raw pvalue), OddsRatio (Odds that the pathway is enriched more than expected), ExpCount (Expected number of genes in path), Count (Number of DE genes in pathway), Size (Number of total genes in path), Term (GO BP full name), GeneRatio (Count/Size), Pvalue.adjusted (BY adjusted pvalue)
pathway_analysis = function(norm.exprs=norm.exprs.filter.category, de_table=de_table, path_pvalue=0.05, DE_p=c("adjusted","raw"), DE_pvalue=0.05, DE_fc=NULL, top_genes=NULL){
  
  # set universe of genes
  univ <- select(get('mogene21sttranscriptcluster.db'), keys=featureNames(norm.exprs), columns=c("ENTREZID"), keytype="PROBEID")
  
  # Match DE pvalue type (adjusted or raw)
  DE_p <- match.arg(DE_p)
  
  DE_p <- switch(DE_p,adjusted = 'p.value.BY', raw = "p.value")
  
  # subset DE pvalue
  de_table[which(de_table[,DE_p] < DE_pvalue),] -> de_table_v2
  
  # subset DE FC
  if(is.null(DE_fc)){
    de_table_v2 -> de_table_v3
  }else{
    de_table_v2[which(abs(de_table_v2[,3]) > DE_fc),] -> de_table_v3
  }
  
  # subset top number of genes if selected, uses all genes by all
  if(is.null(top_genes)){
    
    order(de_table_v3[,DE_p])[1:as.numeric(as.character(top_genes))] -> rows    
    
    de_table_v3[rows,] -> de_table_v4
    
  }else{
    de_table_v3 -> de_table_v4
  }
  
  # selected genes for pathway analysis
  selectedEntrezIds = as.character(de_table_v4[,"ENTREZID"])
  
  # GO enrichment using hypergeometic distribution, note: using over-represented genes in BP = biological processes
  params <- new("GOHyperGParams",
                geneIds=selectedEntrezIds,
                universeGeneIds=univ,
                annotation="mogene21sttranscriptcluster.db",
                ontology="BP",
                pvalueCutoff=path_pvalue,
                conditional=FALSE,
                testDirection="over")
  
  # summary of results
  hpOver = hyperGTest(params)
  
  path_results = summary(hpOver)
  
  # Annotate gene ratio
  GeneRatio = path_results[,"Count"]/path_results[,"Size"]
  
  cbind(path_results, GeneRatio) -> path_results_v2
  
  # add adjusted pvalue
  path_results_v2$Pvalue.adjusted = p.adjust(path_results_v2[,"Pvalue"],method="BY")
  
  head(path_results_v2)
  
  #print("Saving Pathway Results...")
  # Write pathway results to file
  #write.table(file="./go_results.txt", x=path_results_v2, sep="\t",quote=F,row.names=F)
  
  return(path_results_v2)
}

# Pathway Visualization
# Input:
  # path_results: Pathway analysis data frame
  # N: Number of pathways to display in plot
  # parameter: Display OddsRatio or GeneRatio
  # pvalue: Display adjusted or raw pvalue
# Output: Pathway Visualization
pathway_plots = function(path_results=path_results, N=10, parameter=c("OddsRatio","GeneRatio"), pvalue=c('adjusted','raw')){
  
  parameter <- match.arg(parameter, choices=c("OddsRatio","GeneRatio"))
  #print(parameter)
  
  pvalue <- match.arg(pvalue,choices=c('raw', 'adjusted'))
  
  # subset top 10
  res_table_v2 = path_results[1:N, ]
  
  res_table_v3 = res_table_v2[order(res_table_v2[, parameter]),]
  
  if (pvalue == 'raw') {
    p = ggplot(res_table_v3, aes_string(parameter, 'Term'))
    p + geom_point(aes(size = Count, colour = Pvalue)) + scale_colour_gradient(low = "red", high="blue") + labs(y='', x=as.character(parameter))   
  } else {
    p = ggplot(res_table_v3, aes_string(parameter, 'Term'))
    p + geom_point(aes(size = Count, colour = Pvalue.adjusted)) + scale_colour_gradient(low = "red", high="blue") + labs(y='', x=as.character(parameter))
  }
}