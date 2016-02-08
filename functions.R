##########################################################

# DE analysis functions

##########################################################

# setwd("/Users/choonoo/expression_scripts")

# source("http://bioconductor.org/biocLite.R")

#biocLite('BiocInstaller')
require(gdata)
require(plyr)
require(oligo)
require(pd.mogene.2.1.st)
require(Heatplus)
require(reshape2)
require(ReportingTools)
require(hwriter)
require(ggplot2)
require(limma)
require(clusterProfiler)
require(XML)
require(mogene21sttranscriptcluster.db)
require(GOstats)
require(genefilter)
require(dplyr)

##########################################################

# function definitions

##########################################################


#from genefilter:::rowIQRs, just bringing it into this namespace
rowIQRs <- function (eSet) 
{
  numSamp <- ncol(eSet)
  lowQ <- rowQ(eSet, floor(0.25 * numSamp))
  upQ <- rowQ(eSet, ceiling(0.75 * numSamp))
  upQ - lowQ
}

#Computes differential expression between 2 specified categories
#Input:
#       use.exprs: ExpressionSet with a pData slot containing a data.frame with Mating, Timepoint, Inf_Status and optionally Tissue
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

de.analysis <- function(use.exprs, contrast.correction="separate", category="Category_v2")
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
  
  res.list <- list(infecteds=decideTests(fit.2, method=contrast.correction), coefs=fit.2$coefficients)
  
  res.list <- list(res.list)
  
  names(res.list) <- colnames(fit.2$coefficients)
  
  #null.list <- sapply(res.list, is.null)
  
  return(list(results=res.list, fit=fit.2, cont.counts=cont.counts))
}

#Summarizes the result from de.bycross.to.mock
#Input:
#       res.list: A list as provided by the 'results' element of de.bycross.to.mock
#       by: Whether or not to annotated the 'summary.dta' data.frame with gene symbols or probeset IDs.
#       ...: Any filters to apply in the form: filter1=filter.vec, filter2=filter2.vec where filter.vec and filter2.vec are logical vectors named by probeset IDs.
#           TRUE should indicate that the probeset should be removed and FALSE should indicate that the value should be kept.
#           If not all the filters are FALSE for a given probeset, it will be removed from the 'summary.dta' element of the result, however it will only be flagged in the 'de.res' element.
#Output:
#       A list containing the following elements:
#           -de.res: A list named containing a data.frame containg the probeset, gene symbol, coeffecients (log2 scale) in the form: Inf_Timepoint.Inf_name.Mock_Timpoint.Mock_name.logFC
#                    and whether or not the comparison was significant named as before but with .Signif instead of .logFC.
#           -summary.dta: A data.frame containing either symbols or probeset IDs followed by a column or 1s and 0s for each Mating, indicating whether at least one contrast was significant over the timecourse
#                           after application of any supplied filters.

# Describe this function more


contrast.summary.tc <- function(res.list, by=c("symbol", "probe"), ...)
{
  #all rownames should be the same...
  
  by <- match.arg(by)
  
  filt.list <- list(...)
  
  annot.temp <- select(mogene21sttranscriptcluster.db, keys=rownames(res.list[[1]][[1]]), columns=c("SYMBOL"), keytype="PROBEID")
  stopifnot(anyDuplicated(annot.temp$SYMBOL) == 0)
  rownames(annot.temp) <- annot.temp$PROBEID
  
  cont.res.list <- lapply(res.list, function(x)
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
  
  cols <- lapply(cont.res.list, colnames)
  col.lens <- sapply(cols, length)
  col.max <- max(col.lens)
  
  exp.cols <- unique(unlist(cols))
  
  if (any(col.lens < col.max))
  {
    max.cols <- cols[[which.max(col.lens)]]
    
    for(i in which(col.lens < col.max))
    {
      diff.cols <- setdiff(max.cols, cols[[i]])
      fill.mat <- matrix(NA, ncol=length(diff.cols), nrow=nrow(cont.res.list[[i]]), dimnames=list(NULL, diff.cols))
      
      cont.res.list[[i]] <- cbind(cont.res.list[[i]], fill.mat)[,max.cols]
    }
  }
  
  cont.res.list <- cont.res.list[sapply(cont.res.list, nrow) > 0]
  
  cont.summary.dta <- data.frame(do.call("rbind", lapply(1:length(cont.res.list), function(x) cbind(cont.res.list[[x]], Mating=names(cont.res.list)[x]))))
  
  if (length(filt.list) > 0)
  {
    filt.mat <- apply(do.call("cbind", filt.list), 1, function(x) all(x == F))
    cont.summary.dta <- cont.summary.dta[cont.summary.dta$ProbeId %in% names(filt.mat)[filt.mat],]
  }
  
  use.form <- formula(paste0(switch(by, probe="ProbeId", symbol="Symbol"), "~Mating"))
  
  cont.summary <- dcast(use.form, data=cont.summary.dta, fill=0, value.var="Mating", fun.aggregate=length)
  
  return(list(de.res=cont.res.list, summary.dta=cont.summary))
}


normalize_expression = function(raw.exprs=raw.exprs){
  # Compute RMA summarized transcript cluster data using annotated genes according to affy
  
  norm.exprs = rma(raw.exprs, normalize=TRUE, target="core")
  
  # save normalized expression to file
  norm.exprs_file = exprs(norm.exprs)
  
  # insert folder location to save normalized expression
  reports.dir = "./"
  
  # save file
  print("Saving normalized expression to file...")
  write.csv(file=paste(reports.dir,"normalized_expression.csv",sep="/"), x=norm.exprs_file, quote=F)
  
  # return normalized expression
  return(norm.exprs)
}

filter_features = function(norm.exprs=norm.exprs){
  # Map probe IDs to Entrez ID
  suppressWarnings(probe.ents <- select(get('mogene21sttranscriptcluster.db'), keys=featureNames(norm.exprs), columns=c("ENTREZID"), keytype="PROBEID"))
  
  # Filter NA
  use.probe.ents <- probe.ents[!is.na(probe.ents$ENTREZID),]
  
  # Remove duplicated mappings
  sub.eset <- norm.exprs[unique(use.probe.ents$PROBEID),]
  
  # Filter low varying expression values
  unique.probes <- findLargest(featureNames(sub.eset), rowIQRs(sub.eset), 'mogene21sttranscriptcluster.db')
  
  # Code chunk to output mapping of entrez to probe used in de analysis if necessary
  # use.probe.ents[as.vector(unlist(sapply(unique.probes, function(x)which(x==use.probe.ents[,1])))),] -> probe_mapping
  
  # Filter expression set to include unique probes
  norm.exprs.filter <- sub.eset[unique.probes,]
  
  # return filered data
  return(norm.exprs.filter)
}

sample_filter = function(norm.exprs.filter=norm.exprs.filter){
  # Annotate category to use for DE analysis
  pData(norm.exprs.filter)$Category = "n/a"
  
  pData(norm.exprs.filter)[which(pData(norm.exprs.filter)[,"D4_percent"] <.85),"Category"] <- "Sensitive"
  
  pData(norm.exprs.filter)[which(pData(norm.exprs.filter)[,"D4_percent"] >.98),"Category"] <- "Resistant"
  
  # Filter expression set to include annotated samples
  norm.exprs.filter.category <- norm.exprs.filter[, norm.exprs.filter$Category != "n/a"]
  
  # return filtered samples
  return(norm.exprs.filter.category)
}

de_analysis_table = function(norm.exprs.filter.category=norm.exprs.filter.category, category='Category'){
  all.res.list <- de.analysis(use.exprs=norm.exprs.filter.category, contrast.correction="separate", category=category)
  
  html.res <- contrast.summary.tc(all.res.list$results)
  
  # Save DE table of results
  
  html.res$de.res[[1]] -> de_matrix
  
  rows = sapply(row.names(de_matrix),function(x)which(x==row.names(all.res.list$fit$p.value)))
  
  pvalue = all.res.list$fit$p.value
  
  pvalue_v2 = as.matrix(pvalue[as.vector(unlist(rows)),])
  
  colnames(pvalue_v2) = "p.value"
  
  sum(row.names(de_matrix) == row.names(pvalue_v2))
  
  de_table = cbind(de_matrix, pvalue_v2)
  
  # Order by pvalue
  
  de_table[order(de_table[,"p.value"]),] -> de_table_v2
  
  # Annotate category in expression file
  
  sum(colnames(norm.exprs.filter.category) == row.names(pData(norm.exprs.filter.category))) == length(colnames(norm.exprs.filter.category))
  
  colnames(norm.exprs.filter.category) <- paste(colnames(norm.exprs.filter.category),pData(norm.exprs.filter.category)[,category],sep="_")
  
  # Add group means
  resistant_mean = rowMeans(norm.exprs_file[,grep(sapply(strsplit(names(de_matrix)[3],"\\."),"[",1),colnames(norm.exprs_file))])
  
  sensitive_mean = rowMeans(norm.exprs_file[,grep(sapply(strsplit(names(de_matrix)[3],"\\."),"[",2),colnames(norm.exprs_file))])
  
  norm.exprs_file_means = cbind(norm.exprs_file,resistant_mean, sensitive_mean)
  
  norm.exprs_file_means[,c("resistant_mean","sensitive_mean")] -> de_means
  
  rows2 = sapply(row.names(de_table_v2),function(x)which(x==row.names(de_means)))
  
  de_means[as.vector(unlist(rows2)),] -> de_means_v2
  
  sum(row.names(de_table_v2) == row.names(de_means_v2))
  
  de_table_v3 = cbind(de_table_v2, de_means_v2) # downregulated in resistant and upregulated in sensitive
  
  # Add adjusted pvalue
  de_table_v3$p.value.BY = p.adjust(de_table_v3[,"p.value"],method="BY")
  
  # Add entrez symbol
  names(probe_mapping)[1] <- "ProbeId"
  
  merge(de_table_v3, probe_mapping, by="ProbeId",all.x=T) -> de_table_v4
  
  # Write table to file
  print("Saving DE Table to file...")
  write.table(file=paste0(reports.dir,"/de_table.txt"), x=de_table_v4, sep="\t",quote=F)
  
  # return DE table
  return(de_table_v4)
}

# add parameter for top ranking DE genes (pvalue, fc > x, top 100 genes)


pathway_analysis = function(norm.exprs.filter.category=norm.exprs.filter.category, de_table=de_table, pvalue_cutoff=0.05){
  # set universe of genes
  univ <- select(mogene21sttranscriptcluster.db, keys=featureNames(norm.exprs.filter.category), columns=c("ENTREZID"), keytype="PROBEID")
  
  # Get list of DE genes, 
  selectedEntrezIds = de_table[which(de_table[,4] !=0),"ENTREZID"]
  
  # GO enrichment using hypergeometic distribution
  
  # Set significant pvalue < 0.05, Over-represented genes in BP = biological processes
  
  hgCutoff <- pvalue_cutoff
  params <- new("GOHyperGParams",
                geneIds=selectedEntrezIds,
                universeGeneIds=univ,
                annotation="mogene21sttranscriptcluster.db",
                ontology="BP",
                pvalueCutoff=hgCutoff,
                conditional=FALSE,
                testDirection="over")
  
  hpOver = hyperGTest(params)
  
  #save(hpOver, file=paste0(reports.dir,"/go_results.RData"))
  
  path_results = summary(hpOver)
  
  #htmlReport(hpOver, file=paste0(reports.dir,"/go_results.html"))
  
  # Annotate gene ratio
  generatio = path_results[,"Count"]/path_results[,"Size"]
  
  cbind(path_results, generatio) -> path_results_v2
  
  # add adjusted pvalue
  path_results_v2$Pvalue.adjusted = p.adjust(path_results_v2[,"Pvalue"],method="BY")
  
  head(path_results_v2)
  
  print("Saving Pathway Results...")
  # Write pathway results to file
  write.table(file=paste0(reports.dir,"/go_results.txt"), x=path_results_v2, sep="\t",quote=F,row.names=F)
  
  return(path_results_v2)
}

pathway_plots = function(path_results=path_results,parameter=c("OddsRatio","generatio"), pvalue=c('raw','adjusted')){
  
  parameter <- match.arg(parameter)
  
  pvalue <- match.arg(pvalue)
  
  # subset top 10
  res_table_v2 = path_results[1:10,]
  
  res_table_v2[order(res_table_v2[,parameter]),] -> res_table_v3
  
  if(pvalue == 'raw'){
     
    p = ggplot(res_table_v3, aes(res_table_v3[,parameter], reorder(res_table_v3$Term,res_table_v3[,parameter])))

    p + geom_point(aes(size = Count, colour = res_table_v3$Pvalue)) + scale_colour_gradient(low = "red", high="blue") + labs(y='',x=as.character(parameter))

    p + geom_point(aes(size = res_table_v3$Count, colour = res_table_v3$Pvalue)) + scale_colour_gradient(low = "red", high="blue") + labs(y='',x=as.character(parameter))
    
  }else{
    p = ggplot(res_table_v3, aes(res_table_v3[,parameter], reorder(res_table_v3$Term,res_table_v3[,parameter])))
    p + geom_point(aes(size = Count, colour = res_table_v3$Pvalue.adjusted)) + scale_colour_gradient(low = "red", high="blue") + labs(y='',x=as.character(parameter))
  }
  

}


