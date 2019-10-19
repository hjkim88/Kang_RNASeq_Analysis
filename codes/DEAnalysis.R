###
#   File name : DEAnalysis.R
#   Author    : Hyunjin Kim
#   Date      : Oct 19, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Perform DE analysis on the raw counts
#
#   Instruction
#               1. Source("DEAnalysis.R")
#               2. Run the function "dea_kang" - specify raw count path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_DEAnalysis.R/DEAnalysis.R")
#               > dea_kang(rCntPath="./data/raw_counts.rda",
#                          outputDir="./results/differential_expression/")
###

dea_kang <- function(rCntPath="./data/raw_counts.rda",
                     outputDir="./results/differential_expression/") {
  
  ### load library
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  
  ### load raw count data
  load(rCntPath)
  
  ### make sample info
  sampleInfo <- c(rep("ER_High", 3), rep("ER_Low", 3))
  
  #####################################################
  ### A function to perform DE analysis with DESeq2 ###
  #####################################################
  #' @title deseqWithComparisons
  #' @param rCnt raw count matrix
  #' @param grp a character vector of class info of the samples
  #' @param exp_class a string of the experiment group's name
  #' @param ctrl_class a string of the control group's name
  #' @param bat_eff a character vector of batch effect info of the samples
  #' @param thresh numeric. Filters out from the results genes with adjusted
  #' 		p-value larger than this value
  #' @return data.frame
  #' @export
  #' @author Hyunjin Kim
  ####################################################
  deseqWithComparisons <- function(rCnt, grp, exp_class, ctrl_class, bat_eff=NULL, thresh = 1) {
    
    ### load library
    if(!require(DESeq2, quietly = TRUE)) {
      source("https://bioconductor.org/biocLite.R")
      biocLite("DESeq2")
      require(DESeq2, quietly = TRUE)
    }
    
    ### make a design matrix for DE analysis
    sampleType <- as.character(grp)
    
    if(is.null(bat_eff)) {
      Coldata <- data.frame(sampleType)
    } else {
      batch_eff <- as.character(bat_eff)
      Coldata <- data.frame(sampleType, batch_eff)
    }
    
    rownames(Coldata) <- colnames(rCnt)
    Coldata$sampleType <- relevel(Coldata$sampleType, ref = ctrl_class)
    
    ### data preparation for DE analysis
    if(is.null(bat_eff)) {
      deSeqData <- DESeqDataSetFromMatrix(countData=rCnt, colData=Coldata, design= ~sampleType)
    } else {
      deSeqData <- DESeqDataSetFromMatrix(countData=rCnt, colData=Coldata, design= ~sampleType+batch_eff)
    }
    
    deSeqData <- deSeqData[rowSums(counts(deSeqData))>1,]
    
    ### run DE analysis
    dea <- DESeq(deSeqData)
    deresults <- results(dea, contrast = c("sampleType", exp_class, ctrl_class))
    if(length(which(is.na(deresults$pvalue))) > 0) {
      deresults$pvalue[which(is.na(deresults$pvalue))] <- 1
    }
    if(length(which(is.na(deresults$padj))) > 0) {
      deresults$padj[which(is.na(deresults$padj))] <- 1
    }
    deresults <- deresults[order(deresults$padj, na.last = TRUE), ,drop = FALSE]
    deresults <- deresults[deresults$padj <= thresh, ,drop = FALSE]
    
    ### add baseMean for each group
    nCnt <- counts(dea, normalized=TRUE, replaced=TRUE)
    exp_rowMeans <- apply(nCnt[,which(Coldata$sampleType == exp_class), drop=FALSE], 1, mean)
    ctrl_rowMeans <- apply(nCnt[,which(Coldata$sampleType == ctrl_class), drop=FALSE], 1, mean)
    deresults <- data.frame(baseMean=deresults[,1],
                            V1=exp_rowMeans[rownames(deresults)],
                            V2=ctrl_rowMeans[rownames(deresults)],
                            deresults[,2:6],
                            stringsAsFactors = FALSE, check.names = FALSE)
    colnames(deresults)[2:3] <- c(paste0("baseMean_", exp_class), paste0("baseMean_", ctrl_class))
    
    return(deresults)
  }
  
  ### A function to print volcano plot of DE analysis with DESeq2 result
  volPlotWithDeseq <- function(deresult, outputFilePath, fdr=0.05, lfc=0.6) {
    
    ### load library
    if(!require(ggplot2)) {
      install.packages("ggplot2")
      library(ggplot2)
    }
    
    deresult$padj[which(is.na(deresult$padj))] <- 1
    volcanoData <- as.data.frame(cbind(deresult$log2FoldChange, -log10(deresult$padj), as.character(deresult$padj < fdr & abs(deresult$log2FoldChange) > lfc)))
    colnames(volcanoData) <- c("logFC", "logFDR", "Significance")
    volcanoData$logFC <- as.numeric(as.character(volcanoData$logFC))
    volcanoData$logFDR <- as.numeric(as.character(volcanoData$logFDR))
    
    s <- strsplit(outputFilePath, "/")
    f1 <- s[[1]][length(s[[1]])]
    
    ggplot(data=volcanoData, aes(x=logFC, y=logFDR, colour=Significance)) + xlab("log2_Fold_Change") + ylab("-log10(FDR)") + ggtitle(paste("Significant (adj.p < ", fdr, "& |logFC| > ", lfc, ") DE genes -", sum(volcanoData$Significance == "TRUE"))) + theme_classic(base_size = 16) + geom_point(alpha=0.4)
    ggsave(filename = outputFilePath, width = 12, height = 10)
  }
  
  ### run DE analysis
  deresult <- deseqWithComparisons(rCnt = rawCnt[,-c(1,2)], grp = sampleInfo,
                                   exp_class = "ER_High", ctrl_class = "ER_Low",
                                   bat_eff = NULL, thresh = 1)
  
  ### write out the DE result table
  write.xlsx2(data.frame(Entrez_ID=rownames(deresult), Gene_Symbol=rawCnt[rownames(deresult),"Gene_Symbol"], deresult,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "DE_Result.xlsx"),
              sheetName = "DE_Result", row.names = FALSE)
  
  ### draw a volcano plot
  volPlotWithDeseq(deresult, outputFilePath = paste0(outputDir, "Volcano_Plot.png"), fdr = 0.01, lfc = 1)
  
}
