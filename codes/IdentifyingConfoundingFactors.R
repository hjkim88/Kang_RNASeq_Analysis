###
#   File name : IdentifyingConfoundingFactors.R
#   Author    : Hyunjin Kim
#   Date      : Oct 18, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Run PCA and to find possible confounding factors
#
#   Instruction
#               1. Source("IdentifyingConfoundingFactors.R")
#               2. Run the function "identifyCF" - specify an input file path (raw counts) and output directory
#               3. The PCA plots will be generated under the output directory
#
#   Example
#               > source("The_directory_of_IdentifyingConfoundingFactors.R/IdentifyingConfoundingFactors.R")
#               > identifyCF(rawCntPath="./data/raw_counts.rda",
#                            outputDir="./results/pca/")
###

identifyCF <- function(rawCntPath="./data/raw_counts.rda",
                       outputDir="./results/pca/") {
  
  ### load dataset
  load(rawCntPath)
  
  ### make sample info
  sampleInfo <- c(rep("ER_High", 3), rep("ER_Low", 3))
  
  ### A function to transform RNA-Seq data with VST in DESeq2 package
  normalizeRNASEQwithVST <- function(readCount, filtering=TRUE) {
    
    ### load library
    if(!require(DESeq2)) {
      source("https://bioconductor.org/biocLite.R")
      biocLite("DESeq2")
      library(DESeq2)
    }
    
    ### make a design matrix for DESeq2 data
    condition <- data.frame(factor(rep("OneClass", ncol(readCount))))
    
    ### Data preparation for DESeq2 format
    deSeqData <- DESeqDataSetFromMatrix(countData=readCount, colData=condition, design= ~0)
    
    if(filtering == TRUE) {
      ### Remove rubbish rows - this will decrease the number of rows
      deSeqData <- deSeqData[rowSums(counts(deSeqData))>1,]
    }
    
    ### VST
    vsd <- vst(deSeqData)
    transCnt <- data.frame(assay(vsd), check.names = FALSE)
    
    return (transCnt)
    
  }
  
  ### A function to perform 2D PCA and save a plot
  ### normalizedMat: rows are genes and columns are samples
  ### grp: group information of the samples
  ### num: the number of top genes to be used based on variance (-1 [default]: use all the genes)
  ### component: to draw a plot with PC1 & PC2 or PC2 & PC3
  ### title: title of the plot
  ### outDir: output directory for the plot
  pca_plot <- function(normalizedMat, grp,
                       num = -1, component=c("PC1&PC2", "PC2&PC3"),
                       title="PCA_Plot", outDir="./") {
    
    ### load library
    if(!require(ggfortify)) {
      install.packages("ggfortify")
      library(ggfortify)
    }
    
    ### select the top genes based on variance
    if(num >= 0 && num <= nrow(normalizedMat)) {
      v <- apply(normalizedMat, 1, var)
      v <- v[order(-v)]
      top_genes <- names(v)[1:num]
    } else {
      top_gnes <- rownames(normalizedMat)
    }
    
    ### PCA
    pca_result <- prcomp(t(normalizedMat[top_genes,]))
    pca_group <- data.frame(pca_result$x, group=grp)
    
    colors = topo.colors(length(unique(grp)))
    names(colors) = unique(grp)
    
    ### save as png
    if(component[1] == "PC1&PC2") {
      ggplot(pca_group,aes(x=PC1,y=PC2,col=group)) +
        labs(title=paste0(title, "_PC1-2")) +
        geom_text(aes(label=colnames(normalizedMat)),hjust="inward", vjust="inward") +
        scale_color_manual(values = colors) +
        theme_classic(base_size = 16)
      ggsave(filename = paste0(outDir, title, "_PC1-2_", ".png"), width = 10, height = 8)
    } else if(component[1] == "PC2&PC3") {
      ggplot(pca_group,aes(x=PC2,y=PC3,col=group)) +
        labs(title=paste0(title, "_PC2-3")) +
        geom_text(aes(label=colnames(normalizedMat)),hjust="inward", vjust="inward") +
        scale_color_manual(values = colors) +
        theme_classic(base_size = 16)
      ggsave(filename = paste0(outDir, title, "_PC2-3_", ".png"), width = 10, height = 8)
    } else {
      stop("\"component\" parameter should be \"PC1&PC2\" or \"PC2&PC3\"")
    }
    
  }
  
  ### make a pca plot
  pca_plot(normalizeRNASEQwithVST(rawCnt[,-c(1,2)]), grp = sampleInfo, num = 1000,
           title = "PCA_plot_ER", outDir = outputDir)
  
}
