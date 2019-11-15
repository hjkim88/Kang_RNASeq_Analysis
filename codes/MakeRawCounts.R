###
#   File name : MakeRawCounts.R
#   Author    : Hyunjin Kim
#   Date      : Oct 17, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Do alignment, make Bam files, and make raw count matricies originally from fastq.gz files
#
#   * THIS CODE SHOULD BE RUN ON LINUX
#
#   * THIS IS MICE DATA, HENCE MM10 WILL BE USED AS REFERENCE
#
#   Instruction
#               1. Source("MakeRawCounts.R")
#               2. Run the function "makeRCnt" - specify the input directory (fastq.gz) and output directory
#               3. The Bam files and raw counts will be generated under the output directory
#
#   Example
#               > source("The_directory_of_MakeRawCounts.R/MakeRawCounts.R")
#               > makeRCnt(fastqgzPath="./data/lane_combined/",
#                          referencePath="E:/Reference/mm10.fa",
#                          referenceIdxPath="E:/Reference/mm10.index",
#                          readType=c("single-end", "paired-end"),
#                          outputDir="./data/")
###

makeRCnt <- function(fastqgzPath="/mnt/c/Research/CUMC/Kang_RNASeq_Analysis/data/lane_combined/",
                     referencePath="/mnt/e/Reference/mm10.fa",
                     referenceIdxPath="/mnt/e/Reference/mm10.index",
                     readType=c("single-end", "paired-end"),
                     outputDir="/mnt/c/Research/CUMC/Kang_RNASeq_Analysis/data/") {
  
  ### load library
  if(!require(Rsubread, quietly = TRUE)) {
    if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("Rsubread", version = "3.8")
    require(Rsubread, quietly = TRUE)
  }
  if(!require(org.Mm.eg.db, quietly = TRUE)) {
    if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("org.Mm.eg.db", version = "3.8")
    require(org.Mm.eg.db, quietly = TRUE)
  }
  
  ### build the index of the reference genome
  ### already built, therefore, no need to run this line
  # buildindex(basename=referenceIdxPath, reference=referencePath)
  
  ### get a list of fastq files
  fastqFiles <- list.files(fastqgzPath)
  fastqFiles <- fastqFiles[which(endsWith(fastqFiles, ".fastq.gz"))]
  
  ### create new directory for Bam files
  dir.create(paste0(fastqgzPath, "../bam_files/"), showWarnings = FALSE)
  
  ### iteratively perform alignment
  if(readType[1] == "paired-end") {
    sample_names <- unique(c(sapply(fastqFiles[grep("R1", fastqFiles)],
                                    function(x) {
                                      strsplit(x, "R1", TRUE)[[1]][1]                                      
                                    }),
                             sapply(fastqFiles[grep("R2", fastqFiles)],
                                    function(x) {
                                      strsplit(x, "R2", TRUE)[[1]][1]
                                    })))
    sample_names <- sapply(sample_names, function(x) substr(x, 1, nchar(x)-1), USE.NAMES = FALSE)
    
    for(samp in sample_names) {
      align(index=referenceIdxPath,
            readfile1=paste0(fastqgzPath, samp, "_R1.fastq.gz"),
            readfile2=paste0(fastqgzPath, samp, "_R2.fastq.gz"),
            output_file=paste0(fastqgzPath, "../bam_files/", samp, ".bam"),
            nthreads = 4)
    }
  } else if(readType[1] == "single-end") {
    for(i in 1:length(fastqFiles)) {
      align(index=referenceIdxPath,
            readfile1 = paste0(fastqgzPath, fastqFiles[i]),
            readfile2 = NULL,
            output_file = paste0(fastqgzPath, "../bam_files/", strsplit(fastqFiles[i], ".fastq", fixed = TRUE)[[1]][1], ".bam"),
            nthreads = 4)
    }
  } else {
    stop("The \"readType\" parameter should be \"paired-end\" or \"single-end\".")
  }
  
  ### get a list of created bam files
  bamFiles <- list.files(paste0(fastqgzPath, "../bam_files/"))
  bamFiles <- bamFiles[which(endsWith(bamFiles, ".bam"))]
  
  ### iteratively get counts from the bam files
  rawCnt <- NULL
  for(i in 1:length(bamFiles)) {
    counts <- featureCounts(files = paste0(fastqgzPath, "../bam_files/", bamFiles[i]),
                            annot.inbuilt = "mm10",
                            isPairedEnd = TRUE)
    
    if(is.null(rawCnt)) {
      rawCnt <- counts$counts
    } else {
      rawCnt <- cbind(rawCnt, counts$counts)
    }
  }
  
  ### numerize the raw counts (now they are characters)
  rawCnt <- as.data.frame(rawCnt)
  rawCnt[1:ncol(rawCnt)] <- lapply(rawCnt[1:ncol(rawCnt)], function(x) as.numeric(as.character(x)))
  
  ### set column names for rawCnt and annotate gene symbols
  colnames(rawCnt) <- substr(bamFiles, 1, nchar(bamFiles)-4)
  map_eg_symbol <- mappedkeys(org.Mm.egSYMBOL)
  list_eg2symbol <- as.list(org.Mm.egSYMBOL[map_eg_symbol])
  rawCnt <- cbind(Entrez_ID=rownames(rawCnt),
                  Gene_Symbol=as.character(list_eg2symbol[rownames(rawCnt)]),
                  rawCnt)
  
  ### write out the result
  write.table(rawCnt, file = paste0(outputDir, "raw_counts.txt"),
              sep = "\t", row.names = FALSE)
  save(list = c("rawCnt"), file = paste0(outputDir, "raw_counts.rda"))
  
}
