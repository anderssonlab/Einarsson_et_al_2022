library(ggplot2)
library(GenomicRanges)
library(SummarizedExperiment)
library(GenomicFeatures)
library(BiocParallel)
library(InteractionSet)
library(CAGEfightR)



#Rename data packages
txdb <- loadDb("gencode.v29.annotation.sqlite")
genomeInfo <- keepStandardChromosomes(SeqinfoForUCSCGenome("hg38"),species="Homo_sapiens")
seqlevels(txdb) <- c(seqlevels(genomeInfo))  #Keep only standard chromosomes

## Download CAGEfightR extension
system("wget https://github.com/anderssonlab/CAGEfightR_extensions/archive/v0.1.1.tar.gz")
system("tar xvzf v0.1.1.tar.gz")
system("mv CAGEfightR_extensions-0.1.1 CAGEfightR_extensions")

source("CAGEfightR_extensions/decompose.R")
source("CAGEfightR_extensions/utils.R")



  ### CAGE clustering, quantification and annotation ###

#Format bw files
plus_files <- system("ls /binf-isilon/alab/projects/reg_variation/LCL_full/CAGE/MAPPING_bwa_hg38/WASP_biallelic/POOLED/bw_files/*plus.bw", intern=TRUE)
minus_files <- system("ls /binf-isilon/alab/projects/reg_variation/LCL_full/CAGE/MAPPING_bwa_hg38/WASP_biallelic/POOLED/bw_files/*minus.bw", intern=TRUE)
names(plus_files) <- names(minus_files) <- as.character(sapply(plus_files, function(n) gsub(".plus.bw","",gsub("hg38.","",gsub("_wasp","",gsub("/binf-isilon/alab/projects/reg_variation/LCL_full/CAGE/MAPPING_bwa_hg38/WASP_biallelic/POOLED/bw_files/","",n))))))

keep <- grep("_0",names(plus_files))
plus_files <- plus_files[keep]
minus_files <- minus_files[keep]

bw_plus <- BigWigFileList(plus_files)
bw_minus <- BigWigFileList(minus_files)


#Design matrix
design <- read.table("design.txt", header=T, sep="\t",stringsAsFactors = FALSE)


#Quantify CTSSs
register(MulticoreParam(workers=30))

CTSSs <- quantifyCTSSs(plusStrand = bw_plus,
                       minusStrand = bw_minus,
                       genome=genomeInfo,
                       design = design)

CTSSs <- CTSSs %>% calcPooled(inputAssay="counts")

supported.CTSSs <- subsetBySupport(CTSSs,minSamples=4)

supported.CTSS <- CTSSs %>% calcPooled(inputAssay="counts")

#First clustering step
TCs <- clusterUnidirectionally(supported.CTSSs,
                               pooledCutoff=0,
                               mergeDist=60)

TCs <- quantifyClusters(supported.CTSSs,TCs)


TCs <- assignTxType(TCs, txModels = txdb)

TCs<- subset(TCs, txType %in% c("promoter","proximal","fiveUTR"))



TCs <- subsetBySupport(TCs, inputAssay = "counts", unexpressed = 10, minSamples = 10) 

TCs <- assignGeneID(TCs,geneModels=txdb)

#Decomposing
decomposed.TCs <- decompose(rowRanges(TCs),supported.CTSSs,fn=local_maxima_decompose,smoothPad=1,fraction=0.1,mergeDist=1,maxGap=10) # lack of merging? Issue with ends trimming?! Very few not the same

#Save files
save(TCs, decomposed.TCs, supported.CTSSs,file="../data/Clusters.and.CTSSs.RData")





