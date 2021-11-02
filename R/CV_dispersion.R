require(CAGEfightR)
require(recount)
require(ggplot2)
require(Seurat)
require(SingleCellExperiment)

## Download CAGEfightR extension
system("wget https://github.com/anderssonlab/CAGEfightR_extensions/archive/v0.1.1.tar.gz")
system("tar xvzf v0.1.1.tar.gz")
system("mv CAGEfightR_extensions-0.1.1 CAGEfightR_extensions")

source("CAGEfightR_extensions/dispersion.R")
source("CAGEfightR_extensions/normalize.R")



### CAGE DATA ###


## Load data
load("../data/Clusters.and.CTSSs.RData")


## Quantify expression
supported.CTSSs <- sort(sortSeqlevels(supported.CTSSs))

TC.expr <- quantifyClusters(supported.CTSSs, rowRanges(TCs))
TC.expr$totalTags <- supported.CTSSs$totalTags
TC.expr <- calcTPM(TC.expr,totalTags="totalTags")

decomposed.TC.expr <- quantifyClusters(supported.CTSSs, rowRanges(decomposed.TCs))
decomposed.TC.expr$totalTags <- supported.CTSSs$totalTags
decomposed.TC.expr <- calcTPM(decomposed.TC.expr,totalTags="totalTags")

supported.CTSSs <- calcTPM(supported.CTSSs,totalTags="totalTags")
supported.CTSSs <- sort(sortSeqlevels(supported.CTSSs))


## Calculate TC-level dispersion
TC.expr <- calcMedian(TC.expr,inputAssay="TPM",outputColumn="median")
TC.expr <- DMadjustedCV(TC.expr)


## Calculate decomposed TC-level dispersion
decomposed.TC.expr <- calcMedian(decomposed.TC.expr,inputAssay="TPM",outputColumn="median")
decomposed.TC.expr <- DMadjustedCV(decomposed.TC.expr)


# Save data
save(TC.expr,decomposed.TC.expr,file="../data/promoter_CTSS_adjusted_CV_data.RData")



### GEUVADIS RNA-seq data ###


if( ! file.exists("ERP001942/rse_gene.Rdata") ){
    download_study('ERP001942')
}

load(file.path('ERP001942', 'rse_gene.Rdata'))
geuvadis <- scale_counts(rse_gene)


info = read.table("https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt", header=TRUE, sep="\t")
colnames(info) = gsub("Comment\\.", "", colnames(info))
colnames(info) = gsub("Characteristics\\.", "", colnames(info))
colnames(info) = gsub("\\.$", "", colnames(info))

info = unique(info[,!colnames(info) %in% c('Scan.Name', 'SUBMITTED_FILE_NAME', 'FASTQ_URI')])
info = subset(info, Factor.Value.ancestry.category=="Yoruba")

metadata = merge(colData(geuvadis), info, by.x="run", by.y='ENA_RUN')
metadata = data.frame(metadata)
rownames(metadata) = metadata$run

geuvadis = geuvadis[,colnames(geuvadis) %in% metadata$run]
geuvadis = geuvadis[,match(metadata$run, colnames(geuvadis))]

## Filter expressed (>1 TPM) in at least 10% of samples
assay(geuvadis,"TPM") = getTPM(geuvadis)
isexpr = rowSums(assay(geuvadis,"TPM")>1) >= 0.1*ncol(assay(geuvadis))
geuvadis = geuvadis[isexpr,]

## Calculate dispersion
geuvadis <- calcMedian(geuvadis,inputAssay="TPM",outputColumn="median")
geuvadis <- DMadjustedCV(geuvadis)

## Save data
save(geuvadis,file="../data/GEUVADIS_Yoruba_adjusted_CV_data.RData")



### scRNA-seq data ###


## Acquire data
GM12878 <- Read10X(data.dir = "GSM3596321/GM12878")
GM12878 <- SummarizedExperiment(list(counts=GM12878))
GM12878 <- as(GM12878, "SingleCellExperiment")

## Filter cells
meta <- read.table("GSM3596321/GSM3596321_GM12878_cellQC.tsv",header=TRUE)
cells.keep <- rownames(subset(meta, UMI<2.5*sd(meta$UMI) & meta$mtProportion<0.1))
GM12878 <- GM12878[,paste0(cells.keep,"-1")]

n.cells <- rowSums(assay(GM12878, "counts")>0)
GM12878 <- GM12878[n.cells>=10,]

## Normalize
colData(GM12878)[, "sizeFactors"] <- scran::calculateSumFactors(as.matrix(assay(GM12878,"counts")))
GM12878 <- normalizeBySizeFactors(GM12878, inputAssay="counts", outputAssay="scran_normalized_counts", sizeFactors=GM12878$sizeFactors)

## Calculate dispersion
GM12878 <- calcMedian(GM12878,inputAssay="scran_normalized_counts",outputColumn="median")
GM12878 <- DMadjustedCV(GM12878,inputAssay="scran_normalized_counts")

## Save data
save(GM12878,file="../data/GM12878_adjusted_CV_data.RData")
