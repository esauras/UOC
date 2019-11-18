## Directorio de trabajo
wd <- setwd("C:/Users/Esther/Desktop/Máster Bioinformática y Bioestadística/TFM/R")

## Instalación de TCGAbiolinks desde biocLite
## source("https://bioconductor.org/biocLite.R")
## biocLite("TCGAbiolinks", dependecies=TRUE)
## if (!requireNamespace("BiocManager", quietly = TRUE))
##   install.packages("BiocManager")
## BiocManager::install("biocLite")

## Librerías
require(TCGAbiolinks)
require(dplyr)
require(maftools)
require(data.table)


## Descarga de los datos genómicos para carcinoma renal de células cromófobas
## (cRCC) y carcinoma renal papilar (pRCC).
cRCC <- TCGAbiolinks::GDCquery_Maf("KICH", directory=wd, pipelines = "muse")
pRCC <- TCGAbiolinks::GDCquery_Maf("KIRP", directory=wd, pipelines = "muse")


## KIRP
queryKIRP <- GDCquery_Maf("KIRP", pipelines = "muse") %>% read.maf

data.table(getSampleSummary(queryKIRP),
           filter = 'top',
           options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
           rownames = FALSE)

plotmafSummary(maf = queryKIRP, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()

oncoplot(maf = queryKIRP, top = 10, removeNonMutated = TRUE)
titv <- titv(maf = queryKIRP, plot = FALSE, useSyn = TRUE)
## plot titv summary
plotTiTv(res = titv)

## 281 muestras
barcodeKIRP <- getSampleSummary(queryKIRP)$Tumor_Sample_Barcode


## KICH
queryKICH <- GDCquery_Maf("KICH", pipelines = "muse") %>% read.maf

data.table(getSampleSummary(queryKICH),
           filter = 'top',
           options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
           rownames = FALSE)

plotmafSummary(maf = queryKICH, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()

oncoplot(maf = queryKICH, top = 10, removeNonMutated = TRUE)
titv <- titv(maf = queryKICH, plot = FALSE, useSyn = TRUE)
## plot titv summary
plotTiTv(res = titv)

## 66 muestras
barcodeKICH <- getSampleSummary(queryKICH)$Tumor_Sample_Barcode


## Instalación de los paquetes necesarios antes de utilizar MuSiCa.
if(!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
if(!require(MutationalPatterns)) {BiocManager::install("MutationalPatterns")}
if(!require(VariantAnnotation)) {BiocManager::install("VariantAnnotation")}
if(!require(BSgenome.Hsapiens.UCSC.hg38)) {BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")}
if(!require(BSgenome.Hsapiens.UCSC.hg19)) {BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")}
if(!require(BSgenome.Hsapiens.1000genomes.hs37d5)) {BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")}

if(!require(ggplot2)) install.packages("ggplot2")
if(!require(heatmaply)) install.packages("heatmaply")
if(!require(gplots)) install.packages("gplots")
if(!require(reshape2)) install.packages("reshape2")
if(!require(data.table)) install.packages("data.table")
if(!require(readxl)) install.packages("readxl")
if(!require(openxlsx)) install.packages("openxlsx")

if(!require(shiny)) install.packages("shiny")
if(!require(shinyBS)) install.packages("shinyBS")
if(!require(devtools)) install.packages("devtools")
if(!require(shinysky)) {library(devtools); devtools::install_github("AnalytixWare/ShinySky")}
if(!require(shinyjs)) install.packages("shinyjs")
if(!require(V8)) install.packages("V8")
if(!require(shinythemes)) install.packages("shinythemes")
if(!require(plotly)) install.packages("plotly")
if(!require(webshot)) install.packages("webshot")
webshot::install_phantomjs()

library(shiny)
runGitHub("musica", "marcos-diazg")


## Descarga de los datos de expresión para carcinoma renal de células 
## cromófobas (cRCC) y carcinoma renal papilar (pRCC).
## Método 1: FIREBROWSER
## install.packages("FirebrowseR")
## devtools::install_github("mariodeng/FirebrowseR")
require(FirebrowseR)
require(ggplot2)
genes <- c("TP53","PTEN","TTN","ZAN","AFF3","AICDA","ARFGAP3","ARHGEF2",
           "ATM","ATP10B")
par(mfrow=c(1,3)) 

## TP53
KICHexp1 <- Samples.mRNASeq(format = "csv",
                            gene = genes[1],
                            cohort = "KICH")
p1 <- ggplot(KICHexp1, aes(factor(gene), z.score))
p1 +
  geom_boxplot(aes(fill = factor(sample_type))) +
  # we drop some outlier, so plot looks nicer, this also causes the warning
  scale_y_continuous(limits = c(-1, 5)) +
  scale_fill_discrete(name = "Tissue")

## PTEN
KICHexp2 <- Samples.mRNASeq(format = "csv",
                            gene = genes[2],
                            cohort = "KICH")
p2 <- ggplot(KICHexp2, aes(factor(gene), z.score))
p2 +
  geom_boxplot(aes(fill = factor(sample_type))) +
  # we drop some outlier, so plot looks nicer, this also causes the warning
  scale_y_continuous(limits = c(-1, 5)) +
  scale_fill_discrete(name = "Tissue")

## TTN
KICHexp3 <- Samples.mRNASeq(format = "csv",
                            gene = genes[3],
                            cohort = "KICH")
p3 <- ggplot(KICHexp3, aes(factor(gene), z.score))
p3 +
  geom_boxplot(aes(fill = factor(sample_type))) +
  # we drop some outlier, so plot looks nicer, this also causes the warning
  scale_y_continuous(limits = c(-1, 5)) +
  scale_fill_discrete(name = "Tissue")

## ZAN
KICHexp4 <- Samples.mRNASeq(format = "csv",
                            gene = genes[4],
                            cohort = "KICH")
p4 <- ggplot(KICHexp4, aes(factor(gene), z.score))
p4 +
  geom_boxplot(aes(fill = factor(sample_type))) +
  # we drop some outlier, so plot looks nicer, this also causes the warning
  scale_y_continuous(limits = c(-1, 5)) +
  scale_fill_discrete(name = "Tissue")

## AFF3
KICHexp5 <- Samples.mRNASeq(format = "csv",
                            gene = genes[5],
                            cohort = "KICH")
p5 <- ggplot(KICHexp5, aes(factor(gene), z.score))
p5 +
  geom_boxplot(aes(fill = factor(sample_type))) +
  # we drop some outlier, so plot looks nicer, this also causes the warning
  scale_y_continuous(limits = c(-1, 5)) +
  scale_fill_discrete(name = "Tissue")

## AICDA
KICHexp6 <- Samples.mRNASeq(format = "csv",
                            gene = genes[6],
                            cohort = "KICH")
p6 <- ggplot(KICHexp6, aes(factor(gene), z.score))
p6 +
  geom_boxplot(aes(fill = factor(sample_type))) +
  # we drop some outlier, so plot looks nicer, this also causes the warning
  scale_y_continuous(limits = c(-1, 5)) +
  scale_fill_discrete(name = "Tissue")

## ARFGAP3
KICHexp7 <- Samples.mRNASeq(format = "csv",
                            gene = genes[7],
                            cohort = "KICH")
p7 <- ggplot(KICHexp7, aes(factor(gene), z.score))
p7 +
  geom_boxplot(aes(fill = factor(sample_type))) +
  # we drop some outlier, so plot looks nicer, this also causes the warning
  scale_y_continuous(limits = c(-1, 5)) +
  scale_fill_discrete(name = "Tissue")

## ARHGEF2
KICHexp8 <- Samples.mRNASeq(format = "csv",
                            gene = genes[8],
                            cohort = "KICH")
p8 <- ggplot(KICHexp8, aes(factor(gene), z.score))
p8 +
  geom_boxplot(aes(fill = factor(sample_type))) +
  # we drop some outlier, so plot looks nicer, this also causes the warning
  scale_y_continuous(limits = c(-1, 5)) +
  scale_fill_discrete(name = "Tissue")

## ATM
KICHexp9 <- Samples.mRNASeq(format = "csv",
                            gene = genes[9],
                            cohort = "KICH")
p9 <- ggplot(KICHexp9, aes(factor(gene), z.score))
p9 +
  geom_boxplot(aes(fill = factor(sample_type))) +
  # we drop some outlier, so plot looks nicer, this also causes the warning
  scale_y_continuous(limits = c(-1, 5)) +
  scale_fill_discrete(name = "Tissue")

## ATP10B
KICHexp10 <- Samples.mRNASeq(format = "csv",
                             gene = genes[10],
                             cohort = "KICH")
p10 <- ggplot(KICHexp10, aes(factor(gene), z.score))
p10 +
  geom_boxplot(aes(fill = factor(sample_type))) +
  # we drop some outlier, so plot looks nicer, this also causes the warning
  scale_y_continuous(limits = c(-1, 5)) +
  scale_fill_discrete(name = "Tissue")

## Método 2: TCGAANALYZE
queryKIRPexp <- TCGAbiolinks::GDCquery(project = "TCGA-KIRP", 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts",
                  legacy = FALSE)

samplesKIRP <- getResults(queryKIRPexp, cols=c("cases"))

## Muestras tumor (TP)
dataSmTP <- TCGAquery_SampleTypes(barcode = samplesKIRP,
                                  typesample = "TP")

## Muestras control (NP)
dataSmNT <- TCGAquery_SampleTypes(barcode = samplesKIRP,
                                  typesample = "NT")
dataSmTP_short <- dataSmTP[1:10]
dataSmNT_short <- dataSmNT[1:10]

## Hay que seleccionar las muestras para las que se tenga tanto datos
## genómicos como datos de expresión. 
list <- c()
for (i in as.vector(barcodeKIRP)){
  x <- unlist(strsplit(i,"-"))
  list <- c(list,x[3])
}

list2 <- c()
for (j in dataSmTP){
  y <- unlist(strsplit(j,"-"))
  if (y[3] %in% list){
    list2 <- c(list2,y[3])
  }
}
list <- unique(list)
list2 <- unique(list2)
list %in% list2
## Hay 281 muestras en list (barcodeKIRP) y 278 muestras en list2 (dataSmTP).
list[38] %in% list2 ## FALSE
list[75] %in% list2 ## FALSE
list[253] %in% list2 ## FALSE

library(SummarizedExperiment)
library(dplyr)
library(DT)

queryDown <- GDCquery(project = "TCGA-KIRP", 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - Counts", 
                      barcode = c(dataSmTP_short, dataSmNT_short))

GDCdownload(queryDown)

dataPrep <- GDCprepare(queryDown)
dataPrep <- TCGAanalyze_Preprocessing(object = dataPrep, 
                                      cor.cut = 0.6,
                                      datatype = "HTSeq - Counts")                      

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfoHT,
                                      method = "gcContent") 

boxplot(dataPrep, outline = FALSE)
boxplot(dataNorm, outline = FALSE)

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)   

dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmTP_short],
                            mat2 = dataFilt[,dataSmNT_short],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")  

## Filtrar muestras
unique(KICHexp$tcga_participant_barcode)
KICHexpTP <- KICHexp[which(KICHexp$sample_type == "TP"),]
sort(KICHexpTP$tcga_participant_barcode)
KICHexpTP[which(KICHexpTP$tcga_participant_barcode == "TCGA-KL-8323"),]

## ESTIMATE
library(utils)
## rforge <- "http://r-forge.r-project.org"
## install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
help(package="estimate")
