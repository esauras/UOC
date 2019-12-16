## Directorio de trabajo
wd <- "C:/Users/Esther/Desktop/Máster Bioinformática y Bioestadística/TFM/R"
setwd(wd)

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
barcodeKIRP <- sort(getSampleSummary(queryKIRP)$Tumor_Sample_Barcode)


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
barcodeKICH <- sort(getSampleSummary(queryKICH)$Tumor_Sample_Barcode)


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

## PTEN, TP53, TTN
KICHexp1 <- Samples.mRNASeq(format = "csv",
                            gene = c("TP53","PTEN","TTN"),
                            cohort = "KICH")
p1 <- ggplot(KICHexp1, aes(factor(gene), z.score))
p1 +
  geom_boxplot(aes(fill = factor(sample_type))) +
  # we drop some outlier, so plot looks nicer, this also causes the warning
  scale_y_continuous(limits = c(-1, 5)) +
  scale_fill_discrete(name = "Tissue")

## AFF3, AICDA, ZAN
KICHexp2 <- Samples.mRNASeq(format = "csv",
                            gene = c("ZAN","AFF3","AICDA"),
                            cohort = "KICH")
p2 <- ggplot(KICHexp2, aes(factor(gene), z.score))
p2 +
  geom_boxplot(aes(fill = factor(sample_type))) +
  # we drop some outlier, so plot looks nicer, this also causes the warning
  scale_y_continuous(limits = c(-1, 5)) +
  scale_fill_discrete(name = "Tissue")

## ARFGAP3, ARHGEF2
KICHexp3 <- Samples.mRNASeq(format = "csv",
                            gene = c("ARFGAP3","ARHGEF2"),
                            cohort = "KICH")
p3 <- ggplot(KICHexp3, aes(factor(gene), z.score))
p3 +
  geom_boxplot(aes(fill = factor(sample_type))) +
  # we drop some outlier, so plot looks nicer, this also causes the warning
  scale_y_continuous(limits = c(-1, 5)) +
  scale_fill_discrete(name = "Tissue")

## ATM, ATP10B
KICHexp4 <- Samples.mRNASeq(format = "csv",
                            gene = c("ATM","ATP10B"),
                            cohort = "KICH")
p4 <- ggplot(KICHexp4, aes(factor(gene), z.score))
p4 +
  geom_boxplot(aes(fill = factor(sample_type))) +
  # we drop some outlier, so plot looks nicer, this also causes the warning
  scale_y_continuous(limits = c(-1, 5)) +
  scale_fill_discrete(name = "Tissue")

## Método 2: TCGAANALYZE KIRP
queryKIRPexp <- TCGAbiolinks::GDCquery(project = "TCGA-KIRP", 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts",
                  legacy = FALSE)

samplesKIRP <- sort(getResults(queryKIRPexp, cols=c("cases")))

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
    list2 <- c(list2,paste0(y[1],"-",y[2],"-",y[3]))
  }
}

list <- c()
for (i in as.vector(barcodeKIRP)){
  x <- unlist(strsplit(i,"-"))
  list <- c(list,paste0(x[1],"-",x[2],"-",x[3]))
}

list <- unique(list)
list2 <- unique(list2)
list %in% list2
## Hay 281 muestras en list (barcodeKIRP) y 278 muestras en list2 (dataSmTP).
list[94] %in% list2 ## FALSE
list[155] %in% list2 ## FALSE
list[238] %in% list2 ## FALSE

library(SummarizedExperiment)
library(dplyr)
library(DT)

queryDown <- GDCquery(project = "TCGA-KIRP", 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - Counts", 
                      barcode = c(dataSmTP, dataSmNT))

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

dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmTP],
                            mat2 = dataFilt[,dataSmNT],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")  
## length 5301 
library(org.Hs.eg.db)
ann <- select(org.Hs.eg.db,keytype="ENSEMBL",keys=rownames(dataDEGs),columns=c("ENSEMBL","SYMBOL"))
ann <- ann %>% 
  group_by(ENSEMBL) %>% 
  dplyr::slice(1)
ann <- data.frame(ann)
dataDG <- as.data.frame(dataDEGs[order(rownames(dataDEGs)),])
resultKIRP <- cbind(ann,dataDG)

## Filtrar muestras
unique(KICHexp$tcga_participant_barcode)
KICHexpTP <- KICHexp[which(KICHexp$sample_type == "TP"),]
sort(KICHexpTP$tcga_participant_barcode)
KICHexpTP[which(KICHexpTP$tcga_participant_barcode == "TCGA-KL-8323"),]

## Método 2: TCGAANALYZE KICH
queryKICHexp <- TCGAbiolinks::GDCquery(project = "TCGA-KICH", 
                                       data.category = "Transcriptome Profiling",
                                       data.type = "Gene Expression Quantification",
                                       workflow.type = "HTSeq - Counts",
                                       legacy = FALSE)

samplesKICH <- sort(getResults(queryKICHexp, cols=c("cases")))

## Muestras tumor (TP)
dataSmTP2 <- TCGAquery_SampleTypes(barcode = samplesKICH,
                                  typesample = "TP")

## Muestras control (NP)
dataSmNT2 <- TCGAquery_SampleTypes(barcode = samplesKICH,
                                  typesample = "NT")
dataSmTP2_short <- dataSmTP2[1:10]
dataSmNT2_short <- dataSmNT2[1:10]

## Hay que seleccionar las muestras para las que se tenga tanto datos
## genómicos como datos de expresión. 
listt <- c()
for (i in as.vector(barcodeKICH)){
  x2 <- unlist(strsplit(i,"-"))
  listt <- c(listt,x2[3])
}

listt2 <- c()
for (j in dataSmTP2){
  y2 <- unlist(strsplit(j,"-"))
  if (y2[3] %in% listt){
    listt2 <- c(listt2,paste0(y2[1],"-",y2[2],"-",y2[3]))
  }
}

listt <- c()
for (i in as.vector(barcodeKICH)){
  x2 <- unlist(strsplit(i,"-"))
  listt <- c(listt,paste0(x2[1],"-",x2[2],"-",x2[3]))
}

listt <- unique(listt)
listt2 <- unique(listt2)
listt %in% listt2
## Hay 66 muestras en listt (barcodeKICH) y 65 muestras en listt2 (dataSmTP2).
listt[37] %in% listt2 ## FALSE

library(SummarizedExperiment)
library(dplyr)
library(DT)

queryDown2 <- GDCquery(project = "TCGA-KICH", 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - Counts", 
                      barcode = c(dataSmTP2, dataSmNT2))

GDCdownload(queryDown2)

dataPrep2 <- GDCprepare(queryDown2)
dataPrep2 <- TCGAanalyze_Preprocessing(object = dataPrep2, 
                                      cor.cut = 0.6,
                                      datatype = "HTSeq - Counts")                      

dataNorm2 <- TCGAanalyze_Normalization(tabDF = dataPrep2,
                                      geneInfo = geneInfoHT,
                                      method = "gcContent") 

boxplot(dataPrep2, outline = FALSE)
boxplot(dataNorm2, outline = FALSE)

dataFilt2 <- TCGAanalyze_Filtering(tabDF = dataNorm2,
                                  method = "quantile", 
                                  qnt.cut =  0.25)   

dataDEGs2 <- TCGAanalyze_DEA(mat1 = dataFilt2[,dataSmTP2],
                            mat2 = dataFilt2[,dataSmNT2],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")  
## length 5301 
library(org.Hs.eg.db)
ann2 <- select(org.Hs.eg.db,keytype="ENSEMBL",keys=rownames(dataDEGs2),columns=c("ENSEMBL","SYMBOL"))
ann2 <- ann2 %>% 
  group_by(ENSEMBL) %>% 
  dplyr::slice(1)
ann2 <- data.frame(ann2)
dataDG2 <- as.data.frame(dataDEGs2[order(rownames(dataDEGs2)),])
resultKICH <- cbind(ann2,dataDG2)

## ESTIMATE
library(utils)
## rforge <- "http://r-forge.r-project.org"
## install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
help(package="estimate")

## Datos clinicos
## if (!requireNamespace("BiocManager"))
##   install.packages("BiocManager")
## BiocManager::install("RTCGAToolbox")

## KIRP
library(RTCGAToolbox)
kirpData <- getFirehoseData(dataset="KIRP", runDate="20160128", clinical=TRUE)
clinicData <- getData(kirpData,"clinical")
head(clinicData)

smp <- row.names(clinicData)
smp <- gsub("[[:punct:]]","-",smp)
smp <- toupper(smp) ## Para convertirlos en mayusculas y que concuerden con KIRP
row.names(clinicData) <- smp
smp2 <- smp[smp %in% list2]
## Base de datos con las mismas muestras que en RNA-Seq (278).
newclinicData <- subset(clinicData, row.names(clinicData) %in% smp2)
newclinicData <- newclinicData[with(newclinicData, order(row.names(newclinicData))), ]

newclinicData <- newclinicData[, 3:5]
newclinicData[is.na(newclinicData[, 3]), 3] <- newclinicData[is.na(newclinicData[, 3]), 2]
survData <- data.frame(Samples=rownames(newclinicData),
                       Time=as.numeric(newclinicData[, 3]), Censor=as.numeric(newclinicData[, 1]))

getSurvival(dataObject=kirpData, geneSymbols=c("PIK3CA"),
            sampleTimeCensor=survData)

## KM curves FM1
library(survival)
row.names(newCOSM) <- newCOSM[,1]
g1s <- as.factor(row.names(newCOSM)[newCOSM$V1 <= summary(newCOSM$V1)[3]])
g2s <- as.factor(row.names(newCOSM)[newCOSM$V1 > summary(newCOSM$V1)[3]])
time.group1 <- as.numeric(survData[c(g1s,g2s),2])
censor.group1 <- as.numeric(survData[c(g1s,g2s),3])
surv.group1 <- c()
for (name in row.names(newCOSM)){
  if (name %in% g1s) surv.group1 <- c(surv.group1,1)
  if (name %in% g2s) surv.group1 <- c(surv.group1,2)
}
surv.fit1 <- survfit(Surv(time.group1,censor.group1)~as.factor(surv.group1))
surv.diff1 <- survdiff(Surv(time.group1,censor.group1)~as.factor(surv.group1))
par(mfrow=c(2,2))
plot(surv.fit1, xlab="Time", ylab="Survival", main="Mutational signature 1", col=c(2,4))
## No significativo p = 0.1

## KM curves FM6
h1s <- as.factor(row.names(newCOSM)[newCOSM$V6 <= summary(newCOSM$V6)[3]])
h2s <- as.factor(row.names(newCOSM)[newCOSM$V6 > summary(newCOSM$V6)[3]])
time.group6 <- as.numeric(survData[c(h1s,h2s),2])
censor.group6 <- as.numeric(survData[c(h1s,h2s),3])
surv.group6 <- c()
for (name in row.names(newCOSM)){
  if (name %in% h1s) surv.group6 <- c(surv.group6,1)
  if (name %in% h2s) surv.group6 <- c(surv.group6,2)
}
surv.fit6 <- survfit(Surv(time.group6,censor.group6)~as.factor(surv.group6))
surv.diff6 <- survdiff(Surv(time.group6,censor.group6)~as.factor(surv.group6))
plot(surv.fit6, xlab="Time", ylab="Survival", main="Mutational signature 6", col=c(2,4))
## No significativo p = 0.9

## KM curves FM19
i1s <- as.factor(row.names(newCOSM)[newCOSM$V19 <= summary(newCOSM$V19)[3]])
i2s <- as.factor(row.names(newCOSM)[newCOSM$V19 > summary(newCOSM$V19)[3]])
time.group19 <- as.numeric(survData[c(i1s,i2s),2])
censor.group19 <- as.numeric(survData[c(i1s,i2s),3])
surv.group19 <- c()
for (name in row.names(newCOSM)){
  if (name %in% i1s) surv.group19 <- c(surv.group19,1)
  if (name %in% i2s) surv.group19 <- c(surv.group19,2)
}
surv.fit19 <- survfit(Surv(time.group19,censor.group19)~as.factor(surv.group19))
surv.diff19 <- survdiff(Surv(time.group19,censor.group19)~as.factor(surv.group19))
plot(surv.fit19, xlab="Time", ylab="Survival", main="Mutational signature 19", col=c(2,4))
## No significativo p = 0.3

## KM curves FM22
j1s <- as.factor(row.names(newCOSM)[newCOSM$V22 <= summary(newCOSM$V22)[3]])
j2s <- as.factor(row.names(newCOSM)[newCOSM$V22 > summary(newCOSM$V22)[3]])
time.group22 <- as.numeric(survData[c(j1s,j2s),2])
censor.group22 <- as.numeric(survData[c(j1s,j2s),3])
surv.group22 <- c()
for (name in row.names(newCOSM)){
  if (name %in% j1s) surv.group22 <- c(surv.group22,1)
  if (name %in% j2s) surv.group22 <- c(surv.group22,2)
}
surv.fit22 <- survfit(Surv(time.group22,censor.group22)~as.factor(surv.group22))
surv.diff22 <- survdiff(Surv(time.group22,censor.group22)~as.factor(surv.group22))
plot(surv.fit22, xlab="Time", ylab="Survival", main="Mutational signature 22", col=c(2,4))
## No significativo p = 0.5

## KM curves Immune_score
dev.off()
row.names(newKIRP) <- newKIRP[,1]
k1s <- as.factor(row.names(newKIRP)[newKIRP$gr == 0])
k2s <- as.factor(row.names(newKIRP)[newKIRP$gr == 1])
time.groupIS <- as.numeric(survData[c(k1s,k2s),2])
censor.groupIS <- as.numeric(survData[c(k1s,k2s),3])
surv.groupIS <- c()
for (name in row.names(newKIRP)){
  if (name %in% k1s) surv.groupIS <- c(surv.groupIS,1)
  if (name %in% k2s) surv.groupIS <- c(surv.groupIS,2)
}
surv.fitIS <- survfit(Surv(time.groupIS,censor.groupIS)~as.factor(surv.groupIS))
surv.diffIS <- survdiff(Surv(time.groupIS,censor.groupIS)~as.factor(surv.groupIS))
plot(surv.fitIS, xlab="Time", ylab="Survival", main="Immune_score", col=c(2,4))
## Significativo p = 0.04

## Relación con sexo y edad
newclinicData <- subset(clinicData, row.names(clinicData) %in% smp2)
newclinicData <- newclinicData[with(newclinicData, order(row.names(newclinicData))), ]

summary(as.factor(newclinicData$gender))
tbl1 <- table(as.factor(newclinicData$gender),surv.group1) ## 0.5253
tbl1 <- table(as.factor(newclinicData$gender),surv.group6) ## 0.8041
tbl1 <- table(as.factor(newclinicData$gender),surv.group19) ## 0.6307
tbl1 <- table(as.factor(newclinicData$gender),surv.group22) ## 1
chisq.test(tbl1)
## No hay diferencias significativas en el sexo entre grupos.

summary(as.numeric(newclinicData$years_to_birth))
l1s <- as.factor(row.names(newclinicData)[as.numeric(newclinicData$years_to_birth) <= summary(as.numeric(newclinicData$years_to_birth))[3]])
l2s <- as.factor(row.names(newclinicData)[as.numeric(newclinicData$years_to_birth) > summary(as.numeric(newclinicData$years_to_birth))[3]])
l3s <- as.factor(row.names(newclinicData)[is.na(newclinicData$years_to_birth) == TRUE])
surv.groupAGE <- c()
for (name in row.names(newclinicData)){
  if (name %in% l1s) surv.groupAGE <- c(surv.groupAGE,1)
  if (name %in% l2s) surv.groupAGE <- c(surv.groupAGE,2)
  if (name %in% l3s) surv.groupAGE <- c(surv.groupAGE,0)
}
tab <- newclinicData[,c(1:2)]
dhf <- cbind(tab,surv.groupAGE,newKIRP$V19)
dhf <- dhf[dhf$surv.groupAGE != 0,]
boxplot(dhf$`newKIRP$V19`~dhf$surv.groupAGE,xlab="Age",ylab="Level of expression",
        main="Mutational signature 19",col=c("cadetblue","coral1"))
tab <- cbind(tab,surv.groupAGE,surv.group1,surv.group6,surv.group19,surv.group22)
tab <- tab[tab$surv.groupAGE != 0,]
tbl2 <- table(as.factor(tab$surv.groupAGE),tab$surv.group1) ## 0.7934
tbl2 <- table(as.factor(tab$surv.groupAGE),tab$surv.group6) ## 0.3617
tbl2 <- table(as.factor(tab$surv.groupAGE),tab$surv.group19) ## 0.01793**
tbl2 <- table(as.factor(tab$surv.groupAGE),tab$surv.group22) ## 0.7704
chisq.test(tbl2)

## Immune_score
setwd("C:/Users/Esther/Desktop/Máster Bioinformática y Bioestadística/TFM/ESTIMATE")
KIRPscore <- read.csv("kidney_renal_papillary_cell_carcinoma_RNAseqV2.txt",sep="\t")
KIRPscore

## Seleccionar las muestras que coinciden (list2)
idfilt <- c()
for (i in KIRPscore$ID){
  z <- unlist(strsplit(i,"-"))
  idfilt <- c(idfilt,paste0(z[1],"-",z[2],"-",z[3]))
}
idfilt2 <- unique(idfilt[idfilt %in% list2])

KIRPscore$ID <- idfilt
nums <- which(KIRPscore$ID %in% idfilt2)
options(max.print=999999)

## Base de datos IMMUNESCORE depurada
newKIRPscore <- KIRPscore[nums,]
newKIRPscore <- read.csv("KIRPscore.csv") ## Abrir directamente
newKIRPscore <- newKIRPscore[,c(2:5)]

## De aquí seleccionar dos grupos (<0 y >0)
scorepos <- which(newKIRPscore$Immune_score >0) ## 150
scoreneg <- which(newKIRPscore$Immune_score <0) ## 128

## COSMIC
setwd(wd)
COSM <- read.csv("COSMICsignatures.csv", sep=";")
COSM <- t(COSM)
COSM <- COSM[-(1:2),]

nfilt <- c()
for (i in rownames(COSM)){
  a <- unlist(strsplit(i,"[[:punct:]]"))
  nfilt <- c(nfilt,paste0(a[1],"-",a[2],"-",a[3]))
}
nfilt2 <- nfilt[nfilt %in% list2]

rownames(COSM) <- nfilt
num <- which(rownames(COSM) %in% nfilt2)

## Base de datos COSMIC depurada
newCOSM <- COSM[num,]
newCOSM <- read.csv("COSM.csv") ## Abrir directamente

## Correlación
newKIRP <- cbind(newCOSM,newKIRPscore$Immune_score)
colnames(newKIRP)[32] <- "Immune_score"
hist(newKIRP$Immune_score,xlab="Immune_score")
shapiro.test(newKIRP$Immune_score) ## No normalidad - pruebas no paramétricas

gr <- c()
for (i in newKIRP$Immune_score){
  if (i > median(newKIRP$Immune_score)) gr <- c(gr,1)
  else gr <- c(gr,0)
}
newKIRP <- cbind(newKIRP,gr)
gr2 <- c()
for (i in newKIRP$Immune_score){
  if ((quantile(newKIRP$Immune_score)[1] <= i) && (i < quantile(newKIRP$Immune_score)[2])) 
    gr2 <- c(gr2,0)
  if ((quantile(newKIRP$Immune_score)[2] <= i) && (i < quantile(newKIRP$Immune_score)[3])) 
    gr2 <- c(gr2,1)
  if ((quantile(newKIRP$Immune_score)[3] <= i) && (i < quantile(newKIRP$Immune_score)[4])) 
    gr2 <- c(gr2,2)
  if ((quantile(newKIRP$Immune_score)[4] <= i) && (i <= quantile(newKIRP$Immune_score)[5])) 
    gr2 <- c(gr2,3)
}
newKIRP <- cbind(newKIRP,gr2)
par(mfrow=c(2,2))
wilcox.test(newKIRP$V1~newKIRP$gr) ## 0.003298 **
boxplot(newKIRP$V1~newKIRP$gr,xlab="Immune_score",ylab="Level of expression",
        main="Mutational signature 1",col=c("cadetblue","coral1"))
boxplot(newKIRP$V1~newKIRP$gr2,xlab="Immune_score",ylab="Level of expression",
        main="Mutational signature 1",col=c("cadetblue","cadetblue2","coral","coral3"))
wilcox.test(newKIRP$V2~newKIRP$gr) ## 0.2632
wilcox.test(newKIRP$V3~newKIRP$gr) ## 0.3898
wilcox.test(newKIRP$V4~newKIRP$gr) ## 0.2357
wilcox.test(newKIRP$V5~newKIRP$gr) ## 0.6751
wilcox.test(newKIRP$V6~newKIRP$gr) ## 0.007948 **
boxplot(newKIRP$V6~newKIRP$gr,xlab="Immune_score",ylab="Level of expression",
        main="Mutational signature 6",col=c("cadetblue","coral1"))
boxplot(newKIRP$V6~newKIRP$gr2,xlab="Immune_score",ylab="Level of expression",
        main="Mutational signature 6",col=c("cadetblue","cadetblue2","coral","coral3"))
wilcox.test(newKIRP$V7~newKIRP$gr) ## 0.9141
wilcox.test(newKIRP$V8~newKIRP$gr) ## 0.8915
wilcox.test(newKIRP$V9~newKIRP$gr) ## 0.4719
wilcox.test(newKIRP$V10~newKIRP$gr) ## 0.6435
wilcox.test(newKIRP$V11~newKIRP$gr) ## 0.7505
wilcox.test(newKIRP$V12~newKIRP$gr) ## 0.2068
wilcox.test(newKIRP$V13~newKIRP$gr) ## 0.4608
wilcox.test(newKIRP$V14~newKIRP$gr) ## 0.2823
wilcox.test(newKIRP$V15~newKIRP$gr) ## 0.9488
wilcox.test(newKIRP$V16~newKIRP$gr) ## 0.785
wilcox.test(newKIRP$V17~newKIRP$gr) ## 0.7278
wilcox.test(newKIRP$V18~newKIRP$gr) ## 0.3127
wilcox.test(newKIRP$V19~newKIRP$gr) ## 0.002028 **
boxplot(newKIRP$V19~newKIRP$gr,xlab="Immune_score",ylab="Level of expression",
        main="Mutational signature 19",col=c("cadetblue","coral1"))
boxplot(newKIRP$V19~newKIRP$gr2,xlab="Immune_score",ylab="Level of expression",
        main="Mutational signature 19",col=c("cadetblue","cadetblue2","coral","coral3"))
wilcox.test(newKIRP$V20~newKIRP$gr) ## 0.5143
wilcox.test(newKIRP$V21~newKIRP$gr) ## 0.405
wilcox.test(newKIRP$V22~newKIRP$gr) ## 0.03125 **
boxplot(newKIRP$V22~newKIRP$gr,xlab="Immune_score",ylab="Level of expression",
        main="Mutational signature 22",col=c("cadetblue","coral1"))
boxplot(newKIRP$V22~newKIRP$gr2,xlab="Immune_score",ylab="Level of expression",
        main="Mutational signature 22",col=c("cadetblue","cadetblue2","coral","coral3"))
wilcox.test(newKIRP$V23~newKIRP$gr) ## 0.2644
wilcox.test(newKIRP$V24~newKIRP$gr) ## 0.6549
wilcox.test(newKIRP$V25~newKIRP$gr) ## 0.2372
wilcox.test(newKIRP$V26~newKIRP$gr) ## 0.2789
wilcox.test(newKIRP$V27~newKIRP$gr) ## 0.319
wilcox.test(newKIRP$V28~newKIRP$gr) ## 0.1283
wilcox.test(newKIRP$V29~newKIRP$gr) ## 0.8259
wilcox.test(newKIRP$V30~newKIRP$gr) ## 0.9301

## KICH
library(RTCGAToolbox)
kichData <- getFirehoseData(dataset="KICH", runDate="20160128", clinical=TRUE)
clinickichData <- getData(kichData,"clinical")
head(clinickichData)

smpp <- row.names(clinickichData)
smpp <- gsub("[[:punct:]]","-",smpp)
smpp <- toupper(smpp) ## Para convertirlos en mayusculas y que concuerden con KIRP
row.names(clinickichData) <- smpp
smpp2 <- smpp[smpp %in% listt2]
## Base de datos con las mismas muestras que en RNA-Seq (78).
newclinickichData <- subset(clinickichData, row.names(clinickichData) %in% smpp2)
newclinickichData <- newclinickichData[with(newclinickichData, order(row.names(newclinickichData))), ]

newclinickichData <- newclinickichData[, 3:5]
newclinickichData[is.na(newclinickichData[, 3]), 3] <- newclinickichData[is.na(newclinickichData[, 3]), 2]
survkichData <- data.frame(Samples=rownames(newclinickichData),
                       Time=as.numeric(newclinickichData[, 3]), Censor=as.numeric(newclinickichData[, 1]))

getSurvival(dataObject=kichData, geneSymbols=c("PIK3CA"),
            sampleTimeCensor=survkichData)

## Immune_score
setwd("C:/Users/Esther/Desktop/Máster Bioinformática y Bioestadística/TFM/ESTIMATE")
KICHscore <- read.csv("kidney_chromophobe_renal_cell_carcinoma_RNAseqV2.txt",sep="\t")
KICHscore

## Seleccionar las muestras que coinciden (listt2)
idfiltt <- c()
for (i in KICHscore$ID){
  z2 <- unlist(strsplit(i,"-"))
  idfiltt <- c(idfiltt,paste0(z2[1],"-",z2[2],"-",z2[3]))
}
idfiltt2 <- unique(idfiltt[idfiltt %in% listt2])

KICHscore$ID <- idfiltt
numss <- which(KICHscore$ID %in% idfiltt2)
options(max.print=999999)

## Base de datos IMMUNESCORE depurada
newKICHscore <- KICHscore[numss,]
newKICHscore <- read.csv("KICHscore.csv") ## Abrir directamente
newKICHscore <- newKICHscore[,c(2:5)]

## De aquí seleccionar dos grupos (<0 y >0)
scorepos2 <- which(newKICHscore$Immune_score >0) ## 12
scoreneg2 <- which(newKICHscore$Immune_score <0) ## 53

## COSMIC
setwd(wd)
COSMKICH <- read.csv("COSMICKICHsignatures.csv", sep=";")
COSMKICH <- t(COSMKICH)
COSMKICH <- COSMKICH[-(1:2),]

nfiltt <- c()
for (i in rownames(COSMKICH)){
  a2 <- unlist(strsplit(i,"[[:punct:]]"))
  nfiltt <- c(nfiltt,paste0(a2[1],"-",a2[2],"-",a2[3]))
}
nfiltt2 <- nfiltt[nfiltt %in% listt2]

rownames(COSMKICH) <- nfiltt
numm <- which(rownames(COSMKICH) %in% nfiltt2)

## Base de datos COSMIC depurada
newCOSMKICH <- COSMKICH[numm,]
newCOSMKICH <- read.csv("COSMKICH.csv") ## Abrir directamente

## Correlación
newKICH <- cbind(newCOSMKICH,newKICHscore$Immune_score)
colnames(newKICH)[32] <- "Immune_score"
hist(newKICH$Immune_score,xlab="Immune_score")
shapiro.test(newKICH$Immune_score) ## No normalidad - pruebas no paramétricas

grr <- c()
for (i in newKICH$Immune_score){
  if (i > median(newKICH$Immune_score)) grr <- c(grr,1)
  else grr <- c(grr,0)
}
newKICH <- cbind(newKICH,grr)
grr2 <- c()
for (i in newKICH$Immune_score){
  if ((quantile(newKICH$Immune_score)[1] <= i) && (i < quantile(newKICH$Immune_score)[2])) 
    grr2 <- c(grr2,0)
  if ((quantile(newKICH$Immune_score)[2] <= i) && (i < quantile(newKICH$Immune_score)[3])) 
    grr2 <- c(grr2,1)
  if ((quantile(newKICH$Immune_score)[3] <= i) && (i < quantile(newKICH$Immune_score)[4])) 
    grr2 <- c(grr2,2)
  if ((quantile(newKICH$Immune_score)[4] <= i) && (i <= quantile(newKICH$Immune_score)[5])) 
    grr2 <- c(grr2,3)
}
newKICH <- cbind(newKICH,grr2)
par(mfrow=c(2,2))
wilcox.test(newKICH$V1~newKICH$grr) ## 0.6028
wilcox.test(newKICH$V2~newKICH$grr) ## 0.878
wilcox.test(newKICH$V3~newKICH$grr) ## 0.8418
wilcox.test(newKICH$V4~newKICH$grr) ## 0.1972
wilcox.test(newKICH$V5~newKICH$grr) ## 0.1671
wilcox.test(newKICH$V6~newKICH$grr) ## 0.9337
wilcox.test(newKICH$V7~newKICH$grr) ## 0.6244
wilcox.test(newKICH$V8~newKICH$grr) ## 0.4197
wilcox.test(newKICH$V9~newKICH$grr) ## 0.3673
wilcox.test(newKICH$V10~newKICH$grr) ## 0.09989
wilcox.test(newKICH$V11~newKICH$grr) ## 0.421
wilcox.test(newKICH$V12~newKICH$grr) ## 0.2771
wilcox.test(newKICH$V13~newKICH$grr) ## 0.6349
wilcox.test(newKICH$V14~newKICH$grr) ## 0.7014
wilcox.test(newKICH$V15~newKICH$grr) ## 0.3909
wilcox.test(newKICH$V16~newKICH$grr) ## 0.9064
wilcox.test(newKICH$V17~newKICH$grr) ## 0.9417
wilcox.test(newKICH$V18~newKICH$grr) ## 0.8144
wilcox.test(newKICH$V19~newKICH$grr) ## 0.4276
wilcox.test(newKICH$V20~newKICH$grr) ## 0.678
wilcox.test(newKICH$V21~newKICH$grr) ## 0.7357
wilcox.test(newKICH$V22~newKICH$grr) ## 0.6858
wilcox.test(newKICH$V23~newKICH$grr) ## 0.4045
wilcox.test(newKICH$V24~newKICH$grr) ## 0.794
wilcox.test(newKICH$V25~newKICH$grr) ## 0.5701
wilcox.test(newKICH$V26~newKICH$grr) ## 0.08648
wilcox.test(newKICH$V27~newKICH$grr) ## 0.7495
wilcox.test(newKICH$V28~newKICH$grr) ## 0.2331
wilcox.test(newKICH$V29~newKICH$grr) ## 0.5797
wilcox.test(newKICH$V30~newKICH$grr) ## 0.7512

## KM curves Immune_score
row.names(newKICH) <- newKICH[,1]
m1s <- as.factor(row.names(newKICH)[newKICH$grr == 0])
m2s <- as.factor(row.names(newKICH)[newKICH$grr == 1])
time.groupkichIS <- as.numeric(survData[c(m1s,m2s),2])
censor.groupkichIS <- as.numeric(survData[c(m1s,m2s),3])
surv.groupkichIS <- c()
for (name in row.names(newKICH)){
  if (name %in% m1s) surv.groupkichIS <- c(surv.groupkichIS,1)
  if (name %in% m2s) surv.groupkichIS <- c(surv.groupkichIS,2)
}
surv.fitkichIS <- survfit(Surv(time.groupkichIS,censor.groupkichIS)~as.factor(surv.groupkichIS))
surv.diffkichIS <- survdiff(Surv(time.groupkichIS,censor.groupkichIS)~as.factor(surv.groupkichIS))
plot(surv.fitkichIS, xlab="Time", ylab="Survival", main="Immune_score", col=c(2,4))
## No significativo p = 0.5