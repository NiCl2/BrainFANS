##---------------------------------------------------------------------#
##
## Title: Derive models for cell composition in brain
##
## Purpose of script: To train a series of models to predict different combinations of neural cell types
##
## Author: Eilis Hannon
## Adapted: Nicholas Clifton
##
## Date Created: 19/10/2022
## Date Adapted: 10/05/2023
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# identifies CpGs for brain deconvolution models
# uses IDOL method & ANOVA (pickCompProbes)


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

set.seed(100)
numProbes<-100 # number of sites included in model
probeSelect<- "any" 

# args<-commandArgs(trailingOnly = TRUE)
# normData<-args[1]
# refPath<-args[2]
# refPanelPath<-args[3]

setwd("~/OneDrive - University of Exeter/Documents/methylation/RawData_Mouse_DNAm/")
normData <- "MouseArray_CellDeconv_FilteredNormalised_Betas.rdat"
refPath <- "../MouseMethylation-12v1-0_A2.csv"
refPanelPath <- ""


#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#
library(dplyr)
# devtools::install_github("ds420/CETYGO")
library(CETYGO)
# devtools::install_github("immunomethylomics/IDOL")
# library(IDOL)
# library(doParallel)
# nworkers <- 2 ## if pushed too high causes OOM errors


#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#
message("Loading data")
load(normData)
pheno.all <- sample_sheet
norm.all <- normbetas
rm(sample_sheet, normbetas)

probeAnnot <- read.csv("../MouseMethylation-12v1-0_A2.csv", skip = 7) %>%
  mutate(probeID = paste(Name, AddressA_ID, sep = "_")) %>%
  filter(CHR != "") %>%
  dplyr::select(probeID, Name, CHR, MAPINFO) 
rownames(probeAnnot) <- probeAnnot$probeID

colnames(pheno.all)[colnames(pheno.all) == "Nuclei_Fraction"] <- "CellType" # need to rename for use with IDOL functions
pheno.all$CellType[pheno.all$CellType != "NeuN"] <- "NonNeuronal"
pheno.all$CellType<-factor(pheno.all$CellType)

probeAnnot<-probeAnnot[rownames(norm.all),]

## filter to autosome probes
message("Filtering data")
norm.all<-norm.all[!probeAnnot$CHR %in% c("X", "Y"),]

#----------------------------------------------------------------------#
# SELECT PROBES FOR DECONVOLUTION
#----------------------------------------------------------------------#

indexCells <- split(1:nrow(pheno.all), pheno.all$CellType)
exclude<-pheno.all$Basename[unlist(lapply(indexCells, sample, size = 1))]
names(exclude)<-sort(unique(pheno.all$CellType))

cellTypes <- unique(pheno.all$CellType)
cellInd<-pheno.all$CellType

## use ANOVA
compData <- pickCompProbesMatrix(rawbetas=norm.all, cellInd=cellInd, cellTypes = cellTypes, numProbes = numProbes/2, probeSelect = probeSelect)
brainCoefANOVA<-compData$coefEsts
save(brainCoefANOVA, file = paste(tools::file_path_sans_ext(normData), "CoefBrainModel.rdata", sep = "_"))


