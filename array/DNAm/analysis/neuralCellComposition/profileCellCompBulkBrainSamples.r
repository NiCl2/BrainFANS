##---------------------------------------------------------------------#
##
## Title: Profile cell composition in bulk brain samples
##
## Purpose of script: To calculate cellular composition of bulk brain DNAm
## profiles and characterise against biological factors
##
## Author: Eilis Hannon
##
## Date Created: 13/12/2022
##
## Adapted 11/05/2023 by Nicholas Clifton
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# takes bulk brain beta matrix and calculates cellular composition
# data provided as R object with betas matrix and pheno data.frame


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  setwd("~/OneDrive - University of Exeter/Documents/methylation/Setd1a_array/")
  args[1]<-"Setd1a_mouseArray_FilteredNormalised_Betas.rdat"
  args[2]<-"../RawData_Mouse_DNAm/MouseArray_CellDeconv_FilteredNormalised_Betas_CoefBrainModel.rdata"
  args[3]<-"."
} 
bulkPath <- args[1]
modelPath <- args[2]
plotPath <- args[3]

# setwd("~/OneDrive - University of Exeter/Documents/methylation/Setd1a_array/")
# bulkPath <- "Setd1a_mouseArray_FilteredNormalised_Betas.rdat"
# modelPath <- "../RawData_Mouse_DNAm/MouseArray_CellDeconv_FilteredNormalised_Betas_CoefBrainModel.rdata"
# plotPath <- "."


#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(CETYGO)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(corrplot)

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

load(modelPath)
load(bulkPath)

betas <- normbetas
rm(normbetas)
pheno <- read.csv("Setd1a_array_meta.csv")
pheno$basename <- paste(pheno$Chip.ID, pheno$Chip.Location, sep = "_")
pos <- position_dodge(0.9)

#----------------------------------------------------------------------#
# CALCULATE CELL COMPOSITION
#----------------------------------------------------------------------#

sites<-intersect(rownames(betas), rownames(brainCoefANOVA))
predCCANOVA<-projectCellTypeWithError(betas[sites,], brainCoefANOVA[sites,])

save(predCCANOVA, file = paste(tools::file_path_sans_ext(bulkPath), "CellCompEstimates.Rdata", sep = "_"))

## reformat into single data.frame for plotting
sumOut <- cbind(pheno, predCCANOVA[pheno$basename, ])
sumOut$age <- factor(sumOut$age, levels = c("E14", "E18", "P7", "P35", "P70"))
sumOut$genotype <- factor(sumOut$genotype, levels = c("WT", "Het"))

#----------------------------------------------------------------------#
# TEST CETYGO AGAINST GENOTYPE & AGE
#----------------------------------------------------------------------#

model<-lm(CETYGO ~ age + genotype, data = sumOut)

summary(model)$coefficients

#----------------------------------------------------------------------#
# PLOT CETYGO AGAINST AGE
#----------------------------------------------------------------------#

ggplot(sumOut, aes(x=age, y=CETYGO, fill = genotype))  +
  geom_violin(position = pos, scale = 'width')  +
  stat_summary(fun = "mean", 
               geom = "point", 
               position = pos) + geom_vline(xintercept = c(1:length(unique(sumOut$age))), linetype="dotted") +
  ylim(c(0, max(sumOut$CETYGO))) + scale_fill_manual(values = c("blue", "red"))
  
ggsave(filename = file.path(plotPath, "ViolinPlotCETYGOxAge.pdf"),  units = "in", width = 12, height = 8)


#----------------------------------------------------------------------#
# TEST AGAINST AGE & GENOTYPE
#----------------------------------------------------------------------#

model<-lm(NeuN ~ age + genotype, data = sumOut)

summary(model)$coefficients

#----------------------------------------------------------------------#
# PLOT DISTRIBUTION OF PREDICTED CELLULAR COMPOSITION
#----------------------------------------------------------------------#
 
f1 <- ggplot(sumOut, aes(x=age, y=NeuN, fill = genotype))  +
		  geom_violin(position = pos, scale = 'width')  +
		  stat_summary(fun = "mean", 
					   geom = "point", 
					   position = pos) + scale_fill_manual(values = c("blue", "red"))

f2 <- ggplot(sumOut, aes(x=age, y=NonNeuronal, fill = genotype))  +
  geom_violin(position = pos, scale = 'width')  +
  stat_summary(fun = "mean", 
               geom = "point", 
               position = pos) + scale_fill_manual(values = c("blue", "red"))

ggarrange(f1, f2, nrow = 1, ncol = 2)

ggsave(filename = file.path(plotPath, paste0("ViolinPlotPredPropnSamples.pdf")),  units = "in", width = 18, height = 8)


