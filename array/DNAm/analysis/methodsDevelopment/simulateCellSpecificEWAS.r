##---------------------------------------------------------------------#
##
## Title: Simulate Cell-specific EWAS
##
## Purpose of script: simulate EWAS with induced changes between cases and controls for both 
## cell-type specific & ubiquitous DMPs to compare analytical models
##
## Author: Eilis Hannon
##
## Date Created: 2022-07-01
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# DEFINE ANALYSIS FUNCTION
#----------------------------------------------------------------------#

runEWAS<-function(row,QCmetrics, status){

	modelLM<-lm(row ~ status * QCmetrics$Cell.type + QCmetrics$CCDNAmAge + QCmetrics$Sex + QCmetrics$Tissue.Centre)
	nullLM<-lm(row ~ status + QCmetrics$Cell.type + QCmetrics$CCDNAmAge + QCmetrics$Sex + QCmetrics$Tissue.Centre)

	data<-cbind(row,status, QCmetrics)
	modelMLM<-lmer(row ~ status * Cell.type + CCDNAmAge + Sex +  + (1 | Tissue.Centre)  + (1 | Indidivual.ID), REML = FALSE, data = QCmetrics)
	nullMLM<-lmer(row ~ status + Cell.type + CCDNAmAge + Sex +  + (1 | Tissue.Centre)  + (1 | Indidivual.ID), REML = FALSE, data = QCmetrics)

	p.df <- pdata.frame(data.frame("meth" = row, "phenotype" = status, "age" = QCmetrics$CCDNAmAge, "sex" = QCmetrics$Sex, "cell.type" = QCmetrics$Cell.type, "brain.bank" = QCmetrics$Tissue.Centre, "id" = QCmetrics$Indidivual.ID), index = c("id"), drop.index = F)
		
	modelp <- plm(meth ~ phenotype * cell.type + age + sex  + brain.bank, data = p.df, model = "pooling")
	nullCRR <- plm(meth ~ phenotype + cell.type + age + sex  + brain.bank, data = p.df, model = "pooling")
	# compute Stata like df-adjustment
	G <- length(unique(p.df$id))
	N <- length(p.df$id)
	dfa <- (G/(G - 1)) * (N - 1)/modelp$df.residual

	# display with cluster VCE and df-adjustment
	firm_c_vcov <- dfa * vcovHC(modelp, type = "HC0", cluster = "group", adjust = T)
	modelCRR<-coeftest(modelp, vcov = firm_c_vcov)


	return(c(summary(modelLM)$coefficients["status1", 4], anova(modelLM, nullLM)[2,6], 
	summary(modelMLM)$coefficients["status1", 5], anova(modelMLM, nullMLM)[2,8],
	 modelCRR["phenotype1",4],waldtest(modelp, nullCRR, vcov = firm_c_vcov, test = "F")[2,4]))
}


#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(lme4)
library(lmerTest)
library(plm)
library(lmtest)
#library(GenABEL)
library(doParallel)
library(devtools)
devtools::load_all(path = "../functionsR")

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

nSim<-10
nSig.options<-c(10,100,1000)
propCS.options<-seq(0,1,0.2)
sigEffect<-0.05

args<-commandArgs(trailingOnly = TRUE)
dataDir <- args[1]
nChunk<-args[2]

set.seed(nChunk)

normData<-file.path(dataDir, "/3_normalised/normalised.rdata")
resPath<-file.path(dataDir, "/4_analysis/methodsDevelopment")

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(dataDir)
load(normData)

## remove total samples and cell types with less than 20 samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell.type != "Total"),]
nSample<-table(QCmetrics$Cell.type)
QCmetrics<-QCmetrics[QCmetrics$Cell.type %in% names(nSample[which(nSample > 19)]),]

celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]

cellTypes<-unique(QCmetrics$Cell.type)

#----------------------------------------------------------------------#
# INTITATE PARALLEL ENV
#----------------------------------------------------------------------#

nCores<-detectCores()
cl <- makeCluster(nCores-1)
registerDoParallel(cl)
clusterExport(cl, list("runEWAS", "lmer", "pdata.frame", "plm", "vcovHC", "coeftest", "waldtest"))

#----------------------------------------------------------------------#
# RUN SIMULATIONS
#----------------------------------------------------------------------#

sumSim<-matrix(data = NA, nrow = nSim*length(nSig.options)*length(propCS.options), ncol = 2+(6*6))
colnames(sumSim)<-c("nProbes", "nCTspecific", paste("LM_ME", c("TotSig", "nTruePos", "nFalsePos","nTruePosCS", "nFalsePosCS", "lambda"), sep = "_"),
paste("LM_Int", c("TotSig", "nTruePos", "nFalsePos", "nTruePosCS", "nFalsePosCS", "lambda"), sep = "_"),
paste("MLM_ME", c("TotSig", "nTruePos", "nFalsePos", "nTruePosCS", "nFalsePosCS", "lambda"), sep = "_"),
paste("MLM_Int", c("TotSig", "nTruePos", "nFalsePos","nTruePosCS", "nFalsePosCS",  "lambda"), sep = "_"),
paste("CRR_ME", c("TotSig", "nTruePos", "nFalsePos", "nTruePosCS", "nFalsePosCS", "lambda"), sep = "_"),
paste("CRR_Int", c("TotSig", "nTruePos", "nFalsePos", "nTruePosCS", "nFalsePosCS", "lambda"), sep = "_"))

rowNum<-1
nullSim<-NULL
for(simNum in 1:nSim){
	message("Simulation: ", simNum, " of ", nSim)

	# randomly select samples to be cases, fix so numbers per brain bank match actual
	status<-rep(0,nrow(QCmetrics))
	for(centre in unique(QCmetrics$Tissue.Centre)){
		ids<-unique(QCmetrics$Indidivual.ID[which(QCmetrics$Tissue.Centre == centre)])
		cases<-sample(ids, floor(length(ids)/2))
		status[QCmetrics$Indidivual.ID %in% cases]<-1
	}
	status<-as.factor(status)
	
	# Run the null EWAS for this simulation
	outtab.null<-foreach(i=1:nrow(celltypeNormbeta), .combine = "rbind") %dopar% runEWAS(celltypeNormbeta[i,], QCmetrics, status)
	
	rownames(outtab.null)<-rownames(celltypeNormbeta)
	colnames(outtab.null)<-c(paste("LM", c("ME", "Int"), sep = "_"),paste("MLM", c("ME", "Int"), sep = "_"),paste("CRR", c("ME", "Int"), sep = "_"))
	nullSim<-cbind(nullSim, outtab.null)
	
	
	# Retest sites effects are induced at as non significiant probes are unaltered		
	for(nSig in nSig.options){
		for(propCS in propCS.options){
			sumSim[rowNum,1]<-nSig
			sumSim[rowNum,2]<-floor(nSig*propCS)
			# randomly select significant probes
			sigProbes<-sample(1:nrow(celltypeNormbeta), nSig)

			# create matrix with effects to add to sample profiles
			diffs<-rnorm(nrow(QCmetrics)*nSig, sigEffect, 0.005) 
			diffs<-matrix(data = diffs, nrow = nSig, byrow = TRUE)
			diffs[,which(status == 0)]<-0 # only add to cases
			
			# make some of these effects cell type specific
			if(propCS > 0){
				ctSpecific<-sample(1:nSig, floor(nSig*propCS))
				for(each in ctSpecific){
					# randomly select cell type to be significant in
					selectCT<-sample(cellTypes, 1)
					# set other cell types to have no effect
					diffs[each, !QCmetrics$Cell.type %in% selectCT]<-0
				}
				ctSpecific<-sigProbes[ctSpecific]
			} else {
				ctSpecific<-NULL
			}				

			# create matrix of betas with induced effects
			testbetas<-celltypeNormbeta[sigProbes,]+diffs

			outtab.sig<-foreach(i=1:nrow(testbetas), .combine = "rbind") %dopar% runEWAS(testbetas[i,], QCmetrics, status)
			
			# merge signif results with null results to generate summary statistics
			outtab.sim<-outtab.null
			outtab.sim[sigProbes,]<-outtab.sig
			
			sumSim[rowNum,2+seq(1,6*6, 6)]<-colSums(outtab.sim < 9e-8)
			sumSim[rowNum,3+seq(1,6*6, 6)]<-colSums(outtab.sim[sigProbes,] < 9e-8)
			sumSim[rowNum,4+seq(1,6*6, 6)]<-colSums(outtab.sim[-sigProbes,] < 9e-8)
			# handle quirk of R converting 1 row matrix to vector!!
			if(length(ctSpecific) > 2){
				sumSim[rowNum,5+seq(1,6*6, 6)]<-colSums(outtab.sim[ctSpecific,] < 9e-8)
				sumSim[rowNum,6+seq(1,6*6, 6)]<-colSums(outtab.sim[-ctSpecific,] < 9e-8)
			} else {
				if(length(ctSpecific) == 1){
					sumSim[rowNum,5+seq(1,6*6, 6)]<-as.numeric(outtab.sim[ctSpecific,] < 9e-8)
					sumSim[rowNum,6+seq(1,6*6, 6)]<-colSums(outtab.sim[-ctSpecific,] < 9e-8)
				} else{
					sumSim[rowNum,5+seq(1,6*6, 6)]<-NA
					sumSim[rowNum,6+seq(1,6*6, 6)]<-NA
				}
			}
						
			sumSim[rowNum,7+seq(1,6*6, 6)]<-apply(outtab.sim, 2, estlambda)
			rowNum<-rowNum+1
		}
	}
}

save(nullSim, file = file.path(resPath, paste0("nullSimulations_Chunk", nChunk, ".rdata")))
save(sumSim, file = file.path(resPath, paste0("nSigSimulateTrueEffects_Chunk", nChunk, ".rdata")))
