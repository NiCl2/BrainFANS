## run this script to calculate QC metrics
## first it excludes samples with really low intensities values
## requires gdsfile to already be generated and sample sheet to be loaded
##

setwd(dataDir)

gfile<-openfn.gds(gdsFile, readonly = FALSE, allow.fork = TRUE)

## see if any QC data already exists
if(file.exists(qcData)){
	load(qcData)
	print("QC file loaded")
} else{
	QCmetrics<-sampleSheet
}

QCmetrics$Age<-as.numeric(as.character(QCmetrics$Age))

## extract a few useful matrices
rawbetas<-betas(gfile)[,]

# calculate median M & U intensity
if(!"M.median" %in% colnames(QCmetrics)){
	
	m_intensities<-methylated(gfile)
	u_intensities<-unmethylated(gfile)
	M.median<-unlist(apply.gdsn(m_intensities, 2, median))
	U.median<-unlist(apply.gdsn(u_intensities, 2, median))
	intens.ratio<-M.median/U.median
	## exclude really poor intensity samples from beginning so rest of QC is not dominated by them
	intensPASS<-M.median > 500
	QCmetrics<-cbind(QCmetrics,M.median, U.median, intens.ratio, intensPASS)

}


 

# calculate bisulfite conversion statistic
if(!"bisulfCon" %in% colnames(QCmetrics)){	
	bisulfCon<-bscon(gfile)
	QCmetrics<-cbind(QCmetrics,bisulfCon)
}

## PCA of control-probe intensities
if(!"PC1_cp" %in% colnames(QCmetrics)){	
	## exclude really poor intensity samples
	qc.meth<-QCmethylated(gfile)[,QCmetrics$intensPASS]
	qc.unmeth<-QCunmethylated(gfile)[,QCmetrics$intensPASS]
	# remove negative controls
	qc.meth<-qc.meth[grep("Negative", rownames(qc.meth), invert=TRUE),]
	qc.unmeth<-qc.unmeth[grep("Negative", rownames(qc.unmeth), invert=TRUE),]
	ctrl.all<-t(rbind(qc.meth, qc.unmeth))

	pca <- prcomp(na.omit(ctrl.all))
	ctrlprobes.scores = pca$x
	colnames(ctrlprobes.scores) = paste(colnames(ctrlprobes.scores), '_cp', sep='')
	ctrl.pca<-pca$sdev^2/sum(pca$sdev^2)
	ctrlprobes.scores<-ctrlprobes.scores[match(QCmetrics$Basename, QCmetrics$Basename[QCmetrics$intensPASS]),]
	rownames(ctrlprobes.scores)<-QCmetrics$Basename	
	## only save PCs which explain > 1% of the variance
	QCmetrics<-cbind(QCmetrics,ctrlprobes.scores[,which(ctrl.pca > 0.01)])
}

## perform PCA on beta values
if(!"PC1_betas" %in% colnames(QCmetrics)){
	## filter to autosomal only
	auto.probes<-which(fData(gfile)$chr != "chrX" & fData(gfile)$chr != "chrY")

	pca <- prcomp(t(na.omit(rawbetas[auto.probes,QCmetrics$intensPASS])))
	betas.scores = pca$x
	colnames(betas.scores) = paste(colnames(betas.scores), '_betas', sep='')
	betas.pca<-pca$sdev^2/sum(pca$sdev^2)
	betas.scores<-betas.scores[match(QCmetrics$Basename, QCmetrics$Basename[QCmetrics$intensPASS]),]
	rownames(betas.scores)<-QCmetrics$Basename	
	## only save PCs which explain > 1% of the variance
	QCmetrics<-cbind(QCmetrics,betas.scores[,which(betas.pca > 0.01)])

}

## Identify outlier samples 
## excluded as comparable to PC filtering already included
#if(!"iqr" %in% colnames(QCmetrics)){
	#outlierDetect<- outlyx(rawbetas, plot = FALSE)
	#QCmetrics<-cbind(QCmetrics,outlierDetect)
#}

## detection p value filtering at this stage only interested in sample filtering, will repeat later for probe filtering
if(!"pFilter" %in% colnames(QCmetrics)){	

	pFOut<-pfilter.gds(gfile, pn = pvals(gfile), bc = index.gdsn(gfile, "NBeads"))
	QCmetrics<-cbind(QCmetrics,"pFilter"= pFOut$samples)
}


## calc Horvaths epigenetic age
if(!"DNAmAge" %in% colnames(QCmetrics)){	
	data(coef)
	DNAmAge<-agep(gfile, coef=coef)
	QCmetrics<-cbind(QCmetrics,DNAmAge)
}

## check sex
# idenitfy X chromosome probes
if(!"predSex" %in% colnames(QCmetrics)){	
	x.probes<-which(fData(gfile)$chr == "chrX")
	y.probes<-which(fData(gfile)$chr == "chrY")
	ints.auto<-methylated(gfile)[c(x.probes, y.probes),]+unmethylated(gfile)[c(x.probes, y.probes),]
	ints.X<-methylated(gfile)[x.probes,]+unmethylated(gfile)[x.probes,]
	ints.Y<-methylated(gfile)[y.probes,]+unmethylated(gfile)[y.probes,]

	x.cp<-colMeans(ints.X, na.rm = TRUE)/colMeans(ints.auto, na.rm = TRUE)
	y.cp<-colMeans(ints.Y, na.rm = TRUE)/colMeans(ints.auto, na.rm = TRUE)
	
	## base prediction on y chromosome
	predSex.y<-rep(NA, length(y.cp))
	predSex.y[which(y.cp > 1)]<-"M"
	predSex.y[which(y.cp < 1)]<-"F"
	
	## base prediction on x chromosome
	predSex.x<-rep(NA, length(x.cp))
	predSex.x[which(x.cp < 1)]<-"M"
	predSex.x[which(x.cp > 1)]<-"F"
	
	## check for consistent prediction
	predSex<-rep(NA, length(x.cp))
	predSex[which(predSex.x == predSex.y)]<-predSex.x[which(predSex.x == predSex.y)]
	QCmetrics<-cbind(QCmetrics,x.cp,y.cp,predSex.x, predSex.y, predSex)
}

## check duplicate samples using SNPs on array
if(!exists("snpCor")){
	rsbetas<-rawbetas[grep("rs", rownames(rawbetas)),]
	snpCor<-cor(rsbetas, use = "pairwise.complete.obs")
}

# check effect of normalisation
if(!"rmsd" %in% colnames(QCmetrics)){
	dasen(gfile, node="normbeta")
	normbetas<-index.gdsn(gfile, "normbeta")[,]
	qualDat<-qual(rawbetas, normbetas)
	QCmetrics<-cbind(QCmetrics,qualDat)
}

## need to close gds file in order to open in another R session
closefn.gds(gfile)

# save QC metrics and SNP correlations to generate QC report
save(QCmetrics, snpCor, betas.pca, ctrl.pca, pFOut, file = qcData)