
## General 
# With regard to forking a repository to only modify and use two files - 
# As far as I understand it, this is the default way of doing it and at least is version-controlled. 
# There are ways to do a 'sparse' checkout for a few files (i.e. https://stackoverflow.com/questions/2466735/how-to-sparsely-checkout-only-one-single-file-from-a-git-repository) 
# but I would personally either stick with the way you've done it, or perhaps create a new script within the original repo that was specific to mouse 
# i.e. called by a different jobSubmission/_.sh 

# as an aside:
# i think mac always adds a .DS_Store to things? could add to gitignore?

#----------------------------------------#

## fitModelsAll.r
## is refPanelPath is being used?
## if you're not using IDOL do you need to do some of the data manipulation that is going on (i.e. line 69 is a rename for use with IDOL?)

# 31 
## to be able to run on command line or interactively without hashing out parts etc:
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
    setwd("~/OneDrive - University of Exeter/Documents/methylation/RawData_Mouse_DNAm/")
    args[1]<-"MouseArray_CellDeconv_FilteredNormalised_Betas.rdat"
    args[2]<-"../MouseMethylation-12v1-0_A2.csv"
    args[3]<-""
} 
normData<-args[1]
refPath<-args[2]
refPanelPath<-args[3]


#63
probeAnnot <- read.csv(refPath, skip = 7) %>%

# 69
# I assume that this line is simply more specific than the line Eilis has used? - so her script could use this and be more generalisable? 
# (colnames(pheno.all)[3]<-"CellType" # need to rename for use with IDOL functions)

#----------------------------------------#

## profileCellCompBulkBrainSamples.r
## I think the general comment is even more relevant for this one! Seems most of it has been deleted and with some sections as completely written from new.

# 28 
## to be able to run on command line or interactively without hashing out parts etc:
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


# 56
# nitpicking
## could leave as normbetas and change the ##66 and ##67 instances of betas

