#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q mrcq # submit to the serial queue
#PBS -l walltime=10:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC190311 # research project to submit under. 
#PBS -l procs=1 # specify number of processors.
#PBS -e SNPArray/logFiles/QCSNPdata.err # error file
#PBS -o SNPArray/logFiles/QCSNPdata.log # output file


## print start date and time
echo Job started on:
date -u

####### 

## NOTE: Do not store confidential information in this file use the config file

######


source ./$1



sh SNPArray/preprocessing/1_QC.sh

module load R/3.6.3-foss-2020a
sh SNPArray/preprocessing/2_CheckEthnicity.sh


## create sample sheets of samples classed as different populations
## run relatedness check within ethnicites
populations=($(cut -f3 --delim="," ${PROCESSDIR}/PredictedPopulations.csv | tail -n +2 | sort | uniq))
for each in ${populations[@]}
do
   grep ${each} ${PROCESSDIR}/PredictedPopulations.csv | cut -f1-2 --delim="," --output-delimiter=" " > ${PROCESSDIR}/${each}Samples.txt
   
   ## 
   if [[ $(wc -l <${PROCESSDIR}/${each}Samples.txt) -ge 1 ]]
   then
      sh SNPArray/preprocessing/3_CheckRelatedness.sh ${PROCESSDIR}/${each}Samples.txt ${each}
   fi
   
done   

## plot relatedness statistics acroos whole cohort 

Rscript ${SCRIPTDIR}/utilitys/plotKinshipCoeff.r ${PROCESSDIR}/QCoutput/

module purge
module load VCFtools
sh SNPArray/preprocessing/4_formatForImputation.sh ALL ${KGG}/1000GP_Phase3_combined.legend
sh SNPArray/preprocessing/4_formatForImputation.sh EUR ${KGG}/../HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab

## print finish date and time
echo Job finished on:
date -u
