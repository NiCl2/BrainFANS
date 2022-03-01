#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=18:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACSeq/logFiles/%u/ATACQCSummary-%A.o
#SBATCH --error=ATACSeq/logFiles/%u/ATACQCSummary-%A.e
#SBATCH --job-name=ATACQCSummary

## print start date and time
echo Job started on:
date -u
	
## needs to be executed from the scripts folder
echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

## load config file provided on command line when submitting job
echo "Loading config file for project: " $1
export PROJECT=$1

source ./ATACSeq/config/config.txt 
echo "Project directory is: " $DATADIR

module load MultiQC
## use multiqc to collate QC output statistics

mkdir -p ${FASTQCDIR}/multiqc
cd ${FASTQCDIR}/
multiqc . -f -o ${FASTQCDIR}/multiqc

## remove redundant html files
rm -f *.html
rm -f ${TRIMDIR}/fastp_reports/*.html

mkdir -p ${ALIGNEDDIR}/multiqc
cd ${ALIGNEDDIR}/
multiqc . -f -o ${ALIGNEDDIR}/multiqc

## run other bespoke utilty scripts to collate other QC metrics
cd ${SCRIPTDIR}/

./ATACSeq/preprocessing/8_progressReport.sh 
./ATACSeq/preprocessing/9_countMTReads.sh 
./ATACSeq/preprocessing/10_collateFlagStatOutput.sh 

echo 'EXITCODE: ' $?

## move log files into a folder
cd ${SCRIPTDIR}/ATACSeq/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv ATACQCSummary-${SLURM_ARRAY_JOB_ID}* ${SLURM_ARRAY_JOB_ID}/
