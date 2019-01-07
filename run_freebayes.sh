#!/bin/sh

#PBS -N freebayes
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M reubennowell@gmail.com

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=10gb

## NB values for ncpus and mem are allocated
## to each node (specified by select=N)
##
## to direct output to cwd, use $PBS_O_WORKDIR:
## specify LOGFILE found in ~/ during execution then moved to cwd on job completion
##
cd $PBS_O_WORKDIR
JOBNUM=`echo $PBS_JOBID | sed 's/\..*//'`
LOGFILE=${PBS_JOBNAME}.o${JOBNUM}

#########################################
##                                     ##
## Output some useful job information. ##
##                                     ##
#########################################

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: job number is $JOBNUM
echo PBS: logfile is $LOGFILE
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
#echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------

## variables
REF=
BAM=
OUTFILE=${BAM}.FreeBayes.standard-filters.vcf

## freebayes
## --standard-filters equivalent to -m 30 -q 20 -R 0 -S 0
~/software/freebayes/bin/freebayes \
-f $REF \
--ploidy 2 \
--standard-filters \
--min-alternate-fraction 0.25 \
--min-coverage 10 \
--max-coverage 300 \
--no-indels \
--no-mnps \
--no-complex \
$BAM > $OUTFILE

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
