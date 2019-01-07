#!/bin/bash

#PBS -N socrates
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M reubennowell@gmail.com

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:mem=105gb

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

## load java 1.8
##
source ~/.jdk/8

## files
##
REF=
BAM=
THREADS=

## command timed to get mem and wallclock info
##
/usr/bin/time -v \
java -Xmx100g -jar ~/software/socrates/socrates-1.13.1-jar-with-dependencies.jar \
${REF}.bwt \
$BAM \
--threads $THREADS \
--verbose

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
