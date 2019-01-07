#!/bin/sh

#PBS -N repeatmodeler
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M reubennowell@gmail.com

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=100:mem=100gb

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

## build database
/home/rnowell/software/RepeatModeler/BuildDatabase -name rmagnacalcarata -engine ncbi /home/rnowell/bdelloid/assembly/freeze/Rmag-1.3.scaffolds.fna

## run rm
/home/rnowell/software/RepeatModeler/RepeatModeler -database /home/rnowell/bdelloid/annotation/Rmag/repeatmodeler/rmagnacalcarata -engine ncbi -pa 100

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
