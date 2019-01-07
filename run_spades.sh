#!/bin/sh

#PBS -N spades
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M reubennowell@gmail.com

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=64:mem=300gb

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
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------

## define files
READS1=
READS2=
SUFFUX=

## load GOMP libraries needed by SPAdes
##
export LD_LIBRARY_PATH=/home/rnowell/libraries/gcc/gcc-4.9.2-build/gcc-4.9.2/lib64/

## command
##
/home/rnowell/software/SPAdes-3.9.0-Linux/bin/spades.py \
-o spades_${SUFFUX} \
-1 $READS1 \
-2 $READS2 \
-k 21,33,55,77 \
--careful \
--only-assembler \
--threads 64 \
--memory 250

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
