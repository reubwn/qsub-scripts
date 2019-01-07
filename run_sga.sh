#!/bin/bash

#PBS -N sga
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M reubennowell@gmail.com

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=16:mem=200gb

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

## source sga environment
source /home/rnowell/virt_env/python/sga/bin/activate
export BAMTOOLS_PATH=/home/rnowell/libraries/bamtools/bin/
export LD_LIBRARY_PATH=/home/rnowell/libraries/bamtools/lib/:${LD_LIBRARY_PATH}

READS1=
READS2=
PREFIX=
THREADS=

## sga pipeline
/home/rnowell/software/sga/bin/sga preprocess --pe-mode=1 -o ${PREFIX}.shuff.fa ${READS1} ${READS2}
/home/rnowell/software/sga/bin/sga index -a ropebwt -t $THREADS ${PREFIX}.shuff.fa
/home/rnowell/software/sga/bin/sga preqc -t $THREADS ${PREFIX}.shuff.fa > ${PREFIX}.shuff.preqc
/home/rnowell/software/sga/src/bin/sga-preqc-report.py ${PREFIX}.shuff.preqc

## tidy up
rm ${PREFIX}.shuff.fa

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
