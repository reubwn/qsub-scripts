#!/bin/sh

#PBS -N bbnorm
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M foo.bar@email.com

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=64:mem=500gb

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

## load zlib and java
source /home/rnowell/.jdk/8
export LD_LIBRARY_PATH=/home/rnowell/libraries/zlib/zlib-1.2.8/lib/:${LD_LIBRARY_PATH}

READS1=
READS2=
PREFIX=
COV=    # target coverage

## run bbnorm
/home/rnowell/software/bbmap-36.02/bbnorm.sh \
-Xmx400g \
threads=64 \
in=${READS1} \
in2=${READS2} \
out=${PREFIX}_1.diginorm${COV}.fq.gz \
out2=${PREFIX}_2.diginorm${COV}.fq.gz \
hist=${PREFIX}.in.khist \
histout=${PREFIX}.out.khist \
peaks=${PREFIX}.peaks \
prefilter=t \
ecc=t \
fixspikes=t \
target=${COV}

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
