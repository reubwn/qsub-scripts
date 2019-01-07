#!/bin/bash

#PBS -N tadpole
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M reubennowell@gmail.com

#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=16:mem=500gb

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


## define variables
IN1=
IN2=
OUT1=
OUT2=
OUTPREFIX=

## load bbmap environment
source /home/rnowell/.jdk/8
export LD_LIBRARY_PATH=/home/rnowell/libraries/zlib/zlib-1.2.8/lib/:${LD_LIBRARY_PATH}

## get kmer dist of input
~/software/bbmap-37.82/kmercountexact.sh \
-Xmx450g \
threads=64 \
in1=$IN1 \
in2=$IN2 \
khist=${OUTPREFIX}.khist \
peaks=${OUTPREFIX}.peaks

## run tadpole
~/software/bbmap-37.82/tadpole.sh \
-Xmx450g \
threads=64 \
mode=correct \
ecc=t \
in1=$IN1 \
in2=$IN2 \
out1=$OUT1 \
out2=$OUT2 \
tossdepth=1

## get kmer dist of result
~/software/bbmap-37.82/kmercountexact.sh \
-Xmx450g \
threads=64 \
in1=$OUT1 \
in2=$OUT2 \
khist= \
peaks=

## count bases after
zcat $OUT1 $OUT2 | paste - - - - | cut -f2 | tr -d '\n' | wc -c > ${OUTPREFIX}.num_bases.ecc.txt

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
