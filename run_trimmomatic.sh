#!/bin/sh

#PBS -N trimmomatic
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M reubennowell@gmail.com

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=64:mem=100gb

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

source /home/rnowell/.jdk/7
java -version

READS1=
READS2=
PREFIX=

## trimmo command
java -jar /home/rnowell/software/Trimmomatic-0.36/trimmomatic-0.36.jar \
PE \
-threads 64 \
-phred33 \
${READS1} \
${READS2} \
${PREFIX}_1.fq.gz \
${PREFIX}_1u.fq.gz \
${PREFIX}_2.fq.gz \
${PREFIX}_2u.fq.gz \
ILLUMINACLIP:/home/rnowell/adapters.fa:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:25 \
MINLEN:75

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
