#!/bin/sh

#PBS -N skewer
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M foo.bar@email.com

#PBS -l walltime=10:00:00
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

READS1=
READS2=
OUTPREFIX=
## for pe reads
ADAPTERS=/home/rnowell/adapters.fa
## for mp reads
ADAPTERS=/home/rnowell/data/adapters/adapters.fa

## count bases before
zcat $READS1 $READS2 | paste - - - - | cut -f2 | tr -d '\n' | wc -c > num_bases.raw.txt

## skewer
/home/rnowell/software/skewer-0.2.2/skewer \
-x ${ADAPTERS} \
-m pe \
-q 28 \
-Q 30 \
-l 75 \
-n \
-o ${OUTPREFIX} \
-z \
-t 64 \
${READS1} ${READS2}

## fastqc
source ~/.jdk/7
/home/rnowell/software/fastqc-v0.11.5/fastqc --threads 4 --nogroup $READS1 $READS2 ${OUTPREFIX}-trimmed-pair*

## count bases after
zcat ${OUTPREFIX}-trimmed-pair* | paste - - - - | cut -f2 | tr -d '\n' | wc -c > num_bases.trimmed.txt

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
