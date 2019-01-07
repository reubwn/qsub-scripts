#!/bin/bash

#PBS -N diamond
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M foo.bar@email.com

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

####################### diamond

## define query
QUERY=
PREFIX=
DB=~/ncbi/db/uniref90.fasta

## make dmnd db if needs be
#~/bin/diamond makedb \
#--in ${QUERY} \
#-d ${QUERY}

## run diamond
~/bin/diamond blastx \
--max-target-seqs 100 \
--sensitive \
--index-chunks 1 \
-e 1e-25 \
-p 64 \
-q $QUERY \
-d $DB \
-a ${PREFIX}.vs.uniref90.diamond

## convert to .tagc for blobtools
~/scripts/daa_to_tagc.pl ~/ncbi/db/uniref90.taxlist ${PREFIX}.vs.uniref90.diamond.daa

## get tab outfile
~/bin/diamond view \
-a ${PREFIX}.vs.uniref90.diamond.daa \
-o ${PREFIX}.vs.uniref90.diamond.daa.out

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
