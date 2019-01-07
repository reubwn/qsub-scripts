#!/bin/bash

#PBS -N blast_nt
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M foo.bar@email.com

#PBS -l walltime=48:00:00
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

## define vars
QUERY=
PREFIX=
DB=~/ncbi/db/nt

printf "\n[INFO] BLAST: ${QUERY} vs ${DB}"

## run megablast vs nt
/home/rnowell/software/ncbi-blast-2.3.0+/bin/blastn \
-task megablast \
-num_threads 100 \
-evalue 1e-25 \
-culling_limit 5 \
-query $QUERY \
-db $DB \
-outfmt '6 qseqid staxids bitscore std' \
-out ${PREFIX}.vs.nt.megablast.1e25.out

printf "\n[INFO] Finished\n"

## move LOGFILE to cwd
mv $HOME/$LOGFILE $PBS_O_WORKDIR
