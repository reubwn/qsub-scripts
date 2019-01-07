#!/bin/sh

#PBS -N sspace
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M foo.bar@email.com

#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=16:mem=100gb

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

## files
LIBFILE=
FASTA=
SUFFIX=

## sspace
perl ~/software/SSPACE-STANDARD-3.0_linux-x86_64/SSPACE_Standard_v3.0.pl \
-l $LIBFILE \
-s $FASTA \
-k 5 \
-T 64 \
-b sspace_${SUFFIX}

## rm tmp files
rm -r sspace_${SUFFIX}/alignoutput/ sspace_${SUFFIX}/intermediate_results/ sspace_${SUFFIX}/pairinfo/ sspace_${SUFFIX}/reads/
touch cleanup.ok

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
