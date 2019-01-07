#!/bin/sh

#PBS -N velv
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M foo.bar@email.com

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=64:mem=400gb

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

## load common modules as standard
##
module load samtools/1.2 blast+/2.2.28 intel-suite

## set num threads
export OMP_NUM_THREADS=64

## define env vars
OUTDIR=
KMER=
READS1=
READS2=
INS_LENGTH=

qalter -N ${PBS_JOBNAME}_${OUTDIR} $PBS_JOBID

## velveth
/home/rnowell/software/velvet/velveth $OUTDIR $KMER \
-fastq \
-separate \
-shortPaired \
$READS1 \
$READS2

## velveth multi_k
#/home/rnowell/software/velvet/velveth out_dir m,M,s \
#-fastq \
#-separate \
#-shortPaired \
#${READS1} \
#${READS2}

## velvetg
/home/rnowell/software/velvet/velvetg \
$OUTDIR \
-exp_cov auto \
-cov_cutoff auto \
-ins_length $INS_LENGTH

## delete big files to save space
cd $PBS_O_WORKDIR/$OUTDIR
rm Graph2 PreGraph LastGraph Roadmaps Sequences

## move LOGFILE to cwd
LOGFILE_NEW=${PBS_JOBNAME}_${OUTDIR}.o${JOBNUM}
mv $HOME/$LOGFILE $PBS_O_WORKDIR/$LOGFILE_NEW
