#!/bin/sh

#PBS -N RSEM.t
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M reubennowell@gmail.com

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=64:mem=250gb

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

## PATH settings
export PATH=/home/rnowell/software/bowtie/:/home/rnowell/software/kallisto_linux-v0.43.0/:/home/rnowell/software/RSEM-1.2.31/bin/:/home/rnowell/software/eXpress-1.5.1/linux_build/src/:${PATH}

## files
FASTA=Trinity.fasta
LEFT=
RIGHT=
#MAP=Trinity.fasta.gene_trans_map
OUTDIR=rsem_outdir
THREADS=64

## RSEM
/home/rnowell/software/trinityrnaseq-2.2.0/util/align_and_estimate_abundance.pl \
--transcripts $FASTA \
--seqType fq \
--left $LEFT \
--right $RIGHT \
--est_method RSEM \
--aln_method bowtie \
--prep_reference \
--trinity_mode \
--output_dir $OUTDIR \
--thread_count $THREADS

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
