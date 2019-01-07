#!/bin/bash

#PBS -N mcclintock
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M r.nowell@imperial.ac.uk

#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=32:mem=42gb

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

## load conda env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate MCCLINTOCK

sleep 1

echo Using samtools at: which samtools
echo Samtools version:
samtools

## load java
source ~/.jdk/8

## vars
REF=../../scaffolds_final.AK11.fasta
#OUTDIR=`basename $REF`
OUTDIR=AK15
PE1=../../../AK15/AK15_1.trimmed.ecc.filtered.fq.gz
PE2=../../../AK15/AK15_2.trimmed.ecc.filtered.fq.gz
TE=~/medbio/software/RepeatMasker/Libraries/Repbase_custom/metazoa_repeatlib.fa
THREADS=16

echo
echo `date`
echo
echo VARIABLES
echo Refence: $REF
echo TE fasta: $TE
echo PE1: $PE1
echo PE2: $PE2
echo Threads: $THREADS
echo

## link to files
mkdir -p $OUTDIR/fq_dir
ln -s $PWD/$REF $OUTDIR/
ln -s $TE $OUTDIR
ln -s $PWD/$PE1 $OUTDIR/
ln -s $PWD/$PE2 $OUTDIR/

echo "Inflating reads into $OUTDIR/fq_dir/ ..."
parallel ::: "unpigz -p 8 $PE1 > $OUTDIR/fq_dir/reads_1.fastq" "unpigz -p 8 $PE2 > $OUTDIR/fq_dir/reads_2.fastq"

## don't begin until dir structure has been created
while [ ! -d $PWD/$OUTDIR ]; do sleep 1; done

## McClintock command
##
~/medbio/software/mcclintock/mcclintock.sh \
-r $PWD/$REF \
-c $TE \
-1 $PWD/$OUTDIR/fq_dir/reads_1.fastq \
-2 $PWD/$OUTDIR/fq_dir/reads_2.fastq \
-o $PWD/$OUTDIR \
-b \
-i \
-p $THREADS \
-M 16

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR/$OUTDIR
