#!/bin/sh

#PBS -N star
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M reubennowell@gmail.com

#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=8:mem=100gb

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
THREADS=8
GENOMEDIR=genomeDir_STAR
FASTA=
READS1=
READS2=

mkdir $GENOMEDIR

## generate index
/home/rnowell/software/STAR/STAR \
--runThreadN $THREADS \
--runMode genomeGenerate \
--genomeDir $GENOMEDIR \
--genomeFastaFiles $FASTA

## align reads
/home/rnowell/software/STAR/STAR \
--runThreadN $THREADS \
--runMode alignReads \
--genomeDir $GENOMEDIR \
--outSAMtype BAM SortedByCoordinate \
--readFilesCommand zcat \
--readFilesIn $READS1 $READS2 \
--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--alignIntronMin 20 \
--alignIntronMax 100000 \ ## set to max intron
--twopassMode Basic  ## improves SJ annotation apparently

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
