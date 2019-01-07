#!/bin/bash

#PBS -N braker1
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M reubennowell@gmail.com

#PBS -l walltime=15:00:00
#PBS -l select=1:ncpus=8:mem=10gb

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

## load perlbrew environment
source /home/rnowell/perl5/perlbrew/etc/bashrc

## load necessities
export PATH=/home/rnowell/software/augustus-3.2.1/bin/:${PATH}
export PATH=/home/rnowell/software/augustus-3.2.1/scripts/:${PATH}
export PATH=/home/rnowell/software/blat/:${PATH}
export AUGUSTUS_CONFIG_PATH=/home/rnowell/software/augustus-3.2.1/config/
export BAMTOOLS_PATH=/home/rnowell/libraries/bamtools/bin/
export GENEMARK_PATH=/home/rnowell/software/genemark-es-4.32/
export SAMTOOLS_PATH=/home/rnowell/software/samtools-1.2/samtools
export LD_LIBRARY_PATH=/home/rnowell/libraries/bamtools/lib/:${LD_LIBRARY_PATH}

## files
FASTA=
BAM=Aligned.sortedByCoord.out.bam
THREADS=8
SPECIES=

## run braker.pl
/home/rnowell/software/braker-1.9/braker.pl \
--genome=${FASTA} \
--bam=${BAM} \
--cores=${THREADS} \
--species=${SPECIES}

## use --useexisting flag if optimization already done
#--useexisting

## get CDS and exon fastas
cd braker/${SPECIES}
~/software/augustus-3.2.1/scripts/getAnnoFasta.pl --seqfile=genome.fa augustus.gff
~/software/augustus-3.2.1/scripts/gtf2gff.pl <augustus.gff --out=augustus.gff3 --gff3

## tidy up
rm genome.split*
rm Align.bam
rm -r GeneMark-ET/

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
