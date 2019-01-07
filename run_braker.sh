#!/bin/bash

#PBS -N braker
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M reubennowell@gmail.com

#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=8:mem=2gb

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
export PATH=/home/rnowell/software/augustus-3.2.1/bin/:/home/rnowell/software/augustus-3.2.1/scripts/:${PATH}
export PERL5LIB=/home/rnowell/perl5lib/:${PERL5LIB}

## load necessities
export AUGUSTUS_CONFIG_PATH=/home/rnowell/software/augustus-3.2.1/config/
export AUGUSTUS_BIN_PATH=/home/rnowell/software/augustus-3.2.1/bin/
export AUGUSTUS_SCRIPTS_PATH=/home/rnowell/software/augustus-3.2.1/scripts/
export BAMTOOLS_PATH=/home/rnowell/libraries/bamtools/bin/
export GENEMARK_PATH=/home/rnowell/software/gm_et_linux_64/gmes_petap/
export SAMTOOLS_PATH=/home/rnowell/software/samtools-1.2/
export ALIGNMENT_TOOL_PATH=/home/rnowell/software/exonerate/bin/exonerate
export BLAST_PATH=/home/rnowell/software/ncbi-blast-2.3.0+/bin/
export LD_LIBRARY_PATH=/home/rnowell/libraries/bamtools/lib/:${LD_LIBRARY_PATH}

## variables
THREADS=8
FASTA=
BAM=Aligned.sortedByCoord.out.bam
## species name
SPECIES=Sp_${RANDOM}

echo ------
echo BRAKER
echo ------
echo

## run braker.pl
~/software/BRAKER_v2.1.0/braker.pl \
--genome=${FASTA} \
--bam=${BAM} \
--species=${SPECIES} \
--gff3 \
--cores=${THREADS} \
--filterOutShort \
--useexisting

cd braker/${SPECIES}
~/software/BRAKER_v2.1.0/getAnnoFasta.pl --seqfile=../../${FASTA} augustus.gff

## env
source ~/virt_env/python/busco/bin/activate

## ln to lineage
ln -s ~/software/busco/eukaryota_odb9
ln -s ~/software/busco/metazoa_odb9

echo -----
echo BUSCO
echo -----
echo

## BUSCO on proteins
~/software/busco/scripts/run_BUSCO.py \
-i augustus.aa \
-o eukaryota_odb9.augustus.aa \
-l eukaryota_odb9/ \
-m prot \
-c $THREADS \
--limit 4

## BUSCO on proteins
~/software/busco/scripts/run_BUSCO.py \
-i augustus.aa \
-o metazoa_odb9.augustus.aa \
-l metazoa_odb9/ \
-m prot \
-c $THREADS \
--limit 4

echo --------
echo FINISHED
echo --------
echo
date

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
