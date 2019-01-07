#!/bin/bash

#PBS -N vc.AK11
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M foo.bar@email.com

#PBS -l walltime=10:00:00
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

#Java
source ~/.jdk/8
#BBMap to PATH
BBMAP=~/software/bbmap-37.82/
export PATH=${PATH}:${BBMAP}

#Set files
READS1=../AK11_1.trimmed.ecc.filtered.fq.gz
READS2=../AK11_2.trimmed.ecc.filtered.fq.gz
REF=../scaffolds_final.AK11.fasta
THREADS=8
ln -s $REF reference.fa

## interleave reads
echo ''
echo 'Merging reads'
echo '-------------'
echo ''
seqtk mergepe $READS1 $READS2 > interleaved.fq

#Remove duplicates
#Optical deduplication requires standard Illumina read headers and will not work with renamed reads, such as most SRA data.
#To perform PCR-duplicate removal (of all duplicates regardless of location), omit the "optical" flag.
#Deduplication is generally not recommended for quantification experiments such as RNA-seq.
echo ''
echo 'Remove duplicates'
echo '-----------------'
echo ''
clumpify.sh -Xmx100g in=interleaved.fq out=clumped.fq.gz dedupe optical

#Remove low-quality regions
#This step requires standard Illumina read headers and will not work with renamed reads, such as most SRA data.
echo ''
echo 'Remove low-qual'
echo '---------------'
echo ''
filterbytile.sh -Xmx100g in=clumped.fq.gz out=filtered_by_tile.fq.gz

#Trim adapters
echo ''
echo 'Trim adapters'
echo '-------------'
echo ''
bbduk.sh -Xmx100g threads=${THREADS} in=filtered_by_tile.fq.gz out=trimmed.fq.gz ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=100 ref=${BBMAP}/resources/adapters.fa ftm=5 ordered

#Remove synthetic artifacts and spike-ins.  Add "qtrim=r trimq=8" to also perform quality-trimming at this point, but not if quality recalibration will be done later.
bbduk.sh -Xmx100g threads=${THREADS} in=trimmed.fq.gz out=filtered.fq.gz k=27 ref=${BBMAP}/resources/sequencing_artifacts.fa.gz,${BBMAP}/resources/phix174_ill.ref.fa.gz ordered

#Map to reference
echo ''
echo 'Map reads'
echo '---------'
echo ''
bbmap.sh -Xmx100g threads=${THREADS} in=filtered.fq.gz out=mapped.sam.gz pigz unpigz ref=reference.fa nodisk

#Call variants
echo ''
echo 'Call variants'
echo '-------------'
echo ''
callvariants.sh -Xmx100g in=mapped.sam.gz out=variants.txt vcf=variants.vcf.gz ref=reference.fa ploidy=2 prefilter

#Do some things
perl -i.bak -lne '$_=~s/\t\t/\t-\t/;print' variants.txt

#Delete files
echo ''
echo 'Delete files'
echo '------------'
echo ''
rm interleaved.fq mapped.sam.gz *fq.gz

echo `date`

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
