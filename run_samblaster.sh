#!/bin/bash

#PBS -N samblaster
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M foo.bar@email.com

#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=1:mem=8gb

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

## htslib
export LD_LIBRARY_PATH=/home/rnowell/software/htslib/lib

## env vars
REF=../../../../AK11/scaffolds_final.AK11.fasta
READS1=../../../AK16_1.trimmed.ecc.filtered.fq.gz
READS2=../../../AK16_2.trimmed.ecc.filtered.fq.gz
PREFIX=AK16_reads.vs.AK11_ref

if [ ! -f ${REF}.fai ]; then
  samtools faidx $REF
fi

if [ ! -f ${REF}.sa ]; then
  ~/software/bwa-0.7.12/bwa index $REF
fi

## run bwa mem | samblaster
~/software/bwa-0.7.12/bwa mem \
$REF \
$READS1 \
$READS2 \
| samblaster -e -d ${PREFIX}.DISC.sam -s ${PREFIX}.SPLIT.sam | samtools sort -O bam > ${PREFIX}.OUT.sorted.bam

## convert samblaster sam to bam
parallel 'samtools sort -O bam -o {.}.sorted.bam {}' ::: ${PREFIX}.DISC.sam ${PREFIX}.SPLIT.sam

## index bams
parallel 'samtools index {}' ::: ${PREFIX}.DISC.sorted.bam ${PREFIX}.SPLIT.sorted.bam ${PREFIX}.OUT.sorted.bam

## gather some stats
parallel 'samtools stats {} > {}.samstats' ::: ${PREFIX}.DISC.sorted.bam ${PREFIX}.SPLIT.sorted.bam ${PREFIX}.OUT.sorted.bam
parallel 'bamtools stats -in {} -insert > {}.bamstats' ::: ${PREFIX}.DISC.sorted.bam ${PREFIX}.SPLIT.sorted.bam ${PREFIX}.OUT.sorted.bam 

## generate genomecov files
parallel 'bedtools genomecov -ibam {} -d | gzip > {}.cov.gz' ::: ${PREFIX}.DISC.sorted.bam ${PREFIX}.SPLIT.sorted.bam ${PREFIX}.OUT.sorted.bam

## move LOGFILE to cwd
mv $HOME/$LOGFILE $PBS_O_WORKDIR
