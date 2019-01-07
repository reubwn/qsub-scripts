#!/bin/bash

#PBS -N jitterbug
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M foo.bar@email.com

#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=8:mem=20gb

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
conda activate jitterbug

## load java
source ~/.jdk/8

## declare variables
REF=
PE1=
PE2=
BAM=sorted.bam
GFF=
TAG=ID
SAMPLE=
OUTDIR=jitterbug_${SAMPLE}
mkdir $OUTDIR
chmod -R 777 $OUTDIR
THREADS=8

start=`date +%s`

echo "Start: "`date`
echo
echo "VARIABLES"
echo "Refence: $REF"
echo "PE1: $PE1"
echo "PE2: $PE2"
echo "Threads: $THREADS"
echo

## index reference
parellel ::: "samtools faidx $REF" "/home/rnowell/software/bwa-0.7.12/bwa index $REF"

## run bwa mem
/home/rnowell/software/bwa-0.7.12/bwa mem \
-Y \
-t $THREADS \
$REF \
$PE1 \
$PE2 \
| samtools view -@ $THREADS -b - \
| samtools sort -@ $THREADS -O bam -T temp_${RANDOM} - > $OUTDIR/$BAM

## index bam
samtools index $OUTDIR/sorted.bam

## gather some stats
bamtools stats -in $OUTDIR/$BAM -insert > $OUTDIR/${BAM}.bamstat

## jitterbug find insertions
##
~/software/jitterbug/jitterbug.py \
--mem \
-l $SAMPLE \
-o $OUTDIR/$SAMPLE
-n $THREADS \
-b 50000000 \
-t $TAG \
$BAM \
$GFF

## jitterbug filter results
##
~/software/jitterbug/tools/jitterbug_filter_results_func.py \
-g $OUTDIR/${SAMPLE}.TE_insertions_paired_clusters.gff3 \
-c $OUTDIR/${SAMPLE}.filter_config.txt \
-o $OUTDIR/${SAMPLE}.TE_insertions_paired_clusters.filtered.gff3

## annotate genome for Ns
~/software/seqtk/seqtk \
cutN \
-gp10000000 \
-n1 $REF \
> $OUTDIR/${REF}.N_annot.bed

## remove insertions that span Ns
~/miniconda3/envs/jitterbug/bin/bedtools intersect \
-a $OUTDIR/${SAMPLE}.TE_insertions_paired_clusters.filtered.gff3 \
-b $OUTDIR/${REF}.N_annot.bed \
-v \
> $OUTDIR/${SAMPLE}.TE_insertions_paired_clusters.filtered.noNs.gff3

end=`date +%s`
runtime=$((end-start))
echo "End: "`date`
echo "Run time: $runtime"
echo "Done"

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR/$OUTDIR
