#!/bin/bash

#PBS -N teflon
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M r.nowell@imperial.ac.uk

#PBS -l walltime=24:00:00
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
conda activate teflon

## load java
source ~/.jdk/8

## declare variables
TE=~/medbio/software/RepeatMasker/Libraries/Repbase_custom/metazoa_repeatlib.fa
REF=../../scaffolds_final.AK11.fasta
PE1=../../../AK15/AK15_1.trimmed.ecc.filtered.fq.gz
PE2=../../../AK15/AK15_2.trimmed.ecc.filtered.fq.gz
BAM=sorted.bam
SAMPLE=AK15
OUTDIR=teflon_${SAMPLE}
mkdir $OUTDIR
ln -s $TE $OUTDIR/
THREADS=8

start=`date +%s`

echo "Start: "`date`
echo
echo "VARIABLES"
echo "Refence: $REF"
echo "PE1: $PE1"
echo "PE2: $PE2"
echo "Sample: $SAMPLE"
echo "Threads: $THREADS"
echo

## prepare genome
python ~/software/TEFLoN/teflon_prep_custom.py \
-wd $PWD/$OUTDIR \
-e ~/miniconda3/envs/teflon/bin/RepeatMasker \
-g $REF \
-l $OUTDIR/metazoa_repeatlib.fa \
-p $SAMPLE \
-t $THREADS

## make te_annot bed file
perl -F"\t" -lane '@a=split(/#/,$F[3]);print join("\t",$F[0],$F[1],$F[2],$F[5],$F[3],$a[1],$a[2])' $PWD/$OUTDIR/${SAMPLE}.prep_RM/${SAMPLE}.bed > $PWD/$OUTDIR/te_annot.bed

## index reference
bwa index $OUTDIR/${SAMPLE}.prep_MP/${SAMPLE}.mappingRef.fa

## run bwa mem
bwa mem \
-t $THREADS \
-Y \
$OUTDIR/${SAMPLE}.prep_MP/${SAMPLE}.mappingRef.fa \
$PE1 \
$PE2 \
| samtools view -@ $THREADS -b - \
| samtools sort -@ $THREADS -O bam -T temp_${RANDOM} - > $OUTDIR/$BAM

### index bam
samtools index $OUTDIR/sorted.bam

### gather some stats
bamtools stats -in $OUTDIR/$BAM -insert > $OUTDIR/${BAM}.bamstat

## generate config.txt
printf "$PWD/$OUTDIR/sorted.bam\t$SAMPLE" > $OUTDIR/samples.txt

## get insert SD manually and save to $SD
SDEV="$(samtools stats -t $PWD/$OUTDIR/AK15.prep_TF/AK15.genomeSize.txt $PWD/$OUTDIR/sorted.bam | grep -Poe '(?<=insert size standard deviation:\s)(\d+)')"
printf "\nInsert SD is found to be: $SDEV\n\n"

## TEFLoN command
##
python ~/software/TEFLoN/teflon.v0.4.py \
-wd $PWD/$OUTDIR \
-d $PWD/$OUTDIR/${SAMPLE}.prep_TF/ \
-s $PWD/$OUTDIR/samples.txt \
-i $SAMPLE \
-eb ~/miniconda3/envs/teflon/bin/bwa \
-es ~/miniconda3/envs/teflon/bin/samtools \
-l1 family \
-l2 family \
-q 30 \
-sd $SDEV \
-t $THREADS

end=`date +%s`
runtime=$((end-start))
echo "End: "`date`
echo "Run time: $runtime"
echo "Done"

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR/$OUTDIR
