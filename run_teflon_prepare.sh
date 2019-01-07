#!/bin/bash

#PBS -N teflon_prep
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M r.nowell@imperial.ac.uk

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=24:mem=20gb

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
TE=~/medbio/software/RepeatMasker/Libraries/Repbase_custom/te_lib.metazoa.all_TEs.fa
REF=../../scaffolds_final.AK16.fasta
SAMPLE=reference_AK16
OUTDIR=teflon_${SAMPLE}
mkdir $OUTDIR
ln -s $TE $OUTDIR/repeats.fa
THREADS=24

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
-g $PWD/$REF \
-l $PWD/$OUTDIR/repeats.fa \
-p $SAMPLE \
-t $THREADS

## make te_annot file for TEPID
perl -F"\t" -lane '@a=split(/#/,$F[3]);print join("\t",$F[0],$F[1],$F[2],$F[5],$F[3],$a[1],$a[2])' $PWD/$OUTDIR/${SAMPLE}.prep_RM/${SAMPLE}.bed > $PWD/$OUTDIR/te_annot.TEPID.txt

## convert BED to GFF3
perl -F"\t" -lane 'BEGIN{print "##gff-version 3"}{print join("\t",$F[0],"teflon",$F[3],($F[1]+1),$F[2],$F[4],$F[5],".",".")}' $PWD/$OUTDIR/${SAMPLE}.prep_RM/${SAMPLE}.bed > $PWD/$OUTDIR/${SAMPLE}.prep_RM/${SAMPLE}.gff3

## link to useful files
ln -s $PWD/$OUTDIR/${SAMPLE}.prep_RM/${SAMPLE}.bed $PWD/$OUTDIR/
ln -s $PWD/$OUTDIR/${SAMPLE}.prep_RM/${SAMPLE}.gff3 $PWD/$OUTDIR/
ln -s $PWD/$OUTDIR/${SAMPLE}.prep_TF/${SAMPLE}.hier $PWD/$OUTDIR/

end=`date +%s`
runtime=$((end-start))
echo "End: "`date`
echo "Run time: $runtime"
echo "Done"

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR/$OUTDIR
