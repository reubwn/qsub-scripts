#!/bin/bash

#PBS -N tepid
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M foo.bar@email.com

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:mem=50gb

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
conda activate TEPID

## load java
source ~/.jdk/8

## declare variables
REF=../../scaffolds_final.AK16.fasta
declare -a SAMPLES=("AK11" "AK15" "AK27" "RS1")
declare -a INSERTS=("140" "132" "140" "141") ## in same order as above!
FQ_DIR=fq_dir
ANNOT=/home/rnowell/data/socialis/filtered/wga/AK16/te_analysis/1.teflon/teflon_workdir/te_annot.TEPID.txt ##from TEFLoN

start=`date +%s`

echo "Start: "`date`
echo
echo "VARIABLES"
echo
echo "Refence: $REF"
echo "Sample names: ${SAMPLES[@]}"
echo "Insert sizes: ${INSERTS[@]}"
echo

## setup
#mkdir fq_dir
ln -s $PWD/$REF reference.fa
ln -s $ANNOT

## bowtie2 index
bowtie2-build reference.fa reference.fa

## yaha index
yaha -g reference.fa

## tepid-map
##
for (( i=0; i<${#SAMPLES[@]}; i++ )); do
  echo "Working on ${SAMPLES[$i]}";
  echo "Insert size ${INSERTS[$i]}";
  ## link to fastq files
  ln -s $PWD/../../../${SAMPLES[$i]}/${SAMPLES[$i]}_1.trimmed.ecc.filtered.fq.gz ${SAMPLES[$i]}_1.fastq.gz;
  ln -s $PWD/../../../${SAMPLES[$i]}/${SAMPLES[$i]}_2.trimmed.ecc.filtered.fq.gz ${SAMPLES[$i]}_2.fastq.gz;
  ## write shell
  shell="
  tepid-map \
  -x reference.fa \
  -y reference.X15_01_65525S \
  -p 8 \
  -s ${INSERTS[$i]} \
  -n ${SAMPLES[$i]} \
  -1 ${SAMPLES[$i]}_1.fastq.gz \
  -2 ${SAMPLES[$i]}_2.fastq.gz \
  -z
  ";
  echo $shell > map_${SAMPLES[$i]}.sh;
done
## execute them
parallel --dry-run "sh map_{}.sh" ::: ${SAMPLES[@]}
parallel "sh map_{}.sh" ::: ${SAMPLES[@]}

echo "#######"
echo "Mapping finished "`date`
echo "#######"

## tepid-discover
##
parallel --dry-run "tepid-discover -p 8 -n {} -c {}.bam -s {}.split.bam -t $ANNOT" ::: ${SAMPLES[@]}
parallel "tepid-discover -p 8 -n {} -c {}.bam -s {}.split.bam -t $ANNOT" ::: ${SAMPLES[@]}

## finish
end=`date +%s`
runtime=$((end-start))
echo "End: "`date`
echo "Run time: $runtime"
echo "Done"

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
