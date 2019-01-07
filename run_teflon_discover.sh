#!/bin/bash

#PBS -N teflon_discover
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
conda activate teflon

## load java
source ~/.jdk/8

## declare variables
REFNAME=reference_RM15 ## edit me!
WORKDIR=teflon_workdir/ ## should be same as in teflon_prepare.sh
declare -a SAMPLES=("MAG1" "MAG2" "MAG3" "RM9") ## add sample names here
## leave alone
REF=$WORKDIR/${REFNAME}.prep_MP/${REFNAME}.mappingRef.fa

start=`date +%s`

echo "Start: "`date`
echo "VARIABLES"
echo "Working directory: $WORKDIR"
echo "Refence: $REF"
echo "Sample names: ${SAMPLES[@]}"
echo

## index reference if needed
if [ ! -f ${REF}.sa ]; then
  bwa index $REF
fi

## print parallel commands to logfile
parallel --dry-run "mkdir teflon_{}" ::: ${SAMPLES[@]}
## ensure paths and extensions to fq files are correct!
parallel --dry-run "bwa mem -t 8 -Y $REF ../../../{}/{}_1.trimmed.ecc.filtered.fq.gz ../../../{}/{}_2.trimmed.ecc.filtered.fq.gz | samtools view -@ 8 -b - | samtools sort -@ 8 -O bam -T teflon_{}/{}.tmp - > teflon_{}/{}.sorted.bam" ::: ${SAMPLES[@]}
parallel --dry-run "samtools stats teflon_{}/{}.sorted.bam > teflon_{}/{}.sorted.bam.samstats" ::: ${SAMPLES[@]}
parallel --dry-run "bamtools stats -in teflon_{}/{}.sorted.bam -insert > teflon_{}/{}.sorted.bam.bamstats" ::: ${SAMPLES[@]}

## then do mapping
parallel "mkdir teflon_{}" ::: ${SAMPLES[@]}
## ensure paths and extensions to fq files are correct!
parallel "bwa mem -t 8 -Y $REF ../../../{}/{}_1.trimmed.ecc.filtered.fq.gz ../../../{}/{}_2.trimmed.ecc.filtered.fq.gz | samtools view -@ 8 -b - | samtools sort -@ 8 -O bam -T teflon_{}/{}.tmp - > teflon_{}/{}.sorted.bam" ::: ${SAMPLES[@]}
parallel "samtools index teflon_{}/{}.sorted.bam" ::: ${SAMPLES[@]}
parallel "samtools stats teflon_{}/{}.sorted.bam > teflon_{}/{}.sorted.bam.samstats" ::: ${SAMPLES[@]}
parallel "bamtools stats -in teflon_{}/{}.sorted.bam -insert > teflon_{}/{}.sorted.bam.bamstats" ::: ${SAMPLES[@]}

## delete the old file if exists
if [ -f samples.txt ]; then
    rm samples.txt
fi
touch samples.txt

## write samples.txt and TEFLoN shells
##
for i in "${SAMPLES[@]}"; do
  ## write samples
  printf "$PWD/teflon_${i}/${i}.sorted.bam\t${i}\n" >> $WORKDIR/samples.txt;
  ## get SD of insert for overwrite
  SDEV="$(grep -Poe '(?<=insert size standard deviation:\s)(\d+)' teflon_${i}/${i}.sorted.bam.samstats)";
  ## print shell commands to file
  shell="
  python ~/software/TEFLoN/teflon.v0.4.py \
  -wd $PWD/$WORKDIR \
  -d $PWD/$WORKDIR/${REFNAME}.prep_TF/ \
  -s $PWD/$WORKDIR/samples.txt \
  -i $i \
  -eb ~/miniconda3/envs/teflon/bin/bwa \
  -es ~/miniconda3/envs/teflon/bin/samtools \
  -l1 family \
  -l2 family \
  -q 30 \
  -sd $SDEV \
  -t 8
  ";
  echo $shell > teflon_${i}/run_teflon.sh;
done
## execute them
parallel --dry-run "sh teflon_{}/run_teflon.sh" ::: ${SAMPLES[@]}
parallel "sh teflon_{}/run_teflon.sh" ::: ${SAMPLES[@]}

## run TEFLoN collapse on all samples
##
python ~/software/TEFLoN/teflon_collapse.py \
-wd $PWD/$WORKDIR \
-d $PWD/$PREPTF_DIR \
-s $PWD/samples.txt \
-es ~/miniconda3/envs/teflon/bin/samtools \
-n1 1 \
-n2 1 \
-q 30 \
-t 8

## run TEFLoN count for each sample
##
for i in "${SAMPLES[@]}"; do
  echo "Running teflon_count.py on sample ${i}...";
  python ~/software/TEFLoN/teflon_count.py \
  -wd $PWD/$WORKDIR \
  -d $PWD/$WORKDIR/${REFNAME}.prep_TF/ \
  -s $PWD/$WORKDIR/samples.txt \
  -i $i \
  -eb ~/miniconda3/envs/teflon/bin/bwa \
  -es ~/miniconda3/envs/teflon/bin/samtools \
  -l2 family \
  -q 30 \
  -t 8;
  echo "Done ${i}!";
done

## run TEFLoN genotype once
##
python ~/software/TEFLoN/teflon_genotype.py \
-wd $PWD/$WORKDIR \
-d $PWD/$WORKDIR/${REFNAME}.prep_TF/ \
-s $PWD/$WORKDIR/samples.txt \
-dt pooled

## finish
end=`date +%s`
runtime=$((end-start))
echo "End: "`date`
echo "Run time: $runtime"
echo "Done"

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR/$WORKDIR
