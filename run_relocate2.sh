#!/bin/bash

#PBS -N rel2_bwt
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M r.nowell@imperial.ac.uk

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=24:mem=18gb

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
conda activate RelocaTE2
echo Using conda in `which conda`

## variables
TE=~/medbio/software/RepeatMasker/Libraries/Repbase_custom/metazoa_repeatlib.fa
REF=
PE1=
PE2=
RMOUT=
INSERT=
SAMPLE=RelocaTE2
THREADS=24

start=`date +%s`

## echo
echo
echo "$start"
echo
echo "VARIABLES"
echo "Refence: $REF"
echo "TE fasta: $TE"
echo "PE1: $PE1"
echo "PE2: $PE2"
echo "Insert size: $INSERT"
echo "RepeatMasker *.out: $RMOUT"
echo "Sample name: $SAMPLE"
echo

## ln
ln -s $PWD/$REF && REF=`basename $REF`
ln -s $PWD/$RMOUT && RMOUT=`basename $RMOUT`
ln -s $TE && TE=`basename $TE`
mkdir fq_dir
ln -s $PWD/$PE1 fq_dir/reads_1.fastq.gz
ln -s $PWD/$PE2 fq_dir/reads_2.fastq.gz

## RelocaTE2 command
##
relocaTE2.py \
-t $TE \
-d fq_dir/ \
-g $REF \
-r $RMOUT \
-o results_bowtie2 \
-s $INSERT \
--cpu $THREADS \
--sample $SAMPLE \
--aligner bowtie2 \
--len_cut_match 20 \
--len_cut_trim 20 \
--mismatch 0 \
--run \
--split \
--verbose 4

end=`date +%s`
runtime=$((end-start))
echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"
echo "Done"

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
