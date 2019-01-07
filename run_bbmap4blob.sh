#!/bin/bash

#PBS -N bbmap
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M reubennowell@gmail.com

#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=8:mem=50gb

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

## load bbmap environment
source /home/rnowell/.jdk/8
export LD_LIBRARY_PATH=/home/rnowell/libraries/zlib/zlib-1.2.8/lib/:${LD_LIBRARY_PATH}

## vars
REF=
IN1=
IN2=
OUTM=
THREADS=8

ln -s $IN1 pe_1.fq.gz
ln -s $IN2 pe_2.fq.gz

## bbmap
~/software/bbmap-37.82/bbmap.sh \
-Xmx450g \
nodisk=t \
ref=${REF} \
in1=pe_1.fq.gz \
in2=pe_2.fq.gz \
threads=${THREADS} \
outm=${OUTM}.sam \
covstats=${OUTM}.sam.covstats \
covhist=${OUTM}.sam.covhist

## generate covfile for blobtools
perl -lane 'BEGIN{print "# contig_id\tread_cov\tbase_cov"}if(/#/){next}else{$reads=($F[6]+$F[7]);print join("\t",$F[0],$reads,$F[1])}' ${OUTM}.sam.covstats > ${OUTM}.sam.cov

## convert to bam
samtools sort -O bam -T temp -@ $THREADS ${OUTM}.sam > ${OUTM}.bam
samtools sort -O bam -T temp -@ $THREADS -n ${OUTM}.sam > readsorted.bam
rm ${OUTM}.sam

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
