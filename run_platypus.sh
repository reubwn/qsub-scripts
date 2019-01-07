#!/bin/sh

#PBS -N platypus
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M foo.bar@email.com

#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=16:mem=10gb

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

## virt env
source ~/virt_env/python/platypus/bin/activate
export LD_LIBRARY_PATH=/home/rnowell/software/htslib/lib

## files
REF=
BAM=

## platypus
##
~/software/Platypus/bin/Platypus.py callVariants \
--output=${BAM}.platypus.vcf \
--refFile=${REF} \
--bamFiles=${BAM} \
--minReads=5 \
--maxReads=1000 \
--logFileName=${BAM}.platypus.log \
--nCPU=16 \
--minMapQual=30 \
--minBaseQual=20 \
--filterDuplicates=1

## get stats
~/software/vcflib/bin/vcfstats ${BAM}.platypus.vcf > ${BAM}.platypus.vcf.stats

## get freqs
perl -lane '
  BEGIN{ print join("\t","CHROM","POS","REF","ALT","FILTER","TC","TR","FREQREF","FREQALT","MAF") }
  if(length($F[3])==1 && length($F[4])==1){
      $TC = $1 if $F[7] =~ /TC\=(\d+)/;
      $TR = $1 if $F[7] =~ /TR\=(\d+)/;
      print join("\t",$F[0],$F[1],$F[3],$F[4],$F[6],$TC,$TR,(($TC-$TR)/$TC),($TR/$TC),(($TC-$TR)/$TC)<($TR/$TC) ? (($TC-$TR)/$TC) : ($TR/$TC));
  }' ${BAM}.platypus.vcf > ${BAM}.platypus.vcf.freqs

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
