#!/bin/bash

#PBS -N bowtie.Rg.uncoll
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M foo.bar@email.com

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=64:mem=100gb

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

## variables
source ~/virt_env/python/platypus/bin/activate
export LD_LIBRARY_PATH=/home/rnowell/software/htslib/lib

REF=nRg.uncollapsed_gapClosed.fa
READS1=~/data/magnacalcarata/raw/matepaired/Rmagna-MP_160629_L006_R1.fastq.gz
READS2=~/data/magnacalcarata/raw/matepaired/Rmagna-MP_160629_L006_R2.fastq.gz
#READSSE=mp_reads.se.fq.gz
OUTBAM=Rmagna-MP.raw.vs.nRc.uncollapsed_gapClosed
THREADS=64

## dump reads into single file:
#zcat $READS1 $READS2 > $READSSE
#echo Reads $READS1 and $READS2 dumped to file $READSSE

## index
samtools faidx $REF
/home/rnowell/software/bowtie2-2.2.6/bowtie2-build $REF $REF

## align
/home/rnowell/software/bowtie2-2.2.6/bowtie2 \
--threads 64 \
--end-to-end \
--very-sensitive \
--no-unal \
-x $REF \
-q \
-U $READS1 \
-U $READS2 \
| samtools view -@ $THREADS -b - \
| samtools sort -@ $THREADS -O bam -T temp_${RANDOM} - > ${OUTBAM}.bowtie2.sorted.bam

## index
samtools index ${OUTBAM}.bowtie2.sorted.bam

## bamstats
bamtools stats -in ${OUTBAM}.bowtie2.sorted.bam -insert > ${OUTBAM}.bowtie2.sorted.bam.bamstat

## get genomeCov
getseqlengths.pl -n $REF > ${REF}.genome.txt
~/software/bedtools2/bin/genomeCoverageBed \
-ibam ${OUTBAM}.bowtie2.sorted.bam \
-g ${REF}.genome.txt \
| perl -lne 'print if /^genome/' > ${OUTBAM}.bowtie2.sorted.bam.genomeCov

## platypus
~/software/Platypus/bin/Platypus.py callVariants \
--output=${OUTBAM}.bowtie2.sorted.bam.platypus.vcf \
--refFile=${REF} \
--bamFiles=${OUTBAM}.bowtie2.sorted.bam \
--minReads=5 \
--logFileName=${OUTBAM}.bowtie2.sorted.bam.platypus.log \
--nCPU=16 \
--minMapQual=30 \
--minBaseQual=20 \
--filterDuplicates=1

## get stats
~/software/vcflib/bin/vcfstats ${OUTBAM}.bowtie2.sorted.bam.platypus.vcf > ${OUTBAM}.bowtie2.sorted.bam.platypus.vcf.stats

## get freqs
perl -lane '
  BEGIN{ print join("\t","CHROM","POS","REF","ALT","FILTER","TC","TR","FREQREF","FREQALT","MAF") }
  if(length($F[3])==1 && length($F[4])==1){
      $TC = $1 if $F[7] =~ /TC\=(\d+)/;
      $TR = $1 if $F[7] =~ /TR\=(\d+)/;
      print join("\t",$F[0],$F[1],$F[3],$F[4],$F[6],$TC,$TR,(($TC-$TR)/$TC),($TR/$TC),(($TC-$TR)/$TC)<($TR/$TC) ? (($TC-$TR)/$TC) : ($TR/$TC));
  }' ${OUTBAM}.bowtie2.sorted.bam.platypus.vcf > ${OUTBAM}.bowtie2.sorted.bam.platypus.vcf.freqs

## delete bamfile to cleanup
rm $OUTBAM.bowtie2.sorted.bam ${OUTBAM}.bowtie2.sorted.bam.bai

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
