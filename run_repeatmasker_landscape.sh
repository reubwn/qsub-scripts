#!/bin/bash

#PBS -N rmscape
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M foo.bar@email.com

#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=16:mem=25gb

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

THREADS=16
LIB=~/medbio/software/RepeatMasker/Libraries/Repbase_custom/te_lib.metazoa.all_TEs.fa
FASTA=
ln -s $LIB
ln -s $FASTA

## get basename
BASE=`basename $FASTA`
printf "\nFile is: $BASE\n\n"

## hard mask repeatlib
/med-bio/rnowell/software/RepeatMasker/RepeatMasker \
-pa $THREADS \
-lib $LIB \
-dir . \
-alignments \
-gff \
-no_is \
-nolow \
-norna \
$FASTA

## convert .out to GFF3
perl -lne '
  BEGIN{
    print "##gff-version 3";
  };
  if($.<4){
    next;
  }else{
    $_=~s/^\s*(.*?)\s*$/$1/;
    @F=split(m/\s+/,$_);
    @a=split(m/\//,$F[10]);
    print join("\t",
              $F[4],
              "RepeatMasker",
              "repeat",
              $F[5],
              $F[6],
              "$F[1]",
              ($F[8]eq"C"?"-":"+"),
              ".",
              (join(";",
                    "ID=te$F[14]",
                    "Name=$F[9]",
                    (exists($a[1])?"Family=$a[1]":"Family=Unspecified"),
                    (exists($a[0])?"Class=$a[0]":"Class=Unspecified")
                   )
              ));
  }' ${BASE}.out > ${BASE}.out.gff3

## create divergence summary
~/medbio/software/RepeatMasker/util/calcDivergenceFromAlign.pl \
-s ${BASE}.divsum \
-a ${BASE}.new_align \
${BASE}.align

## get genome length
LENGTH=`perl -lne 'if(/>/){next}else{chomp;$total+=length}END{print $total}' $FASTA`

## create landscape
~/medbio/software/RepeatMasker/util/createRepeatLandscape.pl \
-div ${BASE}.divsum \
-g $LENGTH > ${BASE}.repeatlandscape.html

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
