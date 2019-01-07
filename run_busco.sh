#!/bin/bash

#PBS -N busco
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M reubennowell@gmail.com

#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=8:mem=5gb

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

## load perlbrew environment
source /home/rnowell/perl5/perlbrew/etc/bashrc

## load necessities
export PATH=/home/rnowell/software/augustus-3.2.1/bin/:${PATH}
export PATH=/home/rnowell/software/augustus-3.2.1/scripts/:${PATH}
export PATH=/home/rnowell/software/blat/:${PATH}
export AUGUSTUS_CONFIG_PATH=/home/rnowell/software/augustus-3.2.1/config/
export BAMTOOLS_PATH=/home/rnowell/libraries/bamtools/bin/
export GENEMARK_PATH=/home/rnowell/software/genemark-es-4.32/
export SAMTOOLS_PATH=/home/rnowell/software/samtools-1.2/samtools
export LD_LIBRARY_PATH=/home/rnowell/libraries/bamtools/lib/:${LD_LIBRARY_PATH}

## env
source ~/virt_env/python/busco/bin/activate

## ln to lineage
ln -s ~/software/busco/eukaryota_odb9
ln -s ~/software/busco/metazoa_odb9

## files
GENOME=
THREADS=8
LINEAGE=eukaryota_odb9/

## busco
~/software/busco/scripts/run_BUSCO.py \
-i $GENOME \
-o eukaryota_odb9.${GENOME} \
-l $LINEAGE \
-m genome \
-c $THREADS \
--limit 4 \
-t ./tmp.${GENOME}

## rm tmp
rm -rf ./tmp.${GENOME}

LINEAGE=metazoa_odb9/

## busco
~/software/busco/scripts/run_BUSCO.py \
-i $GENOME \
-o metazoa_odb9.${GENOME} \
-l $LINEAGE \
-m genome \
-c $THREADS \
--limit 4 \
-t ./tmp.${GENOME}

## rm tmp
rm -rf ./tmp.${GENOME}

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
