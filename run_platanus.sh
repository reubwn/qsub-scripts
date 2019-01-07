#!/bin/bash

#PBS -N platanus
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M reubennowell@gmail.com

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=64:mem=400gb

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

######################## platanus

## define reads
PE1=
PE2=
MP1=
MP2=
PREFIX=
THREADS=

## inflate reads into cwd
printf "\nInflating readfiles into CWD... "
zcat ${PE1} > pe_1.fq
zcat ${PE2} > pe_2.fq
zcat ${MP1} > mp_1.fq
zcat ${MP2} > mp_2.fq
printf "done\n"

/home/rnowell/bin/platanus assemble \
-o ${PREFIX} \
-f pe_1.fq pe_2.fq \
-u 0.1 \
-m 400 \
-t $THREADS

/home/rnowell/bin/platanus scaffold \
-o ${PREFIX} \
-c ${PREFIX}_contig.fa \
-b ${PREFIX}_contigBubble.fa \
-IP1 pe_1.fq pe_2.fq \
-OP2 mp_1.fq mp_2.fq \
-u 0.1 \
-t $THREADS

/home/rnowell/bin/platanus gap_close \
-o ${PREFIX} \
-c ${PREFIX}_scaffold.fa \
-IP1 pe_1.fq pe_2.fq \
-t $THREADS

printf "Deleting readfiles... "
rm pe_1.fq pe_2.fq mp_1.fq mp_2.fq
printf "done\n"

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
