#!/bin/sh

#PBS -N raxml.Rg
#PBS -j oe
#PBS -k oe

#PBS -m ae
#PBS -M foo.bar@email.com

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=16:mem=100gb

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

## load common modules as standard
##
module load raxml/8.2.9

## RAxML commands, first build ML tree + 100 bootstraps, then draw bootstrap values on best tree
##
for file in *mafft;
  do echo $file;
  raxmlHPC-PTHREADS-AVX -f a -p 12345 -x 12345 -# 100 -m PROTGAMMAAUTO -T 16 -s ${file} -n ${file};
  raxmlHPC-PTHREADS-AVX -f b -t RAxML_bestTree.${file} -z RAxML_bootstrap.${file} -m PROTGAMMAAUTO -n RAxML_bestTree.${file}
done

mv RAxML* ../uniref90_trees/

## move LOGFILE to cwd
##
mv $HOME/$LOGFILE $PBS_O_WORKDIR
