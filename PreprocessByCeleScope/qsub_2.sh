#!/bin/bash
#-V:pass all environment variables to the job
#-N: job name
#-l: specify the amount of maximum memory required
#-d: working directory
work_dir=$PWD #default current working directory
echo $work_dir
#qsub -V -N qTest -l h_vmem=5G -d $work_dir test.sh
#qsub -V -N mt -l mem=10000MB -v "arg1=1000,arg2=444" -d $work_dir run1.sh
#qsub -V -N qs -l mem=10000MB,nodes=cu02:ppn=1 -d $work_dir run.sh
qsub -V -N qSCY_6M -l mem=40000MB,nodes=cu03 -d $work_dir run_SC-Y-6M.sh
qsub -V -N qSCY_4M -l mem=40000MB,nodes=cu02 -d $work_dir run_SC-Y-4M.sh
qsub -V -N qSCO_23M -l mem=40000MB,nodes=cu01 -d $work_dir run_SC-O-23M.sh
