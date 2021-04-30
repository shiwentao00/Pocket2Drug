#!/bin/bash
#PBS -q v100
#PBS -l nodes=1:ppn=36
#PBS -l walltime=72:00:00
#PBS -A hpc_michal01
#PBS -j oe

cd /work/derick/siamese-monet-project/Pocket2Drug/

singularity exec --nv -B /work,/project,/usr/lib64 /work/derick/singularities/pytorch171.simg python cross_val_sample.py -fold 5 -result_dir ../p2d_results/cross_val_fold_5/ &> ./cross_validation/sample_logs/sample_fold5.log 2>&1
singularity exec --nv -B /work,/project,/usr/lib64 /work/derick/singularities/pytorch171.simg python cross_val_sample.py -fold 6 -result_dir ../p2d_results/cross_val_fold_6/ &> ./cross_validation/sample_logs/sample_fold6.log 2>&1
singularity exec --nv -B /work,/project,/usr/lib64 /work/derick/singularities/pytorch171.simg python cross_val_sample.py -fold 7 -result_dir ../p2d_results/cross_val_fold_7/ &> ./cross_validation/sample_logs/sample_fold7.log 2>&1
singularity exec --nv -B /work,/project,/usr/lib64 /work/derick/singularities/pytorch171.simg python cross_val_sample.py -fold 8 -result_dir ../p2d_results/cross_val_fold_8/ &> ./cross_validation/sample_logs/sample_fold8.log 2>&1
singularity exec --nv -B /work,/project,/usr/lib64 /work/derick/singularities/pytorch171.simg python cross_val_sample.py -fold 9 -result_dir ../p2d_results/cross_val_fold_9/ &> ./cross_validation/sample_logs/sample_fold9.log 2>&1

