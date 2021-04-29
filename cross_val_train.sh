#!/bin/bash
#PBS -q v100
#PBS -l nodes=1:ppn=36
#PBS -l walltime=72:00:00
#PBS -A hpc_michal01
#PBS -j oe

cd /work/derick/siamese-monet-project/Pocket2Drug/

singularity exec --nv -B /work,/project,/usr/lib64 /work/derick/singularities/pytorch171.simg python cross_val_train.py -val_fold 5 -out_dir ../p2d_results/cross_val_fold_5/ &> ./cross_validation/train_logs/train_fold5.log 2>&1
singularity exec --nv -B /work,/project,/usr/lib64 /work/derick/singularities/pytorch171.simg python cross_val_train.py -val_fold 6 -out_dir ../p2d_results/cross_val_fold_6/ &> ./cross_validation/train_logs/train_fold6.log 2>&1
singularity exec --nv -B /work,/project,/usr/lib64 /work/derick/singularities/pytorch171.simg python cross_val_train.py -val_fold 7 -out_dir ../p2d_results/cross_val_fold_7/ &> ./cross_validation/train_logs/train_fold7.log 2>&1
singularity exec --nv -B /work,/project,/usr/lib64 /work/derick/singularities/pytorch171.simg python cross_val_train.py -val_fold 8 -out_dir ../p2d_results/cross_val_fold_8/ &> ./cross_validation/train_logs/train_fold8.log 2>&1
singularity exec --nv -B /work,/project,/usr/lib64 /work/derick/singularities/pytorch171.simg python cross_val_train.py -val_fold 9 -out_dir ../p2d_results/cross_val_fold_9/ &> ./cross_validation/train_logs/train_fold9.log 2>&1

