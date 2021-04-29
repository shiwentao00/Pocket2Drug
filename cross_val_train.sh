#!/bin/bash
#PBS -q v100
#PBS -l nodes=1:ppn=36
#PBS -l walltime=72:00:00
#PBS -A hpc_michal01
#PBS -j oe

cd /work/derick/siamese-monet-project/Pocket2Drug/

singularity exec --nv -B /work,/project,/usr/lib64 /work/derick/singularities/pytorch171.simg python cross_val_train.py -val_fold 0 -out_dir ../p2d_results/cross_val_fold_0/ &> ./cross_validation/train_logs/train_fold0.log 2>&1
singularity exec --nv -B /work,/project,/usr/lib64 /work/derick/singularities/pytorch171.simg python cross_val_train.py -val_fold 1 -out_dir ../p2d_results/cross_val_fold_1/ &> ./cross_validation/train_logs/train_fold1.log 2>&1
singularity exec --nv -B /work,/project,/usr/lib64 /work/derick/singularities/pytorch171.simg python cross_val_train.py -val_fold 2 -out_dir ../p2d_results/cross_val_fold_2/ &> ./cross_validation/train_logs/train_fold2.log 2>&1
singularity exec --nv -B /work,/project,/usr/lib64 /work/derick/singularities/pytorch171.simg python cross_val_train.py -val_fold 3 -out_dir ../p2d_results/cross_val_fold_3/ &> ./cross_validation/train_logs/train_fold3.log 2>&1
singularity exec --nv -B /work,/project,/usr/lib64 /work/derick/singularities/pytorch171.simg python cross_val_train.py -val_fold 4 -out_dir ../p2d_results/cross_val_fold_4/ &> ./cross_validation/train_logs/train_fold4.log 2>&1

