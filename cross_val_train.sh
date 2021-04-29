#!/bin/bash
#PBS -q v100
#PBS -l nodes=1:ppn=36
#PBS -l walltime=72:00:00
#PBS -A hpc_michal01
#PBS -j oe

cd /work/derick/siamese-monet-project/Pocket2Drug/

singularity exec --nv -B /work,/project,/usr/lib64 /work/derick/singularities/pytorch171.simg python cross_validation_train.py -val_fold 0 -out_dir ../p2d_results/cross_val_fold_0/ &> ../p2d_results/cross_val_fold_0/train.log 2>&1

singularity exec --nv -B /work,/project,/usr/lib64 /work/derick/singularities/pytorch171.simg python cross_validation_train.py -val_fold 1 -out_dir ../p2d_results/cross_val_fold_1/ &> ../p2d_results/cross_val_fold_1/train.log 2>&1

singularity exec --nv -B /work,/project,/usr/lib64 /work/derick/singularities/pytorch171.simg python cross_validation_train.py -val_fold 2 -out_dir ../p2d_results/cross_val_fold_2/ &> ../p2d_results/cross_val_fold_2/train.log 2>&1

singularity exec --nv -B /work,/project,/usr/lib64 /work/derick/singularities/pytorch171.simg python cross_validation_train.py -val_fold 3 -out_dir ../p2d_results/cross_val_fold_3/ &> ../p2d_results/cross_val_fold_3/train.log 2>&1

singularity exec --nv -B /work,/project,/usr/lib64 /work/derick/singularities/pytorch171.simg python cross_validation_train.py -val_fold 4 -out_dir ../p2d_results/cross_val_fold_4/ &> ../p2d_results/cross_val_fold_4/train.log 2>&1

