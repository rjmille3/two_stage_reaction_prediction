#!/bin/bash

python reaction_prediction/atom/scripts/train_atom_model.py \
  --train_h5 output/mc_train_fold0/atom_training/source/train.hdf5 \
  --val_h5 output/mc_train_fold0/atom_training/source/val.hdf5 \
  --model_out output/mc_train_fold0/models/atom/source.h5 \
  --plot_dir output/mc_train_fold0/models/atom/source_plots \
  --config model_configs/atom_config.json