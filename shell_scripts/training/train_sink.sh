#!/bin/bash

python reaction_prediction/atom/scripts/train_atom_model.py \
  --train_h5 output/mc_train_fold0/atom_training/sink/train.hdf5 \
  --val_h5 output/mc_train_fold0/atom_training/sink/val.hdf5 \
  --model_out output/mc_train_fold0/models/atom/sink.h5 \
  --plot_dir output/mc_train_fold0/models/atom/sink_plots \
  --config model_configs/atom_config.json