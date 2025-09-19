#!/bin/bash

python reaction_prediction/atom/scripts/prepare_all_feat_ids.py \
  --input data/mc_train_fold0/with_quotes/train_val_combined.txt \
  --output output/mc_train_fold0/train_val_combined_allids.json \
  --length 3