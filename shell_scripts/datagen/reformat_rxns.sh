#!/bin/bash

python reaction_prediction/util/write_full_reaction_polar.py \
  --input data/mc_train_fold0/with_quotes/test.txt \
  --output data/mc_train_fold0/reformatted/test.txt