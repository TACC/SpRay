#!/bin/bash

if [ -z $SPRAY_HOME_PATH ]
then
  echo "[error] SPRAY_HOME_PATH not found. Do export SPRAY_HOME_PATH=<path_to_spray_home>."

else
  SPRAY_BIN=$SPRAY_HOME_PATH/build
  EXAMPLE_PATH=$SPRAY_HOME_PATH/examples
  WAVELET16_PATH=$SPRAY_HOME_PATH/examples/wavelets64
  MODEL=$WAVELET16_PATH/wavelets64.domain
  PLY_PATH=$EXAMPLE_PATH
  
  mpirun -n 2 $SPRAY_BIN/spray_insitu --ply-path $PLY_PATH --nthreads 1 -w 512 -h 512 --frames 1 --mode film --cache-size -1 --partition insitu --camera 90.172180 84.141418 82.480225 30.000000 28.649426 30.000000 --pixel-samples 1 --ao-samples 1 --bounces 1 --num-partitions 8 --shading blinn --blinn 0.4 0.4 0.4 10 $MODEL

fi
