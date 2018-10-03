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


  echo "Choose an application (1-4):"
  echo "1. spray_insitu"
  echo "2. spray_ooc"
  echo "3. baseline_insitu"
  echo "4. baseline_ooc"
  
  read app
  if [ "$app" == "1" ]
  then
    NUM_MPI_TASKS=1
    NUM_THREADS=1
    COMMAND="$SPRAY_BIN/spray_insitu --ply-path $PLY_PATH --nthreads $NUM_THREADS -w 512 -h 512 --frames -1 --mode glfw --cache-size -1 --partition insitu --camera 90.172180 84.141418 82.480225 30.000000 28.649426 30.000000 --pixel-samples 1 --ao-samples 1 --bounces 1 --num-partitions 8 --blinn 0.4 0.4 0.4 10 $MODEL"
  elif [ "$app" == "4" ]
  then
    NUM_MPI_TASKS=1
    NUM_THREADS=4
    COMMAND="$SPRAY_BIN/baseline_ooc --ply-path $PLY_PATH --nthreads $NUM_THREADS -w 512 -h 512 --frames -1 --mode glfw --cache-size -1 --partition image --camera 90.172180 84.141418 82.480225 30.000000 28.649426 30.000000 --pixel-samples 1 --ao-samples 1 --bounces 1 --blinn 0.4 0.4 0.4 10 $MODEL"
  else
    NUM_MPI_TASKS=1
    NUM_THREADS=1
    COMMAND=""
    echo "[error] invalid input"
  fi

  echo "NUM_MPI_TASKS=$NUM_MPI_TASKS"
  echo "OMP_NUM_THREADS=$NUM_THREADS"
  echo $COMMAND

  export OMP_NUM_THREADS=$NUM_THREADS
  $COMMAND
fi
