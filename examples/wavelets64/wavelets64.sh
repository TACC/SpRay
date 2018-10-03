#!/bin/bash

if [ -z $SPRAY_HOME_PATH ]
then
  echo "[error] SPRAY_HOME_PATH not found. Do export SPRAY_HOME_PATH=<path_to_spray_home>."
  return
fi

MODE=$1

if [ $MODE != "film" ] && [ $MODE != "glfw" ] 
then
  echo "[error] invalid mode: $MODE"
  echo "[syntax] wavelets64.sh MODE"
  echo "         MODE = {film, glfw}"
  return
fi

NUM_FRAMES="-1"
if [ $MODE == "film" ]
then
  NUM_FRAMES="1"
fi

SPRAY_BIN=$SPRAY_HOME_PATH/build
EXAMPLE_PATH=$SPRAY_HOME_PATH/examples
WAVELET16_PATH=$SPRAY_HOME_PATH/examples/wavelets64
MODEL=$WAVELET16_PATH/wavelets64.domain
PLY_PATH=$EXAMPLE_PATH
MPI_BIN="mpirun -n"

echo "Choose an application (1-4):"
echo "1. spray_insitu"
echo "2. spray_ooc"
echo "3. baseline_insitu"
echo "4. baseline_ooc"

read APP
if [ $APP == "1" ]
then
  NUM_MPI_TASKS=2
  NUM_THREADS=1 # threading not supported
  COMMAND="$MPI_BIN $NUM_MPI_TASKS $SPRAY_BIN/spray_insitu --ply-path $PLY_PATH --nthreads 1 -w 512 -h 512 --frames $NUM_FRAMES --mode $MODE --cache-size -1 --partition insitu --camera 90.172180 84.141418 82.480225 30.000000 28.649426 30.000000 --pixel-samples 1 --ao-samples 1 --bounces 1 --num-partitions 8 --blinn 0.4 0.4 0.4 10 $MODEL"
elif [ $APP == "4" ]
then
  NUM_MPI_TASKS=2
  NUM_THREADS=2
  COMMAND="$MPI_BIN $NUM_MPI_TASKS $SPRAY_BIN/baseline_ooc --ply-path $PLY_PATH --nthreads $NUM_THREADS -w 512 -h 512 --frames $NUM_FRAMES --mode $MODE --cache-size -1 --partition image --camera 90.172180 84.141418 82.480225 30.000000 28.649426 30.000000 --pixel-samples 1 --ao-samples 1 --bounces 1 --blinn 0.4 0.4 0.4 10 $MODEL"
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

