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

SPRAY_BIN_PATH=$SPRAY_HOME_PATH/build
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
  SPRAY_BIN=spray_insitu
  PARTITION=insitu

elif [ $APP == "2" ]
then
  NUM_MPI_TASKS=2
  NUM_THREADS=2
  SPRAY_BIN=spray_ooc
  PARTITION=image

elif [ $APP == "4" ]
then
  SPRAY_BIN=baseline_ooc
  NUM_MPI_TASKS=2
  NUM_THREADS=2
  PARTITION=image

else
  NUM_MPI_TASKS=1
  NUM_THREADS=1
  COMMAND=""
  echo "[error] invalid input"
fi

CACHE_SIZE=-1
CAMERA="90.172180 84.141418 82.480225 30.000000 28.649426 30.000000"
NUM_PIXEL_SAMPLES=1
NUM_AO_SAMPLES=1
NUM_BOUNCES=1
BLINN_SPECULAR_SHININESS="0.4 0.4 0.4 10"

COMMAND="$MPI_BIN $NUM_MPI_TASKS $SPRAY_BIN_PATH/$SPRAY_BIN \
         --ply-path $PLY_PATH \
         --nthreads $NUM_THREADS \
         -w 512 -h 512 \
         --frames $NUM_FRAMES \
         --mode $MODE \
         --cache-size $CACHE_SIZE \
         --partition $PARTITION \
         --camera $CAMERA \
         --pixel-samples $NUM_PIXEL_SAMPLES \
         --ao-samples $NUM_AO_SAMPLES \
         --bounces $NUM_BOUNCES \
         --blinn $BLINN_SPECULAR_SHININESS \
         $MODEL"

echo "NUM_MPI_TASKS=$NUM_MPI_TASKS"
echo "OMP_NUM_THREADS=$NUM_THREADS"
echo $COMMAND

export OMP_NUM_THREADS=$NUM_THREADS
$COMMAND

