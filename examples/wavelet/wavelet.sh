#!/bin/bash

if [ -z $SPRAY_HOME_PATH ]
then
  echo "[error] SPRAY_HOME_PATH not found. Do export SPRAY_HOME_PATH=<path_to_spray_home>."
  return
fi

MODE=$1

if [ "$MODE" != "film" ] && [ "$MODE" != "glfw" ]
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

NUM_PARTITIONS=1

SPRAY_BIN_PATH=$SPRAY_HOME_PATH/build
EXAMPLE_PATH=$SPRAY_HOME_PATH/examples
WAVELET_PATH=$SPRAY_HOME_PATH/examples/wavelet
SCENE=$WAVELET_PATH/wavelet.spray
PLY_PATH=$WAVELET_PATH
MPI_BIN="mpirun -n"

echo "Choose an application (1-4):"
echo "1. spray_insitu_singlethread"
echo "2. spray_insitu_multithread"
echo "3. spray_ooc"
echo "4. baseline_insitu"
echo "5. baseline_ooc"
echo "6. visualize domain bounds based on domain IDs"
echo "7. visualize domain bounds based on partition IDs"

read APP

if [ $APP == "1" ] # spray_insitu_singlethread
then
  NUM_MPI_TASKS=2
  NUM_THREADS=1 # threading not supported
  SPRAY_BIN=spray_insitu_singlethread
  PARTITION=insitu

elif [ $APP == "2" ] # spray_insitu_multithread
then
  NUM_MPI_TASKS=2
  NUM_THREADS=2
  SPRAY_BIN=spray_insitu_multithread
  PARTITION=insitu

elif [ $APP == "3" ] # spray_ooc
then
  NUM_MPI_TASKS=2
  NUM_THREADS=2
  SPRAY_BIN=spray_ooc
  PARTITION=image

elif [ $APP == "4" ] # baseline_insitu
then
  NUM_MPI_TASKS=2
  NUM_THREADS=2
  SPRAY_BIN=baseline_insitu
  PARTITION=insitu

elif [ $APP == "5" ] # baseline_ooc
then
  SPRAY_BIN=baseline_ooc
  NUM_MPI_TASKS=2
  NUM_THREADS=2
  PARTITION=image

elif [ $APP == "6" ] # visualize domain bounds based on domain IDs
then
  SPRAY_BIN=spray_insitu_singlethread # any will work
  NUM_MPI_TASKS=1
  NUM_THREADS=1
  PARTITION=insitu
  MODE=domain

elif [ $APP == "7" ] # visualize domain bounds based on partition IDs
then
  SPRAY_BIN=spray_insitu_singlethread # any will work
  NUM_MPI_TASKS=1
  NUM_THREADS=1
  PARTITION=insitu
  MODE=partition
  NUM_PARTITIONS=2

else # undefined
  NUM_MPI_TASKS=1
  NUM_THREADS=1
  COMMAND=""
  echo "[error] invalid input"
  return
fi

CACHE_SIZE=-1
CAMERA="-5 10 15 0 0 0"
BLINN_SPECULAR_SHININESS="0.4 0.4 0.4 10"

echo "Choose shader type (1 or 2):"
echo "1. ambient occlusion"
echo "2. path tracing"

read SHADER

AO_MODE=""
if [ $SHADER == "1" ] # ambient occlusion
then
  NUM_PIXEL_SAMPLES=1
  NUM_AO_SAMPLES=4
  NUM_BOUNCES=1
  SHADER_TYPE="ambient occlusion"
  AO_MODE="--ao-mode"

elif [ $SHADER == "2" ] # path tracing
then
  NUM_PIXEL_SAMPLES=2
  NUM_AO_SAMPLES=2 # number of samples of diffuse light sources
  NUM_BOUNCES=2
  SHADER_TYPE="path tracing"

else # undefined
  echo "[error] invalid shader input"
  return
fi

echo shader type: $SHADER_TYPE
echo number of pixel samples: $NUM_PIXEL_SAMPLES
echo number of bounces: $NUM_PIXEL_SAMPLES
echo number of light samples: $NUM_AO_SAMPLES

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
         $AO_MODE \
         --bounces $NUM_BOUNCES \
         --blinn $BLINN_SPECULAR_SHININESS \
         --num-partitions $NUM_PARTITIONS \
         $SCENE"

echo "NUM_MPI_TASKS=$NUM_MPI_TASKS"
echo "OMP_NUM_THREADS=$NUM_THREADS"
echo $COMMAND

export OMP_NUM_THREADS=$NUM_THREADS
$COMMAND

