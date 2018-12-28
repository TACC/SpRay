#!/bin/bash

if [ -z $SPRAY_HOME_PATH ]
then
  echo "[error] SPRAY_HOME_PATH not found. Do export SPRAY_HOME_PATH=<path_to_spray_home>."
  return
fi

SPRAY_BIN_PATH=$SPRAY_HOME_PATH/build
SHAPES_APP=$SPRAY_BIN_PATH/spray_insitu_shapes
EXAMPLE_PATH=$SPRAY_HOME_PATH/examples
SPHERES_PATH=$SPRAY_HOME_PATH/examples/spheres
SCENE=$SPHERES_PATH/spheres.spray
MPI_BIN="mpirun -n"

echo "Choose a setup (1-2):"
echo "1. simple diffuse, glfw" 
echo "2. simple diffuse, antialiasing, film"

read APP

# cmd="$SPRAY_BIN_PATH/spray_insitu_shapes --nthreads 1 -w 512 -h 512 --frames -1 --mode glfw --cache-size -1 --partition insitu --camera 0 0 5 0 0 0 --pixel-samples 1 --ao-samples 1 --ao-mode --bounces 1 --blinn 0.4 0.4 0.4 10 --num-partitions 1 $SPRAY_HOME_PATH/examples/spheres/spheres.spray"

# common settings
WIDTH=400
HEIGHT=300
NUM_THREADS=1

if [ $APP == "1" ] # simple diffuse, glfw
then
  MODE=glfw
  NUM_FRAMES=-1
  NUM_PIXEL_SAMPLES=1
  cmd="$SHAPES_APP --nthreads $NUM_THREADS -w $WIDTH -h $HEIGHT --frames $NUM_FRAMES --mode $MODE --cache-size -1 --partition insitu --camera 0 0 5 0 0 0 --pixel-samples $NUM_PIXEL_SAMPLES --ao-samples 1 --bounces 1 --blinn 0.4 0.4 0.4 10 --num-partitions 1 $SCENE"
  
elif [ $APP == "2" ] # simple diffuse, antialiasing
then
  # simple diffuse, antialiasing, film
  MODE=film
  NUM_FRAMES=1
  NUM_PIXEL_SAMPLES=16
  cmd="$SHAPES_APP --nthreads $NUM_THREADS -w $WIDTH -h $HEIGHT --frames $NUM_FRAMES --mode $MODE --cache-size -1 --partition insitu --camera 0 0 5 0 0 0 --pixel-samples $NUM_PIXEL_SAMPLES --ao-samples 1 --bounces 1 --blinn 0.4 0.4 0.4 10 --num-partitions 1 $SCENE"

else # undefined
  echo "[error] invalid input"
  return
fi

echo $cmd
$cmd
