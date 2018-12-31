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
echo "2. simple diffuse, antialiasing, glfw"
echo "3. simple diffuse, antialiasing, film"

read APP

# cmd="$SPRAY_BIN_PATH/spray_insitu_shapes --nthreads 1 -w 512 -h 512 --frames -1 --mode glfw --cache-size -1 --partition insitu --camera 0 0 5 0 0 0 --pixel-samples 1 --ao-samples 1 --ao-mode --bounces 1 --blinn 0.4 0.4 0.4 10 --num-partitions 1 $SPRAY_HOME_PATH/examples/spheres/spheres.spray"

# common settings
NUM_THREADS=1

# common settings
CACHE_SIZE=-1
PARTITION=insitu
# CAMERA="0 0 5 0 0 0"
#CAMERA="-5.069664 0.288000 6.636383 -2.660896 0.288000 2.254856"
CAMERA="-3.274734 1.269994 7.118892 -2.660896 0.288000 2.254856"
# CAMERA="13 2 3 0 0 0"
FOV=90
BLINN_PHONG="0.4 0.4 0.4 10"
NUM_PARTITIONS=1

if [ $APP == "1" ] # simple diffuse, glfw
then
  MODE=glfw
  WIDTH=800
  HEIGHT=600
  NUM_FRAMES=-1
  NUM_PIXEL_SAMPLES=1
  NUM_AO_SAMPLES=1
  NUM_BOUNCES=1

elif [ $APP == "2" ] # simple diffuse, antialiasing, glfw
then
  MODE=glfw
  WIDTH=400
  HEIGHT=300
  NUM_FRAMES=-1
  NUM_PIXEL_SAMPLES=4
  NUM_AO_SAMPLES=2
  NUM_BOUNCES=4

elif [ $APP == "3" ] # simple diffuse, antialiasing, film
then
  MODE=film
  WIDTH=800
  HEIGHT=600
  NUM_FRAMES=1
  NUM_PIXEL_SAMPLES=16
  NUM_AO_SAMPLES=4
  NUM_BOUNCES=4

else # undefined
  echo "[error] invalid input"
  return
fi

cmd="$SHAPES_APP --nthreads $NUM_THREADS \
                 -w $WIDTH -h $HEIGHT \
                 --frames $NUM_FRAMES \
                 --mode $MODE \
                 --cache-size $CACHE_SIZE \
                 --partition $PARTITION \
                 --camera $CAMERA \
                 --pixel-samples $NUM_PIXEL_SAMPLES \
                 --ao-samples $NUM_AO_SAMPLES \
                 --bounces $NUM_BOUNCES \
                 --blinn $BLINN_PHONG \
                 --num-partitions $NUM_PARTITIONS \
                 --fov $FOV \
                 $SCENE"

echo $cmd
$cmd
