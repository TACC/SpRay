#!/bin/bash

if [ -z $SPRAY_HOME_PATH ]
then
  echo "[error] SPRAY_HOME_PATH not found. Do export SPRAY_HOME_PATH=<path_to_spray_home>."
  return
fi

SPRAY_BIN_PATH=$SPRAY_HOME_PATH/build
SHAPES_APP_SINGLE_THREAD=$SPRAY_BIN_PATH/spray_insitu_singlethread_shapes
SHAPES_APP_MULTI_THREAD=$SPRAY_BIN_PATH/spray_insitu_multithread_shapes
EXAMPLE_PATH=$SPRAY_HOME_PATH/examples
SPHERES_PATH=$SPRAY_HOME_PATH/examples/spheres
SPHERES_DIFFUSE_LIGHT=$SPHERES_PATH/spheres-diffuse-light.spray
SPHERES_POINT_LIGHT=$SPHERES_PATH/spheres-point-light.spray
ONE_SPHERE_DIFFUSE_LIGHT=$SPHERES_PATH/sphere.spray
MPI_BIN="mpirun -n"

echo "Choose a setup (1-2):"
echo "1. simple diffuse, glfw" 
echo "2. simple diffuse, antialiasing, glfw"
echo "3. simple diffuse, antialiasing, film"

read APP

echo "Choose a scene (1,2,3):"
echo "1. spheres with a point light" 
echo "2. spheres with a diffuse light"
echo "3. a sphere with a diffuse light"

read LIGHT

if [ $LIGHT == "1" ] 
then
  SCENE=$SCENE_POINT_LIGHT

elif [ $LIGHT == "2" ]
then
  SCENE=$SPHERES_DIFFUSE_LIGHT

elif [ $LIGHT == "3" ]
then
  SCENE=$ONE_SPHERE_DIFFUSE_LIGHT

else
  echo "[error] invalid input"
  return
fi

echo "Choose a threading mode (1-2):"
echo "1. single-thread" 
echo "2. multi-thread"

read THREADING

if [ $THREADING == "1" ] 
then
  SHAPES_APP=$SHAPES_APP_SINGLE_THREAD
  NUM_THREADS=1

elif [ $THREADING == "2" ]
then
  SHAPES_APP=$SHAPES_APP_MULTI_THREAD
  NUM_THREADS=4

else
  echo "[error] invalid input"
  return
fi

export OMP_NUM_THREADS=$NUM_THREADS

# common settings
CACHE_SIZE=-1
PARTITION=insitu
# CAMERA="0 0 5 0 0 0"
CAMERA="-4.003263 1.148658 9.021785 -2.584439 1.023140 4.228964"
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
  NUM_BOUNCES=4

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
