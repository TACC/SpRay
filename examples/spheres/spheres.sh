#!/bin/bash

if [ -z $SPRAY_HOME_PATH ]
then
  echo "[error] SPRAY_HOME_PATH not found. Do export SPRAY_HOME_PATH=<path_to_spray_home>."
  return
fi

SPRAY_BIN_PATH=$SPRAY_HOME_PATH/build
# SHAPES_APP_MULTI_THREAD_HYBRIDGEOM=$SPRAY_BIN_PATH/spray_insitu_multithread_hybridgeometry
APP_SPRAY_INSITU_SINGLE_THREAD=$SPRAY_BIN_PATH/spray_insitu_singlethread
APP_SPRAY_INSITU_MULTI_THREAD=$SPRAY_BIN_PATH/spray_insitu_multithread
APP_SPRAY_OOC_MULTI_THREAD=$SPRAY_BIN_PATH/spray_ooc
EXAMPLE_PATH=$SPRAY_HOME_PATH/examples
SPHERES_PATH=$SPRAY_HOME_PATH/examples/spheres
SPHERES_DIFFUSE_LIGHT=$SPHERES_PATH/spheres-diffuse-light.spray
SPHERES_POINT_LIGHT=$SPHERES_PATH/spheres-point-light.spray
ONE_SPHERE_DIFFUSE_LIGHT=$SPHERES_PATH/sphere.spray
SPHERES_ONE_WAVELET_DIFFUSE_LIGHT=$SPHERES_PATH/sphere-wavelet.spray
SPHERES_ONE_WAVELET_DIFFUSE_LIGHT_SEPARATE_DOMAINS=$SPHERES_PATH/sphere-wavelet-domains.spray
MPI_BIN="mpirun -n"

# common settings
CACHE_SIZE=-1
PARTITION=insitu
# CAMERA="0 0 5 0 0 0"
CAMERA="-4.003263 1.148658 9.021785 -2.584439 1.023140 4.228964"
# CAMERA="13 2 3 0 0 0"
FOV=65
NUM_PARTITIONS=1
PLY_PATH=$EXAMPLE_PATH

# scene

echo "Choose a scene (1-4):"
echo "1. many spheres with a point light" 
echo "2. many spheres with a diffuse light"
echo "3. many spheres and one wavelet with a diffuse light"
echo "4. Scene 3 but objects in separate domains"
echo "5. one sphere with a diffuse light"

read SPHERE_SCENE

if [ $SPHERE_SCENE == "1" ] 
then
  SCENE=$SPHERES_POINT_LIGHT

elif [ $SPHERE_SCENE == "2" ]
then
  SCENE=$SPHERES_DIFFUSE_LIGHT

elif [ $SPHERE_SCENE == "3" ]
then
  SCENE=$SPHERES_ONE_WAVELET_DIFFUSE_LIGHT

  # CAMERA="-9.647732 11.591160 42.547039 -8.228907 11.465641 37.754223"
  # CAMERA="15 15 55 5 5 0"
  CAMERA="23 22 48 5 5 0"

elif [ $SPHERE_SCENE == "4" ]
then
  SCENE=$SPHERES_ONE_WAVELET_DIFFUSE_LIGHT_SEPARATE_DOMAINS

  # CAMERA="-9.647732 11.591160 42.547039 -8.228907 11.465641 37.754223"
  # CAMERA="15 15 55 5 5 0"
  CAMERA="23 22 48 5 5 0"

elif [ $SPHERE_SCENE == "5" ]
then
  SCENE=$ONE_SPHERE_DIFFUSE_LIGHT

else
  echo "[error] invalid input"
  return
fi

# app

echo "Choose an app  (1-2):"
echo "1. insitu, single-thread" 
echo "2. insitu, multi-thread"
echo "3. ooc, multi-thread"

read APP

if [ $APP == "1" ] 
then
  APP_BIN=$APP_SPRAY_INSITU_SINGLE_THREAD
  NUM_THREADS=1

elif [ $APP == "2" ] 
then
  APP_BIN=$APP_SPRAY_INSITU_MULTI_THREAD
  NUM_THREADS=4

elif [ $APP == "3" ] 
then
  APP_BIN=$APP_SPRAY_OOC_MULTI_THREAD
  NUM_THREADS=4
  PARTITION=image

else
  echo "[error] invalid input"
  return
fi

# setup

echo "Choose a setup (1-3):"
echo "1. glfw, no antialiasing" 
echo "2. glfw, antialiasing"
echo "3. film, antialiasing"

read SETUP

export OMP_NUM_THREADS=$NUM_THREADS

if [ $SETUP == "1" ] # glfw
then
  MODE=glfw
  WIDTH=800
  HEIGHT=600
  NUM_FRAMES=-1
  NUM_PIXEL_SAMPLES=1
  NUM_AO_SAMPLES=1
  NUM_BOUNCES=4

elif [ $SETUP == "2" ] # antialiasing, glfw
then
  MODE=glfw
  WIDTH=400
  HEIGHT=300
  NUM_FRAMES=-1
  NUM_PIXEL_SAMPLES=4
  NUM_AO_SAMPLES=2
  NUM_BOUNCES=4

elif [ $SETUP == "3" ] # antialiasing, film
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

cmd="$APP_BIN --ply-path $PLY_PATH \
              --nthreads $NUM_THREADS \
              -w $WIDTH -h $HEIGHT \
              --frames $NUM_FRAMES \
              --mode $MODE \
              --cache-size $CACHE_SIZE \
              --partition $PARTITION \
              --camera $CAMERA \
              --pixel-samples $NUM_PIXEL_SAMPLES \
              --ao-samples $NUM_AO_SAMPLES \
              --bounces $NUM_BOUNCES \
              --num-partitions $NUM_PARTITIONS \
              --fov $FOV \
              $SCENE"

echo $cmd
$cmd
