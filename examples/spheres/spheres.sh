# cmd="$SPRAY_BIN_PATH/spray_insitu_shapes --nthreads 1 -w 512 -h 512 --frames -1 --mode glfw --cache-size -1 --partition insitu --camera 0 0 5 0 0 0 --pixel-samples 1 --ao-samples 1 --ao-mode --bounces 1 --blinn 0.4 0.4 0.4 10 --num-partitions 1 $SPRAY_HOME_PATH/examples/spheres/spheres.spray"
cmd="$SPRAY_BIN_PATH/spray_insitu_shapes --nthreads 1 -w 512 -h 512 --frames -1 --mode glfw --cache-size -1 --partition insitu --camera 0 0 5 0 0 0 --pixel-samples 1 --ao-samples 1 --bounces 1 --blinn 0.4 0.4 0.4 10 --num-partitions 1 $SPRAY_HOME_PATH/examples/spheres/spheres.spray"
echo $cmd
$cmd
