# Example 1: rendering a scene of 64 wavelet domains

As you run the bash scripts below, you'll be asked to choose an application from a list of applications. Those starting with `base` are an implementation of the baseline algorithm; those starting with `spray` are an implementation using the speculative technique.

Notice that the scripts launch MPI tasks using the `mpirun` command and specify the number of OpenMP threads by setting the environment variable `OMP_NUM_THREADS`. You may have to modify such settings based on your system requirements.

Additionally, if you wish to use installed binaries, you'll have to modify the variable `SPRAY_BIN_PATH` in the scripts and set runtime search paths, `LD_LIBRARY_PATH` or `DYLD_LIBRARY_PATH`, as needed.

## Running the script

Edit the script file to set environment variables: `examples/env_spray_linux.sh` for Linux or `examples/env_spray_macos.sh` for macOS.

For Linux:
```bash
source examples/env_spray_linux.sh
```

For macOS:
```bash
source examples/env_spray_macos.sh
```

Depending on the mode you choose, run the script as shown below. N and n are optional, where N is the number of MPI tasks and n is the number of OpenMP threads for each MPI task. If these parameters are not set, 2 MPI tasks and 2 OpenMP threads per each task will be launched. For some special cases, the script discards N and n, using interally hardcoded values. For non-rendering modes, 1 MPI task and 1 thread are launched, and for the single-thread rendering mode, only 1 thread is launched with N MPI tasks. For the rendering settings used in this example, refer to the `wavelets64.sh` script.

For film mode:

```bash
source $SPRAY_HOME_PATH/examples/wavelets64/wavelets64.sh film [N] [n]
Type a number from the application list.
display spray.ppm
```

For glfw mode:

```bash
source $SPRAY_HOME_PATH/examples/wavelets64/wavelets64.sh glfw [N] [n]
Type a number from the application list.
Type the q-key to close the window.
```

You should see the following as a result:

![wavelets.jpg](assets/img/wavelets64.jpg)
