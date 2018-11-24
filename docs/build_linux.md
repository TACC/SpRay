# Building SpRay on Linux

Build and install [Embree][2]. Make sure your Embree installation has the `lib` directory. We tested SpRay with Embree v2.17.1.

Check out SpRay with all submodules.

```bash
git clone https://github.com/TACC/SpRay.git
cd SpRay
git submodule init
git submodule update
```
Or simply,

```bash
git clone --recurse-submodules https://github.com/TACC/SpRay.git
cd SpRay
```

Create a directory and build.

```bash
mkdir build
cd build
cmake -DEMBREE_INSTALL_DIR=<path_to_embree_install> ..
make
```

If you wish to install SpRay,

```bash
make install
```

