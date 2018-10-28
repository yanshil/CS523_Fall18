# HW 2

## Problem 1

* Project is build with [Nova](https://github.com/OrionQuest/Nova)
```
cd Nova/Projects
cp -r HW2/Problem1/ .
```

* Config with `ccmake`
```
cd ../
mkdir build
cd build
ccmake ..
```
Set `CMAKE_BUILD_TYPE` to `Release` and turn on the following flags: `ENABLE_GRID_PLUGIN`, `USE_DOUBLE`. Click `c` to config and `g` to generate Makefile.

* Make with `make`

* Run with `./bin/foo`

### Usage
* Dimension (d=2 or d=3) is modified by modifying Line20 in `hw2.cpp`
* Grid size can be set by parsing arguments (when d=2): `./bin/foo -size 500 500`

## Problem 2

* Project dependency: 
	* OpenGL
	* Glut
	* `Eigen`( No need to install): Files has been included in the project folder.
* `Makefile` has been included. 

### Running Guide
```
cd HW2/Problem2/
make
./foo
```

* Initial velocity and omega can be set by changing values in `rb.setVelocity()` and `rb.setOmega()` in `main.cpp`

