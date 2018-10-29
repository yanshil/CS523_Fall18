# HW 2

## Submision Files

* Source Code

```
# find . -maxdepth 2
.
./README.md
./Problem1
./Problem1/CMakeLists.txt
./Problem1/hw2.cpp
./Problem1/.vscode
./Problem2
./Problem2/RigidBody.h
./Problem2/Contact.h
./Problem2/include
./Problem2/RigidBody.cpp
./Problem2/Makefile
./Problem2/CMakeLists.txt
./Problem2/.vscode
./Problem2/Contact.cpp
./Problem2/main.cpp
```

* Output

```
./Output/Problem1_circle.out
./Output/Problem1_ellipse.out
./Output/Problem1_square.out
./Output/Problem1_shpere.out
./Output/Problem1_ellipsolid.out
./Output/Problem1_cubic.out

./Output/Problem2_vedio1.mp4
./Output/Problem2_vedio2.mp4
```

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
* Grid size can be set by parsing arguments (when d=2): `./bin/foo -size 500 500` (default is 500, 500)

```
# Example output file generating
# (When d = 2)
./bin/foo >> ~/Workspace/CS523/HW2/Output/Problem1_circle.out
1
0 0

./bin/foo >> ~/Workspace/CS523/HW2/Output/Problem1_ellipse.out
2
-10 -10

./bin/foo >> ~/Workspace/CS523/HW2/Output/Problem1_square.out
3
-5.62 1.2424

# Change (d=3) before building the project for 3D case.
# Example output file generated
./bin/foo >> ~/Workspace/CS523/HW2/Output/Problem1_shpere.out
1
0 0 0

./bin/foo >> ~/Workspace/CS523/HW2/Output/Problem1_ellipsolid.out
2
0.001 0.245 0.9983

./bin/foo >> ~/Workspace/CS523/HW2/Output/Problem1_cubic.out
3
5 5 5
```

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

