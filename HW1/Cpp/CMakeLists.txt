cmake_minimum_required(VERSION 2.3)
project(foo CXX)

############################################################
#
#	Build Options
#
############################################################

set(CMAKE_CXX_STANDARD 11)
add_executable(${PROJECT_NAME} main.cpp)

############################################################
#
#	Auxiliary Libraries
#
############################################################
add_library(RigidBody RigidBody.cpp)
target_link_libraries(${PROJECT_NAME} RigidBody)

############################################################
#
#	Libraries
#
############################################################

### Eigen:

find_package (Eigen3 3.3 NO_MODULE)
if (TARGET Eigen3::Eigen)
    # Use the imported target
    include_directories(include)
    include_directories(${CMAKE_SOURCE_DIR})
endif (TARGET Eigen3::Eigen)


## Build-in Libraries
find_package(OCMakeFiles/foo.dir/allpenGL)
find_package(GCMakeFiles/foo.dir/allLUT)
CMakeFiles/foo.dir/all
include_directCMakeFiles/foo.dir/allories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS})

target_link_libraries(
    ${PROJECT_NAME}
    ${OPENGL_gl_LIBRARY}
    ${OPENGL_glu_LIBRARY}
    ${GLUT_LIBRARIES})