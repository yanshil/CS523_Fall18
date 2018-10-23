cmake_minimum_required(VERSION 2.6)
project(foo CXX)

############################################################
#
#	Build Options
#
############################################################

set (CMAKE_CXX_STANDARD 11)
add_executable(${PROJECT_NAME} test.cpp)

### Eigen:

#find_package (Eigen3 3.3 NO_MODULE) 
#if (TARGET Eigen3::Eigen)
    # Use the imported target
    include_directories(include)
#    include_directories(set (CMAKE_CXX_STANDARD 11)${CMAKE_SOURCE_DIR})
# endif (TARGET Eigen3::Eigen)

############################################################
#
#	Libraries
#
############################################################

# Target for Tools
add_library( _nova_tools
  src/Grids/Grid.cpp
  src/Krylov_Solvers/Conjugate_Gradient.cpp
  src/Krylov_Solvers/Krylov_System_Base.cpp
  src/Log/Debug_Utilities.cpp
  src/Log/Log.cpp
  src/Log/Log_Entry.cpp
  src/Parsing/Parse_Args.cpp
  src/Random_Numbers/MT19937.cpp
  src/Utilities/Driver.cpp
  src/Utilities/Example.cpp
  src/Utilities/File_Utilities.cpp
  src/Utilities/Pthread_Queue.cpp
  src/Utilities/Timer.cpp
  src/Vectors/Vector_2D.cpp
  src/Vectors/Vector_3D.cpp
  )

find_package(Boost 1.54 REQUIRED COMPONENTS filesystem)

if(USE_DOUBLES)
  target_compile_definitions(_nova_tools PUBLIC COMPILE_WITH_DOUBLE_SUPPORT=1)
endif()

target_include_directories( _nova_tools PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src )

target_compile_features( _nova_tools
  PUBLIC
  cxx_variadic_templates
  cxx_nullptr
  )

target_link_libraries( _nova_tools
  PUBLIC
  ${Boost_LIBRARIES}
  )

add_library( nova::tools ALIAS _nova_tools )
