# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/yanshi/Workspace/CS523/HW2/Problem1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/yanshi/Workspace/CS523/HW2/Problem1/build

# Include any dependencies generated for this target.
include CMakeFiles/_nova_tools.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/_nova_tools.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/_nova_tools.dir/flags.make

CMakeFiles/_nova_tools.dir/src/Grids/Grid.cpp.o: CMakeFiles/_nova_tools.dir/flags.make
CMakeFiles/_nova_tools.dir/src/Grids/Grid.cpp.o: ../src/Grids/Grid.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yanshi/Workspace/CS523/HW2/Problem1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/_nova_tools.dir/src/Grids/Grid.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/_nova_tools.dir/src/Grids/Grid.cpp.o -c /home/yanshi/Workspace/CS523/HW2/Problem1/src/Grids/Grid.cpp

CMakeFiles/_nova_tools.dir/src/Grids/Grid.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/_nova_tools.dir/src/Grids/Grid.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/yanshi/Workspace/CS523/HW2/Problem1/src/Grids/Grid.cpp > CMakeFiles/_nova_tools.dir/src/Grids/Grid.cpp.i

CMakeFiles/_nova_tools.dir/src/Grids/Grid.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/_nova_tools.dir/src/Grids/Grid.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/yanshi/Workspace/CS523/HW2/Problem1/src/Grids/Grid.cpp -o CMakeFiles/_nova_tools.dir/src/Grids/Grid.cpp.s

CMakeFiles/_nova_tools.dir/src/Krylov_Solvers/Conjugate_Gradient.cpp.o: CMakeFiles/_nova_tools.dir/flags.make
CMakeFiles/_nova_tools.dir/src/Krylov_Solvers/Conjugate_Gradient.cpp.o: ../src/Krylov_Solvers/Conjugate_Gradient.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yanshi/Workspace/CS523/HW2/Problem1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/_nova_tools.dir/src/Krylov_Solvers/Conjugate_Gradient.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/_nova_tools.dir/src/Krylov_Solvers/Conjugate_Gradient.cpp.o -c /home/yanshi/Workspace/CS523/HW2/Problem1/src/Krylov_Solvers/Conjugate_Gradient.cpp

CMakeFiles/_nova_tools.dir/src/Krylov_Solvers/Conjugate_Gradient.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/_nova_tools.dir/src/Krylov_Solvers/Conjugate_Gradient.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/yanshi/Workspace/CS523/HW2/Problem1/src/Krylov_Solvers/Conjugate_Gradient.cpp > CMakeFiles/_nova_tools.dir/src/Krylov_Solvers/Conjugate_Gradient.cpp.i

CMakeFiles/_nova_tools.dir/src/Krylov_Solvers/Conjugate_Gradient.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/_nova_tools.dir/src/Krylov_Solvers/Conjugate_Gradient.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/yanshi/Workspace/CS523/HW2/Problem1/src/Krylov_Solvers/Conjugate_Gradient.cpp -o CMakeFiles/_nova_tools.dir/src/Krylov_Solvers/Conjugate_Gradient.cpp.s

CMakeFiles/_nova_tools.dir/src/Krylov_Solvers/Krylov_System_Base.cpp.o: CMakeFiles/_nova_tools.dir/flags.make
CMakeFiles/_nova_tools.dir/src/Krylov_Solvers/Krylov_System_Base.cpp.o: ../src/Krylov_Solvers/Krylov_System_Base.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yanshi/Workspace/CS523/HW2/Problem1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/_nova_tools.dir/src/Krylov_Solvers/Krylov_System_Base.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/_nova_tools.dir/src/Krylov_Solvers/Krylov_System_Base.cpp.o -c /home/yanshi/Workspace/CS523/HW2/Problem1/src/Krylov_Solvers/Krylov_System_Base.cpp

CMakeFiles/_nova_tools.dir/src/Krylov_Solvers/Krylov_System_Base.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/_nova_tools.dir/src/Krylov_Solvers/Krylov_System_Base.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/yanshi/Workspace/CS523/HW2/Problem1/src/Krylov_Solvers/Krylov_System_Base.cpp > CMakeFiles/_nova_tools.dir/src/Krylov_Solvers/Krylov_System_Base.cpp.i

CMakeFiles/_nova_tools.dir/src/Krylov_Solvers/Krylov_System_Base.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/_nova_tools.dir/src/Krylov_Solvers/Krylov_System_Base.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/yanshi/Workspace/CS523/HW2/Problem1/src/Krylov_Solvers/Krylov_System_Base.cpp -o CMakeFiles/_nova_tools.dir/src/Krylov_Solvers/Krylov_System_Base.cpp.s

CMakeFiles/_nova_tools.dir/src/Log/Debug_Utilities.cpp.o: CMakeFiles/_nova_tools.dir/flags.make
CMakeFiles/_nova_tools.dir/src/Log/Debug_Utilities.cpp.o: ../src/Log/Debug_Utilities.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yanshi/Workspace/CS523/HW2/Problem1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/_nova_tools.dir/src/Log/Debug_Utilities.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/_nova_tools.dir/src/Log/Debug_Utilities.cpp.o -c /home/yanshi/Workspace/CS523/HW2/Problem1/src/Log/Debug_Utilities.cpp

CMakeFiles/_nova_tools.dir/src/Log/Debug_Utilities.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/_nova_tools.dir/src/Log/Debug_Utilities.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/yanshi/Workspace/CS523/HW2/Problem1/src/Log/Debug_Utilities.cpp > CMakeFiles/_nova_tools.dir/src/Log/Debug_Utilities.cpp.i

CMakeFiles/_nova_tools.dir/src/Log/Debug_Utilities.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/_nova_tools.dir/src/Log/Debug_Utilities.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/yanshi/Workspace/CS523/HW2/Problem1/src/Log/Debug_Utilities.cpp -o CMakeFiles/_nova_tools.dir/src/Log/Debug_Utilities.cpp.s

CMakeFiles/_nova_tools.dir/src/Log/Log.cpp.o: CMakeFiles/_nova_tools.dir/flags.make
CMakeFiles/_nova_tools.dir/src/Log/Log.cpp.o: ../src/Log/Log.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yanshi/Workspace/CS523/HW2/Problem1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/_nova_tools.dir/src/Log/Log.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/_nova_tools.dir/src/Log/Log.cpp.o -c /home/yanshi/Workspace/CS523/HW2/Problem1/src/Log/Log.cpp

CMakeFiles/_nova_tools.dir/src/Log/Log.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/_nova_tools.dir/src/Log/Log.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/yanshi/Workspace/CS523/HW2/Problem1/src/Log/Log.cpp > CMakeFiles/_nova_tools.dir/src/Log/Log.cpp.i

CMakeFiles/_nova_tools.dir/src/Log/Log.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/_nova_tools.dir/src/Log/Log.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/yanshi/Workspace/CS523/HW2/Problem1/src/Log/Log.cpp -o CMakeFiles/_nova_tools.dir/src/Log/Log.cpp.s

CMakeFiles/_nova_tools.dir/src/Log/Log_Entry.cpp.o: CMakeFiles/_nova_tools.dir/flags.make
CMakeFiles/_nova_tools.dir/src/Log/Log_Entry.cpp.o: ../src/Log/Log_Entry.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yanshi/Workspace/CS523/HW2/Problem1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/_nova_tools.dir/src/Log/Log_Entry.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/_nova_tools.dir/src/Log/Log_Entry.cpp.o -c /home/yanshi/Workspace/CS523/HW2/Problem1/src/Log/Log_Entry.cpp

CMakeFiles/_nova_tools.dir/src/Log/Log_Entry.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/_nova_tools.dir/src/Log/Log_Entry.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/yanshi/Workspace/CS523/HW2/Problem1/src/Log/Log_Entry.cpp > CMakeFiles/_nova_tools.dir/src/Log/Log_Entry.cpp.i

CMakeFiles/_nova_tools.dir/src/Log/Log_Entry.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/_nova_tools.dir/src/Log/Log_Entry.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/yanshi/Workspace/CS523/HW2/Problem1/src/Log/Log_Entry.cpp -o CMakeFiles/_nova_tools.dir/src/Log/Log_Entry.cpp.s

CMakeFiles/_nova_tools.dir/src/Parsing/Parse_Args.cpp.o: CMakeFiles/_nova_tools.dir/flags.make
CMakeFiles/_nova_tools.dir/src/Parsing/Parse_Args.cpp.o: ../src/Parsing/Parse_Args.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yanshi/Workspace/CS523/HW2/Problem1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/_nova_tools.dir/src/Parsing/Parse_Args.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/_nova_tools.dir/src/Parsing/Parse_Args.cpp.o -c /home/yanshi/Workspace/CS523/HW2/Problem1/src/Parsing/Parse_Args.cpp

CMakeFiles/_nova_tools.dir/src/Parsing/Parse_Args.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/_nova_tools.dir/src/Parsing/Parse_Args.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/yanshi/Workspace/CS523/HW2/Problem1/src/Parsing/Parse_Args.cpp > CMakeFiles/_nova_tools.dir/src/Parsing/Parse_Args.cpp.i

CMakeFiles/_nova_tools.dir/src/Parsing/Parse_Args.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/_nova_tools.dir/src/Parsing/Parse_Args.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/yanshi/Workspace/CS523/HW2/Problem1/src/Parsing/Parse_Args.cpp -o CMakeFiles/_nova_tools.dir/src/Parsing/Parse_Args.cpp.s

CMakeFiles/_nova_tools.dir/src/Random_Numbers/MT19937.cpp.o: CMakeFiles/_nova_tools.dir/flags.make
CMakeFiles/_nova_tools.dir/src/Random_Numbers/MT19937.cpp.o: ../src/Random_Numbers/MT19937.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yanshi/Workspace/CS523/HW2/Problem1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/_nova_tools.dir/src/Random_Numbers/MT19937.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/_nova_tools.dir/src/Random_Numbers/MT19937.cpp.o -c /home/yanshi/Workspace/CS523/HW2/Problem1/src/Random_Numbers/MT19937.cpp

CMakeFiles/_nova_tools.dir/src/Random_Numbers/MT19937.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/_nova_tools.dir/src/Random_Numbers/MT19937.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/yanshi/Workspace/CS523/HW2/Problem1/src/Random_Numbers/MT19937.cpp > CMakeFiles/_nova_tools.dir/src/Random_Numbers/MT19937.cpp.i

CMakeFiles/_nova_tools.dir/src/Random_Numbers/MT19937.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/_nova_tools.dir/src/Random_Numbers/MT19937.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/yanshi/Workspace/CS523/HW2/Problem1/src/Random_Numbers/MT19937.cpp -o CMakeFiles/_nova_tools.dir/src/Random_Numbers/MT19937.cpp.s

CMakeFiles/_nova_tools.dir/src/Utilities/Driver.cpp.o: CMakeFiles/_nova_tools.dir/flags.make
CMakeFiles/_nova_tools.dir/src/Utilities/Driver.cpp.o: ../src/Utilities/Driver.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yanshi/Workspace/CS523/HW2/Problem1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/_nova_tools.dir/src/Utilities/Driver.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/_nova_tools.dir/src/Utilities/Driver.cpp.o -c /home/yanshi/Workspace/CS523/HW2/Problem1/src/Utilities/Driver.cpp

CMakeFiles/_nova_tools.dir/src/Utilities/Driver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/_nova_tools.dir/src/Utilities/Driver.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/yanshi/Workspace/CS523/HW2/Problem1/src/Utilities/Driver.cpp > CMakeFiles/_nova_tools.dir/src/Utilities/Driver.cpp.i

CMakeFiles/_nova_tools.dir/src/Utilities/Driver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/_nova_tools.dir/src/Utilities/Driver.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/yanshi/Workspace/CS523/HW2/Problem1/src/Utilities/Driver.cpp -o CMakeFiles/_nova_tools.dir/src/Utilities/Driver.cpp.s

CMakeFiles/_nova_tools.dir/src/Utilities/Example.cpp.o: CMakeFiles/_nova_tools.dir/flags.make
CMakeFiles/_nova_tools.dir/src/Utilities/Example.cpp.o: ../src/Utilities/Example.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yanshi/Workspace/CS523/HW2/Problem1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/_nova_tools.dir/src/Utilities/Example.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/_nova_tools.dir/src/Utilities/Example.cpp.o -c /home/yanshi/Workspace/CS523/HW2/Problem1/src/Utilities/Example.cpp

CMakeFiles/_nova_tools.dir/src/Utilities/Example.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/_nova_tools.dir/src/Utilities/Example.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/yanshi/Workspace/CS523/HW2/Problem1/src/Utilities/Example.cpp > CMakeFiles/_nova_tools.dir/src/Utilities/Example.cpp.i

CMakeFiles/_nova_tools.dir/src/Utilities/Example.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/_nova_tools.dir/src/Utilities/Example.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/yanshi/Workspace/CS523/HW2/Problem1/src/Utilities/Example.cpp -o CMakeFiles/_nova_tools.dir/src/Utilities/Example.cpp.s

CMakeFiles/_nova_tools.dir/src/Utilities/File_Utilities.cpp.o: CMakeFiles/_nova_tools.dir/flags.make
CMakeFiles/_nova_tools.dir/src/Utilities/File_Utilities.cpp.o: ../src/Utilities/File_Utilities.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yanshi/Workspace/CS523/HW2/Problem1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/_nova_tools.dir/src/Utilities/File_Utilities.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/_nova_tools.dir/src/Utilities/File_Utilities.cpp.o -c /home/yanshi/Workspace/CS523/HW2/Problem1/src/Utilities/File_Utilities.cpp

CMakeFiles/_nova_tools.dir/src/Utilities/File_Utilities.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/_nova_tools.dir/src/Utilities/File_Utilities.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/yanshi/Workspace/CS523/HW2/Problem1/src/Utilities/File_Utilities.cpp > CMakeFiles/_nova_tools.dir/src/Utilities/File_Utilities.cpp.i

CMakeFiles/_nova_tools.dir/src/Utilities/File_Utilities.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/_nova_tools.dir/src/Utilities/File_Utilities.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/yanshi/Workspace/CS523/HW2/Problem1/src/Utilities/File_Utilities.cpp -o CMakeFiles/_nova_tools.dir/src/Utilities/File_Utilities.cpp.s

CMakeFiles/_nova_tools.dir/src/Utilities/Pthread_Queue.cpp.o: CMakeFiles/_nova_tools.dir/flags.make
CMakeFiles/_nova_tools.dir/src/Utilities/Pthread_Queue.cpp.o: ../src/Utilities/Pthread_Queue.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yanshi/Workspace/CS523/HW2/Problem1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/_nova_tools.dir/src/Utilities/Pthread_Queue.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/_nova_tools.dir/src/Utilities/Pthread_Queue.cpp.o -c /home/yanshi/Workspace/CS523/HW2/Problem1/src/Utilities/Pthread_Queue.cpp

CMakeFiles/_nova_tools.dir/src/Utilities/Pthread_Queue.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/_nova_tools.dir/src/Utilities/Pthread_Queue.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/yanshi/Workspace/CS523/HW2/Problem1/src/Utilities/Pthread_Queue.cpp > CMakeFiles/_nova_tools.dir/src/Utilities/Pthread_Queue.cpp.i

CMakeFiles/_nova_tools.dir/src/Utilities/Pthread_Queue.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/_nova_tools.dir/src/Utilities/Pthread_Queue.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/yanshi/Workspace/CS523/HW2/Problem1/src/Utilities/Pthread_Queue.cpp -o CMakeFiles/_nova_tools.dir/src/Utilities/Pthread_Queue.cpp.s

CMakeFiles/_nova_tools.dir/src/Utilities/Timer.cpp.o: CMakeFiles/_nova_tools.dir/flags.make
CMakeFiles/_nova_tools.dir/src/Utilities/Timer.cpp.o: ../src/Utilities/Timer.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yanshi/Workspace/CS523/HW2/Problem1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object CMakeFiles/_nova_tools.dir/src/Utilities/Timer.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/_nova_tools.dir/src/Utilities/Timer.cpp.o -c /home/yanshi/Workspace/CS523/HW2/Problem1/src/Utilities/Timer.cpp

CMakeFiles/_nova_tools.dir/src/Utilities/Timer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/_nova_tools.dir/src/Utilities/Timer.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/yanshi/Workspace/CS523/HW2/Problem1/src/Utilities/Timer.cpp > CMakeFiles/_nova_tools.dir/src/Utilities/Timer.cpp.i

CMakeFiles/_nova_tools.dir/src/Utilities/Timer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/_nova_tools.dir/src/Utilities/Timer.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/yanshi/Workspace/CS523/HW2/Problem1/src/Utilities/Timer.cpp -o CMakeFiles/_nova_tools.dir/src/Utilities/Timer.cpp.s

CMakeFiles/_nova_tools.dir/src/Vectors/Vector_2D.cpp.o: CMakeFiles/_nova_tools.dir/flags.make
CMakeFiles/_nova_tools.dir/src/Vectors/Vector_2D.cpp.o: ../src/Vectors/Vector_2D.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yanshi/Workspace/CS523/HW2/Problem1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building CXX object CMakeFiles/_nova_tools.dir/src/Vectors/Vector_2D.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/_nova_tools.dir/src/Vectors/Vector_2D.cpp.o -c /home/yanshi/Workspace/CS523/HW2/Problem1/src/Vectors/Vector_2D.cpp

CMakeFiles/_nova_tools.dir/src/Vectors/Vector_2D.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/_nova_tools.dir/src/Vectors/Vector_2D.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/yanshi/Workspace/CS523/HW2/Problem1/src/Vectors/Vector_2D.cpp > CMakeFiles/_nova_tools.dir/src/Vectors/Vector_2D.cpp.i

CMakeFiles/_nova_tools.dir/src/Vectors/Vector_2D.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/_nova_tools.dir/src/Vectors/Vector_2D.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/yanshi/Workspace/CS523/HW2/Problem1/src/Vectors/Vector_2D.cpp -o CMakeFiles/_nova_tools.dir/src/Vectors/Vector_2D.cpp.s

CMakeFiles/_nova_tools.dir/src/Vectors/Vector_3D.cpp.o: CMakeFiles/_nova_tools.dir/flags.make
CMakeFiles/_nova_tools.dir/src/Vectors/Vector_3D.cpp.o: ../src/Vectors/Vector_3D.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yanshi/Workspace/CS523/HW2/Problem1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Building CXX object CMakeFiles/_nova_tools.dir/src/Vectors/Vector_3D.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/_nova_tools.dir/src/Vectors/Vector_3D.cpp.o -c /home/yanshi/Workspace/CS523/HW2/Problem1/src/Vectors/Vector_3D.cpp

CMakeFiles/_nova_tools.dir/src/Vectors/Vector_3D.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/_nova_tools.dir/src/Vectors/Vector_3D.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/yanshi/Workspace/CS523/HW2/Problem1/src/Vectors/Vector_3D.cpp > CMakeFiles/_nova_tools.dir/src/Vectors/Vector_3D.cpp.i

CMakeFiles/_nova_tools.dir/src/Vectors/Vector_3D.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/_nova_tools.dir/src/Vectors/Vector_3D.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/yanshi/Workspace/CS523/HW2/Problem1/src/Vectors/Vector_3D.cpp -o CMakeFiles/_nova_tools.dir/src/Vectors/Vector_3D.cpp.s

# Object files for target _nova_tools
_nova_tools_OBJECTS = \
"CMakeFiles/_nova_tools.dir/src/Grids/Grid.cpp.o" \
"CMakeFiles/_nova_tools.dir/src/Krylov_Solvers/Conjugate_Gradient.cpp.o" \
"CMakeFiles/_nova_tools.dir/src/Krylov_Solvers/Krylov_System_Base.cpp.o" \
"CMakeFiles/_nova_tools.dir/src/Log/Debug_Utilities.cpp.o" \
"CMakeFiles/_nova_tools.dir/src/Log/Log.cpp.o" \
"CMakeFiles/_nova_tools.dir/src/Log/Log_Entry.cpp.o" \
"CMakeFiles/_nova_tools.dir/src/Parsing/Parse_Args.cpp.o" \
"CMakeFiles/_nova_tools.dir/src/Random_Numbers/MT19937.cpp.o" \
"CMakeFiles/_nova_tools.dir/src/Utilities/Driver.cpp.o" \
"CMakeFiles/_nova_tools.dir/src/Utilities/Example.cpp.o" \
"CMakeFiles/_nova_tools.dir/src/Utilities/File_Utilities.cpp.o" \
"CMakeFiles/_nova_tools.dir/src/Utilities/Pthread_Queue.cpp.o" \
"CMakeFiles/_nova_tools.dir/src/Utilities/Timer.cpp.o" \
"CMakeFiles/_nova_tools.dir/src/Vectors/Vector_2D.cpp.o" \
"CMakeFiles/_nova_tools.dir/src/Vectors/Vector_3D.cpp.o"

# External object files for target _nova_tools
_nova_tools_EXTERNAL_OBJECTS =

lib_nova_tools.a: CMakeFiles/_nova_tools.dir/src/Grids/Grid.cpp.o
lib_nova_tools.a: CMakeFiles/_nova_tools.dir/src/Krylov_Solvers/Conjugate_Gradient.cpp.o
lib_nova_tools.a: CMakeFiles/_nova_tools.dir/src/Krylov_Solvers/Krylov_System_Base.cpp.o
lib_nova_tools.a: CMakeFiles/_nova_tools.dir/src/Log/Debug_Utilities.cpp.o
lib_nova_tools.a: CMakeFiles/_nova_tools.dir/src/Log/Log.cpp.o
lib_nova_tools.a: CMakeFiles/_nova_tools.dir/src/Log/Log_Entry.cpp.o
lib_nova_tools.a: CMakeFiles/_nova_tools.dir/src/Parsing/Parse_Args.cpp.o
lib_nova_tools.a: CMakeFiles/_nova_tools.dir/src/Random_Numbers/MT19937.cpp.o
lib_nova_tools.a: CMakeFiles/_nova_tools.dir/src/Utilities/Driver.cpp.o
lib_nova_tools.a: CMakeFiles/_nova_tools.dir/src/Utilities/Example.cpp.o
lib_nova_tools.a: CMakeFiles/_nova_tools.dir/src/Utilities/File_Utilities.cpp.o
lib_nova_tools.a: CMakeFiles/_nova_tools.dir/src/Utilities/Pthread_Queue.cpp.o
lib_nova_tools.a: CMakeFiles/_nova_tools.dir/src/Utilities/Timer.cpp.o
lib_nova_tools.a: CMakeFiles/_nova_tools.dir/src/Vectors/Vector_2D.cpp.o
lib_nova_tools.a: CMakeFiles/_nova_tools.dir/src/Vectors/Vector_3D.cpp.o
lib_nova_tools.a: CMakeFiles/_nova_tools.dir/build.make
lib_nova_tools.a: CMakeFiles/_nova_tools.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/yanshi/Workspace/CS523/HW2/Problem1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_16) "Linking CXX static library lib_nova_tools.a"
	$(CMAKE_COMMAND) -P CMakeFiles/_nova_tools.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/_nova_tools.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/_nova_tools.dir/build: lib_nova_tools.a

.PHONY : CMakeFiles/_nova_tools.dir/build

CMakeFiles/_nova_tools.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/_nova_tools.dir/cmake_clean.cmake
.PHONY : CMakeFiles/_nova_tools.dir/clean

CMakeFiles/_nova_tools.dir/depend:
	cd /home/yanshi/Workspace/CS523/HW2/Problem1/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/yanshi/Workspace/CS523/HW2/Problem1 /home/yanshi/Workspace/CS523/HW2/Problem1 /home/yanshi/Workspace/CS523/HW2/Problem1/build /home/yanshi/Workspace/CS523/HW2/Problem1/build /home/yanshi/Workspace/CS523/HW2/Problem1/build/CMakeFiles/_nova_tools.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/_nova_tools.dir/depend

