# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/felipe/Desktop/Rotation/poisonM

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/felipe/Desktop/Rotation/poisonM/build

# Include any dependencies generated for this target.
include CMakeFiles/poisson.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/poisson.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/poisson.dir/flags.make

CMakeFiles/poisson.dir/poisson.cpp.o: CMakeFiles/poisson.dir/flags.make
CMakeFiles/poisson.dir/poisson.cpp.o: ../poisson.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/felipe/Desktop/Rotation/poisonM/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/poisson.dir/poisson.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/poisson.dir/poisson.cpp.o -c /home/felipe/Desktop/Rotation/poisonM/poisson.cpp

CMakeFiles/poisson.dir/poisson.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/poisson.dir/poisson.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/felipe/Desktop/Rotation/poisonM/poisson.cpp > CMakeFiles/poisson.dir/poisson.cpp.i

CMakeFiles/poisson.dir/poisson.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/poisson.dir/poisson.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/felipe/Desktop/Rotation/poisonM/poisson.cpp -o CMakeFiles/poisson.dir/poisson.cpp.s

CMakeFiles/poisson.dir/poisson.cpp.o.requires:

.PHONY : CMakeFiles/poisson.dir/poisson.cpp.o.requires

CMakeFiles/poisson.dir/poisson.cpp.o.provides: CMakeFiles/poisson.dir/poisson.cpp.o.requires
	$(MAKE) -f CMakeFiles/poisson.dir/build.make CMakeFiles/poisson.dir/poisson.cpp.o.provides.build
.PHONY : CMakeFiles/poisson.dir/poisson.cpp.o.provides

CMakeFiles/poisson.dir/poisson.cpp.o.provides.build: CMakeFiles/poisson.dir/poisson.cpp.o


# Object files for target poisson
poisson_OBJECTS = \
"CMakeFiles/poisson.dir/poisson.cpp.o"

# External object files for target poisson
poisson_EXTERNAL_OBJECTS =

poisson: CMakeFiles/poisson.dir/poisson.cpp.o
poisson: CMakeFiles/poisson.dir/build.make
poisson: /usr/lib/x86_64-linux-gnu/libmpfr.so
poisson: /usr/lib/x86_64-linux-gnu/libgmp.so
poisson: /usr/lib/x86_64-linux-gnu/libCGAL.so.11.0.1
poisson: /usr/lib/x86_64-linux-gnu/libboost_thread.so
poisson: /usr/lib/x86_64-linux-gnu/libboost_system.so
poisson: /usr/lib/x86_64-linux-gnu/libpthread.so
poisson: /usr/lib/x86_64-linux-gnu/libCGAL.so.11.0.1
poisson: /usr/lib/x86_64-linux-gnu/libboost_thread.so
poisson: /usr/lib/x86_64-linux-gnu/libboost_system.so
poisson: /usr/lib/x86_64-linux-gnu/libpthread.so
poisson: CMakeFiles/poisson.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/felipe/Desktop/Rotation/poisonM/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable poisson"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/poisson.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/poisson.dir/build: poisson

.PHONY : CMakeFiles/poisson.dir/build

CMakeFiles/poisson.dir/requires: CMakeFiles/poisson.dir/poisson.cpp.o.requires

.PHONY : CMakeFiles/poisson.dir/requires

CMakeFiles/poisson.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/poisson.dir/cmake_clean.cmake
.PHONY : CMakeFiles/poisson.dir/clean

CMakeFiles/poisson.dir/depend:
	cd /home/felipe/Desktop/Rotation/poisonM/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/felipe/Desktop/Rotation/poisonM /home/felipe/Desktop/Rotation/poisonM /home/felipe/Desktop/Rotation/poisonM/build /home/felipe/Desktop/Rotation/poisonM/build /home/felipe/Desktop/Rotation/poisonM/build/CMakeFiles/poisson.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/poisson.dir/depend

