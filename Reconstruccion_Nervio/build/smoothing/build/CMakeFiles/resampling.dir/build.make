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
CMAKE_SOURCE_DIR = /home/felipe/Desktop/Reconstrucción/build/smoothing

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/felipe/Desktop/Reconstrucción/build/smoothing/build

# Include any dependencies generated for this target.
include CMakeFiles/resampling.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/resampling.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/resampling.dir/flags.make

CMakeFiles/resampling.dir/resampling.cpp.o: CMakeFiles/resampling.dir/flags.make
CMakeFiles/resampling.dir/resampling.cpp.o: ../resampling.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/felipe/Desktop/Reconstrucción/build/smoothing/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/resampling.dir/resampling.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/resampling.dir/resampling.cpp.o -c /home/felipe/Desktop/Reconstrucción/build/smoothing/resampling.cpp

CMakeFiles/resampling.dir/resampling.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/resampling.dir/resampling.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/felipe/Desktop/Reconstrucción/build/smoothing/resampling.cpp > CMakeFiles/resampling.dir/resampling.cpp.i

CMakeFiles/resampling.dir/resampling.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/resampling.dir/resampling.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/felipe/Desktop/Reconstrucción/build/smoothing/resampling.cpp -o CMakeFiles/resampling.dir/resampling.cpp.s

CMakeFiles/resampling.dir/resampling.cpp.o.requires:

.PHONY : CMakeFiles/resampling.dir/resampling.cpp.o.requires

CMakeFiles/resampling.dir/resampling.cpp.o.provides: CMakeFiles/resampling.dir/resampling.cpp.o.requires
	$(MAKE) -f CMakeFiles/resampling.dir/build.make CMakeFiles/resampling.dir/resampling.cpp.o.provides.build
.PHONY : CMakeFiles/resampling.dir/resampling.cpp.o.provides

CMakeFiles/resampling.dir/resampling.cpp.o.provides.build: CMakeFiles/resampling.dir/resampling.cpp.o


# Object files for target resampling
resampling_OBJECTS = \
"CMakeFiles/resampling.dir/resampling.cpp.o"

# External object files for target resampling
resampling_EXTERNAL_OBJECTS =

resampling: CMakeFiles/resampling.dir/resampling.cpp.o
resampling: CMakeFiles/resampling.dir/build.make
resampling: /usr/lib/x86_64-linux-gnu/libboost_system.so
resampling: /usr/lib/x86_64-linux-gnu/libboost_filesystem.so
resampling: /usr/lib/x86_64-linux-gnu/libboost_thread.so
resampling: /usr/lib/x86_64-linux-gnu/libboost_date_time.so
resampling: /usr/lib/x86_64-linux-gnu/libboost_iostreams.so
resampling: /usr/lib/x86_64-linux-gnu/libboost_serialization.so
resampling: /usr/lib/x86_64-linux-gnu/libboost_chrono.so
resampling: /usr/lib/x86_64-linux-gnu/libboost_atomic.so
resampling: /usr/lib/x86_64-linux-gnu/libboost_regex.so
resampling: /usr/lib/x86_64-linux-gnu/libpthread.so
resampling: /usr/lib/libpcl_common.so
resampling: /usr/lib/x86_64-linux-gnu/libflann_cpp.so
resampling: /usr/lib/libpcl_kdtree.so
resampling: /usr/lib/libpcl_octree.so
resampling: /usr/lib/libpcl_io.so
resampling: /usr/lib/libpcl_stereo.so
resampling: /usr/lib/libpcl_search.so
resampling: /usr/lib/libpcl_sample_consensus.so
resampling: /usr/lib/libpcl_ml.so
resampling: /usr/lib/libpcl_filters.so
resampling: /usr/lib/libpcl_features.so
resampling: /usr/lib/libpcl_visualization.so
resampling: /usr/lib/libpcl_segmentation.so
resampling: /usr/lib/libpcl_people.so
resampling: /usr/lib/x86_64-linux-gnu/libqhull.so
resampling: /usr/lib/libpcl_surface.so
resampling: /usr/lib/libpcl_outofcore.so
resampling: /usr/lib/libpcl_registration.so
resampling: /usr/lib/libpcl_recognition.so
resampling: /usr/lib/libpcl_tracking.so
resampling: /usr/lib/libpcl_keypoints.so
resampling: /usr/lib/libpcl_apps.so
resampling: /usr/lib/x86_64-linux-gnu/libboost_system.so
resampling: /usr/lib/x86_64-linux-gnu/libboost_filesystem.so
resampling: /usr/lib/x86_64-linux-gnu/libboost_thread.so
resampling: /usr/lib/x86_64-linux-gnu/libboost_date_time.so
resampling: /usr/lib/x86_64-linux-gnu/libboost_iostreams.so
resampling: /usr/lib/x86_64-linux-gnu/libboost_serialization.so
resampling: /usr/lib/x86_64-linux-gnu/libboost_chrono.so
resampling: /usr/lib/x86_64-linux-gnu/libboost_atomic.so
resampling: /usr/lib/x86_64-linux-gnu/libboost_regex.so
resampling: /usr/lib/x86_64-linux-gnu/libpthread.so
resampling: /usr/lib/x86_64-linux-gnu/libqhull.so
resampling: /usr/lib/x86_64-linux-gnu/libflann_cpp.so
resampling: /home/felipe/VTK-build/lib/libvtkPoissonReconstruction-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkPowercrust-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkIOInfovis-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkRenderingContextOpenGL2-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkTestingRendering-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkViewsContext2D-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkFiltersGeneric-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkTestingGenericBridge-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkDomainsChemistryOpenGL2-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkIOAMR-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkIOExodus-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkRenderingVolumeOpenGL2-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkFiltersFlowPaths-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkFiltersHyperTree-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkImagingStencil-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkFiltersParallelImaging-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkFiltersProgrammable-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkFiltersSMP-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkFiltersSelection-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkFiltersVerdict-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkIOParallel-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkFiltersTexture-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkFiltersTopology-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkGeovisCore-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkIOEnSight-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkIOExportOpenGL2-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkInteractionImage-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkIOImport-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkIOLSDyna-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkIOMINC-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkIOMovie-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkIOPLY-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkIOParallelXML-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkIOSQL-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkTestingIOSQL-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkIOTecplotTable-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkIOVideo-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkImagingStatistics-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkRenderingImage-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkImagingMorphological-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkRenderingLOD-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkViewsInfovis-8.1.so.1
resampling: /usr/lib/libpcl_common.so
resampling: /usr/lib/libpcl_kdtree.so
resampling: /usr/lib/libpcl_octree.so
resampling: /usr/lib/libpcl_io.so
resampling: /usr/lib/libpcl_stereo.so
resampling: /usr/lib/libpcl_search.so
resampling: /usr/lib/libpcl_sample_consensus.so
resampling: /usr/lib/libpcl_ml.so
resampling: /usr/lib/libpcl_filters.so
resampling: /usr/lib/libpcl_features.so
resampling: /usr/lib/libpcl_visualization.so
resampling: /usr/lib/libpcl_segmentation.so
resampling: /usr/lib/libpcl_people.so
resampling: /usr/lib/libpcl_surface.so
resampling: /usr/lib/libpcl_outofcore.so
resampling: /usr/lib/libpcl_registration.so
resampling: /usr/lib/libpcl_recognition.so
resampling: /usr/lib/libpcl_tracking.so
resampling: /usr/lib/libpcl_keypoints.so
resampling: /usr/lib/libpcl_apps.so
resampling: /home/felipe/VTK-build/lib/libvtkFiltersPoints-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtklibxml2-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkDomainsChemistry-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkFiltersAMR-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkImagingMath-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkverdict-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkIOGeometry-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkexoIIc-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkFiltersParallel-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkIONetCDF-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtknetcdfcpp-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkjsoncpp-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkproj4-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkIOExport-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkRenderingGL2PSOpenGL2-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkRenderingOpenGL2-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkglew-8.1.so.1
resampling: /usr/lib/x86_64-linux-gnu/libSM.so
resampling: /usr/lib/x86_64-linux-gnu/libICE.so
resampling: /usr/lib/x86_64-linux-gnu/libX11.so
resampling: /usr/lib/x86_64-linux-gnu/libXext.so
resampling: /usr/lib/x86_64-linux-gnu/libXt.so
resampling: /home/felipe/VTK-build/lib/libvtkgl2ps-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtklibharu-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkNetCDF-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkhdf5_hl-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkhdf5-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkoggtheora-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkParallelCore-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkIOLegacy-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtksqlite-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkChartsCore-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkRenderingContext2D-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkViewsCore-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkInteractionWidgets-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkFiltersHybrid-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkInteractionStyle-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkRenderingAnnotation-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkImagingColor-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkRenderingVolume-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkIOXML-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkIOXMLParser-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkIOCore-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtklz4-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkexpat-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkFiltersImaging-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkImagingGeneral-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkImagingSources-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkRenderingLabel-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkRenderingFreeType-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkRenderingCore-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkCommonColor-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkFiltersGeometry-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkfreetype-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkInfovisLayout-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkInfovisCore-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkFiltersExtraction-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkFiltersStatistics-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkImagingFourier-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkalglib-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkFiltersModeling-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkFiltersSources-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkFiltersGeneral-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkCommonComputationalGeometry-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkFiltersCore-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkImagingHybrid-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkImagingCore-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkIOImage-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkCommonExecutionModel-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkCommonDataModel-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkCommonMisc-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkCommonSystem-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtksys-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkCommonTransforms-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkCommonMath-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkCommonCore-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkDICOMParser-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkmetaio-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkpng-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtktiff-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkzlib-8.1.so.1
resampling: /home/felipe/VTK-build/lib/libvtkjpeg-8.1.so.1
resampling: /usr/lib/x86_64-linux-gnu/libm.so
resampling: CMakeFiles/resampling.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/felipe/Desktop/Reconstrucción/build/smoothing/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable resampling"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/resampling.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/resampling.dir/build: resampling

.PHONY : CMakeFiles/resampling.dir/build

CMakeFiles/resampling.dir/requires: CMakeFiles/resampling.dir/resampling.cpp.o.requires

.PHONY : CMakeFiles/resampling.dir/requires

CMakeFiles/resampling.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/resampling.dir/cmake_clean.cmake
.PHONY : CMakeFiles/resampling.dir/clean

CMakeFiles/resampling.dir/depend:
	cd /home/felipe/Desktop/Reconstrucción/build/smoothing/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/felipe/Desktop/Reconstrucción/build/smoothing /home/felipe/Desktop/Reconstrucción/build/smoothing /home/felipe/Desktop/Reconstrucción/build/smoothing/build /home/felipe/Desktop/Reconstrucción/build/smoothing/build /home/felipe/Desktop/Reconstrucción/build/smoothing/build/CMakeFiles/resampling.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/resampling.dir/depend

