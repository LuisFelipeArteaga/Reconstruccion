
cmake_minimum_required(VERSION 2.8)
 
PROJECT(main_Reconstruccion)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(main_Reconstruccion main_Reconstruccion.cxx)
 
target_link_libraries(main_Reconstruccion ${VTK_LIBRARIES})
