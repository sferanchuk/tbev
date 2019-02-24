
cmake_minimum_required(VERSION 2.8)

PROJECT(TubesFromSplines)

find_package(VTK REQUIRED)
find_package(Qt5Widgets)
include(${VTK_USE_FILE})

add_executable(vtkivew MACOSX_BUNDLE vtkview.cxx)

if(VTK_LIBRARIES)
  target_link_libraries(vtkview ${VTK_LIBRARIES})
else()
  target_link_libraries(vtkview vtkHybrid vtkWidgets)
endif()
