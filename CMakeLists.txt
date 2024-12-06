cmake_minimum_required(VERSION 3.16)
project(pompeiu_hausdorff)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
SET(CMAKE_BUILD_TYPE Release)

# option to build executable (default true)
option(BUILD_EXECUTABLE "Build executable" ON)

# Libigl
include(libigl)

file(GLOB SRCFILES "${PROJECT_SOURCE_DIR}/*.cpp")

# library called phd (pompeiu hausdorff distance)
set(LIBRARY_NAME phd)
add_library(${LIBRARY_NAME} STATIC
  src/pompeiu_hausdorff.cpp
  src/upper_bounds.cpp
  src/kang_upper_bound.cpp
  src/kang_intersect_edge_and_bisector.cpp
  src/bisector_of_two_points.cpp
  src/gettimeofday.cpp)
target_link_libraries(${LIBRARY_NAME} igl::core)

if(BUILD_EXECUTABLE)
  # executable called pompeiu_hausdorff
  add_executable(${PROJECT_NAME} main.cpp)
  target_link_libraries(${PROJECT_NAME} ${LIBRARY_NAME} igl::core)
  # change name to pompeiu_hausdorff
  set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME pompeiu_hausdorff)
endif()
  
  
