cmake_minimum_required(VERSION 3.16)
project(pompeiu_hausdorff)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
SET(CMAKE_BUILD_TYPE Release)

# Libigl
include(libigl)
igl_include(copyleft/cgal)

file(GLOB SRCFILES "${PROJECT_SOURCE_DIR}/*.cpp")

add_executable(${PROJECT_NAME}
  main.cpp
  src/upper_bounds.cpp
  src/kang_upper_bound.cpp
  src/kang_intersect_edge_and_bisector.cpp
  src/bisector_of_two_points.cpp
  src/gettimeofday.cpp)

target_link_libraries(${PROJECT_NAME} igl::core igl_copyleft::cgal)
