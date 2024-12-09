cmake_minimum_required(VERSION 3.16)
project(pompeiu_hausdorff)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
SET(CMAKE_BUILD_TYPE Release)

# option to build executable (default true)
option(BUILD_EXECUTABLE "Build executable" ON)
option(BUILD_PYTHON_BINDINGS "Build python bindings" OFF)

# Libigl
include(libigl)


set(CMAKE_POSITION_INDEPENDENT_CODE ON)


file(GLOB SRCFILES "${PROJECT_SOURCE_DIR}/*.cpp")

# library called phd (pompeiu hausdorff distance)
set(CMAKE_CXX_STANDARD 11)
set(LIBRARY_NAME phd)
add_library(${LIBRARY_NAME} STATIC
  src/pompeiu_hausdorff.cpp
  src/upper_bounds.cpp
  src/kang_upper_bound.cpp
  src/kang_intersect_edge_and_bisector.cpp
  src/bisector_of_two_points.cpp)
target_link_libraries(${LIBRARY_NAME} igl::core)

if(BUILD_EXECUTABLE)
  # executable called pompeiu_hausdorff
  set(EXECUTABLE_TARGET pompeiu_hausdorff_exe)
  add_executable(${EXECUTABLE_TARGET} main.cpp)
  target_link_libraries(${EXECUTABLE_TARGET} ${LIBRARY_NAME} igl::core)
  # change name to pompeiu_hausdorff
  set_target_properties(${EXECUTABLE_TARGET} PROPERTIES OUTPUT_NAME pompeiu_hausdorff)
endif()
  
  
# Download and set up nanobind
if(BUILD_PYTHON_BINDINGS)
  set(CMAKE_CXX_STANDARD 17)
  
  if (CMAKE_VERSION VERSION_LESS 3.18)
    set(DEV_MODULE Development)
  else()
    set(DEV_MODULE Development.Module)
  endif()

  find_package(Python 3.8 
    REQUIRED COMPONENTS Interpreter ${DEV_MODULE} 
    OPTIONAL_COMPONENTS Development.SABIModule)
  include(FetchContent)
  FetchContent_Declare(
    nanobind
    GIT_REPOSITORY https://github.com/wjakob/nanobind.git
    GIT_TAG        v2.2.0
  )
  FetchContent_MakeAvailable(nanobind)

  nanobind_add_module(
    cascading_upper_bounds_ext 
    STABLE_ABI
    NB_STATIC
    src/bindings.cpp)
  target_link_libraries(
    cascading_upper_bounds_ext PRIVATE ${LIBRARY_NAME} igl::core)
  #nanobind_add_stub(
  #  ${binding_name}_stub
  #  MODULE ${binding_name}
  #  OUTPUT "${binding_name}.pyi"
  #  PYTHON_PATH $<TARGET_FILE_DIR:${binding_name}>
  #  DEPENDS ${binding_name}
  #)

  # Install directive for scikit-build-core
  install(TARGETS cascading_upper_bounds_ext LIBRARY DESTINATION cascading_upper_bounds)

endif()
