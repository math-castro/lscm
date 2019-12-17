#cmake version 
cmake_minimum_required(VERSION 3.1)

# CMake
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# project declarations
project(lscm)

# C++ and Compiler
# add_compile_options(-Wall -Wextra)

# set(Eigen3_DIR "~/opt/include/eigen3/Eigen")
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
find_package(LIBIGL REQUIRED QUIET)

# Source Files
include_directories(src)
file(GLOB_RECURSE SRC_FILES "src/*.cpp")
add_executable(lscm ${SRC_FILES})
target_link_libraries(lscm igl::core igl::opengl_glfw)
