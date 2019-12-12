#cmake version 
cmake_minimum_required(VERSION 3.1)

# project declarations
project(lscm)

# C++ and Compiler
add_compile_options(-Wall -Wextra)

# CMake
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
find_package(LIBIGL REQUIRED QUIET)

# Source Files
include_directories(src)
file(GLOB_RECURSE SRC_FILES "src/*.cpp")
add_executable(lscm ${SRC_FILES})
target_link_libraries(lscm igl::core igl::opengl_glfw)
