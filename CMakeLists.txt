cmake_minimum_required(VERSION 2.6)
project(mcgp)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

### libIGL options: choose between header only and compiled static library
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_EMBREE            "Use Embree"         OFF)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          OFF)

# i need to use an embree build with avx2
# find_package(embree 3 REQUIRED)
# find_package(LIBIGL REQUIRED)
add_subdirectory(external/my-libigl)

# Add your project files
include_directories("include/")
file(GLOB SRCFILES src/*.cpp)


add_library(core ${SRCFILES})
if(LIBIGL_WITH_EMBREE)
    target_link_libraries(core igl::core igl::opengl igl::opengl_glfw igl::embree)
else()
    target_link_libraries(core igl::core igl::opengl igl::opengl_glfw)
endif()


add_subdirectory(tutorials)