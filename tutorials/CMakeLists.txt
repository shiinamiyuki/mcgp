add_executable(lap3d_simple "lap3d_simple.cpp")
target_link_libraries(lap3d_simple core igl::core igl::opengl igl::opengl_glfw igl::png)

add_executable(lap3d_adaptive "lap3d_adaptive.cpp")
target_link_libraries(lap3d_adaptive core igl::core igl::opengl igl::opengl_glfw igl::png)


add_executable(poi3d_simple "poi3d_simple.cpp")
target_link_libraries(poi3d_simple core igl::core igl::opengl igl::opengl_glfw igl::png)