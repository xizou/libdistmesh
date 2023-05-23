include_directories(
        "${PROJECT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

set(SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/scripts/plot_mesh.py")

if(WIN32)
    file(COPY ${SCRIPT} DESTINATION "${CMAKE_BINARY_DIR}")
else()
    file(COPY ${SCRIPT} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
endif()

add_executable(ellipse src/ellipse.cpp)
add_executable(polygon src/polygon.cpp)
add_executable(polygons_x_symmetric src/polygons_x_symmetric.cpp)
add_executable(rectangle_with_circular_hole src/rectangle_with_circular_hole.cpp)
add_executable(square src/square.cpp)
add_executable(unit_circle src/unit_circle.cpp)

target_link_libraries(ellipse distmesh)
target_link_libraries(polygon distmesh)
target_link_libraries(polygons_x_symmetric distmesh)
target_link_libraries(rectangle_with_circular_hole distmesh)
target_link_libraries(square distmesh)
target_link_libraries(unit_circle distmesh)
