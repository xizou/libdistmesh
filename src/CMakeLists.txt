include_directories("${PROJECT_SOURCE_DIR}/include")

add_library(distmesh SHARED
    distmesh.cpp
    distance_function.cpp
    functional.cpp
    triangulation.cpp
    utils.cpp
)

find_package(Eigen3 REQUIRED)
find_package(Qhull REQUIRED)

target_link_libraries(distmesh Eigen3::Eigen qhull)
