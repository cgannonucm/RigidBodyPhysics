cmake_minimum_required(VERSION 3.22.1)

project(RigidBodyPhysics)

add_executable(main.exe main.cpp)

add_subdirectory(extern/Eigen)

set(AUTODIFF_BUILD_TESTS false)

set(AUTODIFF_BUILD_PYTHON false)

set(AUTODIFF_BUILD_EXAMPLES false)

add_subdirectory(extern/autodiff)

#add_subdirectory(core/)

target_include_directories(main.exe
    PRIVATE extern/autodiff
    PUBLIC extern/Eigen
    )

target_link_libraries(main.exe Eigen autodiff:autodiff)
#
#target_include_directories(core
#    PUBLIC extern/autodiff
#    PUBLIC extern/Eigen
#    )
