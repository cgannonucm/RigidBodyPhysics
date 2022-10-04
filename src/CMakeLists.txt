cmake_minimum_required(VERSION 3.22.1)

project(RigidBodyPhysics)

add_executable(main.exe main.cpp)

add_subdirectory(core)

target_link_directories(main.exe PUBLIC core/)

target_link_libraries(main.exe core)

#add_subdirectory(dependencies/autodiff)
#
#target_link_directories(main.exe PUBLIC dependencies/autodiff/)
#
#target_link_libraries(main.exe autodiff)
#
#target_link_libraries(core autodiff)