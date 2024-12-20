cmake_minimum_required(VERSION 3.24)
project(fft2-tests VERSION 1.0.0)
set(CMAKE_CXX_STANDARD 11)

add_subdirectory(".." ${CMAKE_CURRENT_BINARY_DIR}/fft2)

add_executable(fft2-tests ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
target_link_libraries(fft2-tests PUBLIC fft2)