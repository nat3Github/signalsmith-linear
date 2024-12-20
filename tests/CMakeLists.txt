cmake_minimum_required(VERSION 3.24)
project(fft2-tests VERSION 1.0.0)
set(CMAKE_CXX_STANDARD 11)

include(FetchContent)

# DSP library
FetchContent_Declare(
	signalsmith_dsp
	GIT_REPOSITORY https://github.com/Signalsmith-Audio/dsp.git
	GIT_TAG main
	GIT_SHALLOW ON
	SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/others/dsp"
)
FetchContent_MakeAvailable(signalsmith_dsp)

# Plot
FetchContent_Declare(
	signalsmith_plot
	GIT_REPOSITORY https://github.com/Signalsmith-Audio/plot.git
	GIT_TAG main
	GIT_SHALLOW ON
	SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/plot"
)
FetchContent_MakeAvailable(signalsmith_plot)

add_subdirectory(".." ${CMAKE_CURRENT_BINARY_DIR}/fft2)

add_executable(fft2-tests ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
target_link_libraries(fft2-tests INTERFACE fft2)