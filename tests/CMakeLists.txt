cmake_minimum_required(VERSION 3.24)
project(linear-tests VERSION 1.0.0)
set(CMAKE_CXX_STANDARD 11)

include(FetchContent)

# DSP library
FetchContent_Declare(
	signalsmith_dsp
	GIT_REPOSITORY https://github.com/Signalsmith-Audio/dsp.git
	GIT_TAG 5c7d1f3eb375b4862b682d310ccfbaafb6a4477e # main
	GIT_SHALLOW ON
	SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/others/dsp"
)
FetchContent_MakeAvailable(signalsmith_dsp)

# Plot
FetchContent_Declare(
	signalsmith_plot
	GIT_REPOSITORY https://github.com/Signalsmith-Audio/plot.git
	GIT_TAG fd68b5c6dd928d5b18b6d401c88bd7957021d405 # main
	GIT_SHALLOW ON
	SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/plot"
)
FetchContent_MakeAvailable(signalsmith_plot)

add_subdirectory(".." ${CMAKE_CURRENT_BINARY_DIR}/signalsmith-linear)

add_executable(tests ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
target_link_libraries(tests PRIVATE signalsmith-linear)
