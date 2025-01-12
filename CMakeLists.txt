add_library(signalsmith-linear INTERFACE)
set_target_properties(signalsmith-linear PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR})

#find_package(BLAS)
#if(BLAS_FOUND)
#	message(STATUS "found BLAS (${BLAS_LIBRARIES})")
#	target_link_options(fft2 INTERFACE ${BLAS_LINKER_FLAGS})
#	target_link_libraries(fft2 INTERFACE BLAS)
#	target_compile_definitions(fft2 INTERFACE SIGNALSMITH_USE_CBLAS)
#else()
#	message(FATAL "no BLAS found")
##	message(STATUS "using OpenBLAS")
##	include(FetchContent)
##	FetchContent_Declare(
##		openblas
##		GIT_REPOSITORY https://github.com/OpenMathLib/OpenBLAS.git
##		GIT_TAG 5ef8b1964658f9cb6a6324a06f6a1a022609b0c5 # 0.3.28
##		GIT_SHALLOW ON
##	)
##	FetchContent_MakeAvailable(openblas)
##	target_link_libraries(fft2 INTERFACE openblas)
##	target_compile_definitions(fft2 INTERFACE SIGNALSMITH_USE_CBLAS)
#endif()

if(WIN32)
	find_package(IPP)
	if(IPP_FOUND)
		message(STATUS "using IPP")
		target_link_libraries(signalsmith-linear INTERFACE IPP::ippcore IPP::ipps IPP::ippi)
		target_compile_definitions(signalsmith-linear INTERFACE SIGNALSMITH_USE_IPP)
	endif()
elseif(APPLE)
	#find_library(Accelerate_LIBRARY NAMES Accelerate)
	message(STATUS "using Accelerate")
	target_link_libraries(signalsmith-linear INTERFACE "-framework Accelerate")
	target_compile_definitions(signalsmith-linear INTERFACE SIGNALSMITH_USE_ACCELERATE)
endif()
