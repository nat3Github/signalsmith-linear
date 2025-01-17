test: out/main
	@cd out && ./main

dev: out-cmake/tests
	cd out-cmake && ./tests benchmark 0.001 fft

benchmark: out/main
	@cd out && ./main benchmark

benchmark-fast: out/main
	@cd out && ./main benchmark 0.001

# MAC config
out/main: *.cpp *.h ../*.h ../platform/*.h
	mkdir -p out
	g++ -std=c++11 -Ofast \
		-Wall -Wextra -Wfatal-errors -Wpedantic -pedantic-errors \
		-framework Accelerate -DSIGNALSMITH_USE_ACCELERATE -DSIGNALSMITH_USE_CBLAS \
		main.cpp -o out/main

out-cmake/tests: CMakeLists.txt *.cpp *.h ../*.h ../platform/*.h
	cmake -B build -G Xcode
	cmake --build build --config Release
	mkdir -p out-cmake
	cp build/Release/tests out-cmake

cmake: out-cmake/tests
	cd out-cmake && ./tests

cmake-benchmark: out-cmake/tests
	cd out-cmake && ./tests benchmark

cmake-benchmark-fast: out-cmake/tests
	cd out-cmake && ./tests benchmark 0.001

clean:
	rm -rf out
	rm -rf out-cmake
	rm -rf build
