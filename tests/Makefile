test: out/main
	@cd out && ./main

benchmark: out/main
	@cd out && ./main benchmark

benchmark-fast: out/main
	@cd out && ./main benchmark 0.005

# MAC config
out/main: *.cpp *.h ../*.h ../platform/*.h
	mkdir -p out
	g++ -std=c++11 -Ofast \
		-Wall -Wextra -Wfatal-errors -Wpedantic -pedantic-errors \
		-framework Accelerate -DSIGNALSMITH_USE_ACCELERATE -DSIGNALSMITH_USE_CBLAS \
		main.cpp -o out/main

cmake:
	cmake -B build -G Xcode
	cmake --build build --config Release
	mkdir -p out-cmake
	cp build/Release/fft2-tests out-cmake
	cd out-cmake && ./fft2-tests benchmark

clean:
	rm -rf out
