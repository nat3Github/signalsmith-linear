main: out/main
	@cd out && ./main
	
out/main: *.cpp *.h ../*.h
	mkdir -p out
	g++ -std=c++11 -Ofast \
		-Wall -Wextra -Wfatal-errors -Wpedantic -pedantic-errors \
		-framework Accelerate -DSIGNALSMITH_USE_ACCELERATE \
		main.cpp -o out/main

cmake:
	cmake -B build -G Xcode
	cmake --build build --config Release
	mkdir -p out
	cp build/fft2-tests out
	cd out && ./fft2-tests

clean:
	rm -rf out
