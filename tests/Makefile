main: out/main
	@cd out && ./main
	
out/main: *.cpp *.h ../*.h
	mkdir -p out
	g++ -std=c++11 -Ofast \
		-Wall -Wextra -Wfatal-errors -Wpedantic -pedantic-errors \
		-framework Accelerate \
		others/kissfft/kiss_fft.c \
		main.cpp -o out/main

clean:
	rm -rf out
