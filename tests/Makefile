main: out/main
	@cd out && ./main
	
out/main: *.cpp *.h ../*.h
	mkdir -p out
	g++ -std=c++11 -Ofast \
		-Wall -Wextra -Wfatal-errors -Wpedantic -pedantic-errors \
		-framework Accelerate -DSIGNALSMITH_USE_ACCELERATE \
		main.cpp -o out/main

clean:
	rm -rf out
