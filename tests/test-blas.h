#include "./test-runner.h"

#include "../cppblas.h"

void testBlas() {
	RunPlot plot("blas");

	int maxSize = 65536*8;
	auto runSize = [&](int n){
		double refTime = 1e-8*n; // heuristic for expected computation time, just to compare different sizes
	};
	for (int n = 1; n <= maxSize; n *= 2) {
		if (n/8) runSize(n/2 - 1);
		if (n/8) runSize(n/2 + 3);
		runSize(n);
	}
}
