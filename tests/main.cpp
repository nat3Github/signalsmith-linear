#include "./test-linear.h"
#include "./test-ffts.h"

#include <cstdlib>
#include <iostream>

int main(int argc, char *argv[]) {
#ifdef SIGNALSMITH_USE_ACCELERATE
	std::cout << u8"✅ SIGNALSMITH_USE_ACCELERATE\n";
#else
	std::cout << u8"❌ SIGNALSMITH_USE_ACCELERATE\n";
#endif
#ifdef SIGNALSMITH_USE_IPP
	std::cout << u8"✅ SIGNALSMITH_USE_IPP\n";
#else
	std::cout << u8"❌ SIGNALSMITH_USE_IPP\n";
#endif
#ifdef SIGNALSMITH_USE_CBLAS
	std::cout << u8"✅ SIGNALSMITH_USE_CBLAS\n";
#else
	std::cout << u8"❌ SIGNALSMITH_USE_CBLAS\n";
#endif
#ifdef __FAST_MATH__
	std::cout << u8"✅ __FAST_MATH__\n";
#else
	std::cout << u8"❌ __FAST_MATH__\n";
#endif

	int maxSize = 8192;
	double benchmarkSeconds = 0;
	if (argc > 1 && !std::strcmp(argv[1], "benchmark")) {
		maxSize = 65536*8;
		benchmarkSeconds = 0.05;
		if (argc > 2) {
			benchmarkSeconds = std::strtod(argv[2], nullptr);
		}
		
	}
	testLinear(maxSize, benchmarkSeconds);
	testFfts(maxSize, benchmarkSeconds);
}
