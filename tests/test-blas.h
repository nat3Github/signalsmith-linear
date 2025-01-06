#include "../cppblas.h"

#include "./test-runner.h"

static constexpr double blasTestSeconds = 0.005;

template<template<int, int, bool, bool> class Op, int strideIn=1, int strideOut=1>
void testBlasOp(std::string opName, std::string strideName, signalsmith::plot::Plot2D &comparisonPlot) {
	RunPlot plot("blas-" + opName + "-stride-" + strideName, blasTestSeconds);
	auto &comparisonLine = comparisonPlot.line();

	auto blasDouble = plot.runner<Op<strideIn, strideOut, true, true>>("BLAS (double)");
	auto blasFloat = plot.runner<Op<strideIn, strideOut, true, true>>("BLAS (float)");
	auto forDouble = plot.runner<Op<strideIn, strideOut, false, false>>("for-dynamic (double)");
	auto forFloat = plot.runner<Op<strideIn, strideOut, false, false>>("for-dynamic (float)");
	auto forCDouble = plot.runner<Op<strideIn, strideOut, true, false>>("for-constexpr (double)");
	auto forCFloat = plot.runner<Op<strideIn, strideOut, true, false>>("for-constexpr (float)");

	int maxSize = 65536*8;
	auto runSize = [&](int n){
		int strideN = n/std::abs(blasDouble.wrapper.strideMax);
		double refTime = 1e-8*strideN;

		RunData<double> dataDouble(n);
		RunData<float> dataFloat(n);

		RunData<double> refDataDouble = dataDouble;
		RunData<float> refDataFloat = dataFloat;
		forDouble.wrapper.run(refDataDouble);
		forFloat.wrapper.run(refDataFloat);

		blasDouble.run(dataDouble, refTime, &refDataDouble);
		double comparisonRate = blasFloat.run(dataFloat, refTime, &refDataFloat);
		forDouble.run(dataDouble, refTime, &refDataDouble);
		forFloat.run(dataFloat, refTime, &refDataFloat);
		forCDouble.run(dataDouble, refTime, &refDataDouble);
		forCFloat.run(dataFloat, refTime, &refDataFloat);
		
		comparisonLine.add(std::log(n), comparisonRate);
	};
	for (int n = 1; n <= maxSize; n *= 2) {
		if (n >= 8) runSize(int(n*2/M_PI));
		if (n >= 4) runSize(int(n*M_PI/4));
		runSize(n);
		plot.tick(n);
	}
}

template<int strideIn=1, int strideOut=1, bool constexprStride=false, bool useBlas=true>
struct Copy {
	static constexpr int strideMax = (strideIn > strideOut) ? strideIn : strideOut;

	void prepare(int /*size*/, int /*maxSize*/) {}
	
	template<class V>
	void run(RunData<V> &data) {
		if (constexprStride) {
			signalsmith::blas::copy<V, useBlas>(data.size/strideMax, data.rpA, strideIn, data.rpB, strideOut);
		} else {
			volatile int si = strideIn;
			volatile int so = strideOut;
			signalsmith::blas::copy<V, useBlas>(data.size/strideMax, data.rpA, si, data.rpB, so);
		}
	}
};

template<int strideIn=1, int strideOut=1, bool constexprStride=false, bool useBlas=true>
struct RealNorm2 {
	static constexpr int strideMax = strideIn;

	void prepare(int /*size*/, int /*maxSize*/) {}
	
	template<class V>
	void run(RunData<V> &data) {
		if (constexprStride) {
			data.rA = signalsmith::blas::norm2<V, useBlas>(data.size/strideMax, data.rpA, strideIn);
		} else {
			volatile int si = strideIn;
			data.rA = signalsmith::blas::norm2<V, useBlas>(data.size/strideMax, data.rpA, si);
		}
	}
};
template<int strideIn=1, int strideOut=1, bool constexprStride=false, bool useBlas=true>
struct ComplexNorm2 {
	static constexpr int strideMax = strideIn;

	void prepare(int /*size*/, int /*maxSize*/) {}
	
	template<class V>
	void run(RunData<V> &data) {
		if (constexprStride) {
			data.rA = signalsmith::blas::norm2<V, useBlas>(data.size/strideMax, data.cpA, strideIn);
		} else {
			volatile int si = strideIn;
			data.rA = signalsmith::blas::norm2<V, useBlas>(data.size/strideMax, data.cpA, si);
		}
	}
};

template<template<int, int, bool, bool> class Op>
void testBlasOpStrides(std::string name, signalsmith::plot::Plot2D &plot) {
	plot.y.major(0);
	plot.x.major(0).label("copy");
	
	testBlasOp<Op, 1, 1>(name, "1", plot);
	testBlasOp<Op, 2, 2>(name, "2", plot);
	testBlasOp<Op, -3, 5>(name, "m3-5", plot);
	testBlasOp<Op, 13, 7>(name, "13-7", plot);
	testBlasOp<Op, 4, -4>(name, "4-m4", plot);
}

void testBlas() {
	signalsmith::plot::Figure figure;
	testBlasOpStrides<Copy>("copy", figure(0, 0).plot(200, 100));
	testBlasOpStrides<RealNorm2>("norm2r", figure(1, 0).plot(200, 100));
	testBlasOpStrides<ComplexNorm2>("norm2c", figure(1, 0).plot(200, 100));
	
	figure.write("blas-comparison.svg");
}
