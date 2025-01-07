#include "../linear.h"

#include "./test-runner.h"

#include <type_traits>

// like std::abs() and std::max(std::abs()...), but constexpr in C++11
static constexpr int maxAbs(int a) {
	return (a >= 0) ? a : -a;
}
static constexpr int maxAbs(int a, int b) {
	return (a >= 0) ? (
		(b >= 0) ? (
			(a > b) ? a : b
		) : (
			(a > -b) ? a : -b
		)
	) : (
		(b >= 0) ? (
			(-a > b) ? -a : b
		) : (
			(-a > -b) ? -a : -b
		)
	);
}


template<typename V, int strideIn=1, int strideOut=1, bool constexprStride=false, bool useLinear=true>
struct BaseOp {
	using Linear = signalsmith::linear::Linear<V, !useLinear>;
	Linear linear;

	void prepare(int /*size*/, int maxSize) {
		linear.reserve(maxSize);
	}
};

template<typename V, int strideIn=1, int strideOut=1, bool constexprStride=false, bool useLinear=true>
struct Copy : public BaseOp<V, strideIn, strideOut, constexprStride, useLinear> {
	static constexpr int strideMax = maxAbs(strideIn, strideOut);

	void run(RunData<V> &data) {
		if (constexprStride) {
			this->linear.copy(data.size/strideMax, data.rpA, strideIn, data.rpB, strideOut);
		} else {
			volatile int si = strideIn;
			volatile int so = strideOut;
			this->linear.copy(data.size/strideMax, data.rpA, si, data.rpB, so);
		}
	}
};

template<typename V, int strideIn=1, int strideOut=1, bool constexprStride=false, bool useLinear=true>
struct RealNorm2 : public BaseOp<V, strideIn, strideOut, constexprStride, useLinear> {
	static constexpr int strideMax = maxAbs(strideIn);

	void run(RunData<V> &data) {
		if (constexprStride) {
			data.rA = this->linear.norm2(data.size/strideMax, data.rpA, strideIn);
		} else {
			volatile int si = strideIn;
			data.rA = this->linear.norm2(data.size/strideMax, data.rpA, si);
		}
	}
};
template<typename V, int strideIn=1, int strideOut=1, bool constexprStride=false, bool useLinear=true>
struct ComplexNorm2 : public BaseOp<V, strideIn, strideOut, constexprStride, useLinear> {
	static constexpr int strideMax = maxAbs(strideIn);

	void run(RunData<V> &data) {
		if (constexprStride) {
			data.rA = this->linear.norm2(data.size/strideMax, data.cpA, strideIn);
		} else {
			volatile int si = strideIn;
			data.rA = this->linear.norm2(data.size/strideMax, data.cpA, si);
		}
	}
};

template<typename V, int strideIn=1, int strideOut=1, bool constexprStride=false, bool useLinear=true>
struct ComplexNorm2Real : public BaseOp<V, strideIn, strideOut, constexprStride, useLinear> {
	static constexpr int strideMax = maxAbs(strideIn, strideOut);

	void run(RunData<V> &data) {
		if (constexprStride) {
			this->linear.norm2(data.size/strideMax, data.cpA, strideIn, data.rpA, strideOut);
		} else {
			volatile int si = strideIn;
			this->linear.norm2(data.size/strideMax, data.cpA, si, data.rpA, strideOut);
		}
	}
};

struct TestLinear {
	int maxSize;
	double benchmarkSeconds;

	int plotIndex = 0;
	int plotColumns = 3;
	signalsmith::plot::Figure figure;
	
	TestLinear(int maxSize, double benchmarkSeconds) : maxSize(maxSize), benchmarkSeconds(benchmarkSeconds) {}
	~TestLinear() {
		if (benchmarkSeconds) {
			figure.write("linear-comparison.svg");
		}
	}

	template<template<typename, int, int, bool, bool> class Op>
	void opStrides(std::string name) {
		auto &plot = figure(plotIndex%plotColumns, plotIndex/plotColumns).plot(200, 100);
		++plotIndex;
		
		plot.y.major(0);
		plot.x.major(0).label(name);
		
		testLinearOp<Op, 1, 1>(name, "1", plot);
		testLinearOp<Op, 2, 2>(name, "2", plot);
		testLinearOp<Op, -3, 5>(name, "m3-5", plot);
		testLinearOp<Op, 13, 7>(name, "13-7", plot);
		testLinearOp<Op, 4, -4>(name, "4-m4", plot);
		testLinearOp<Op, 16, 1>(name, "16-1", plot);
	}

	template<template<typename, int, int, bool, bool> class Op, int strideIn=1, int strideOut=1>
	void testLinearOp(std::string opName, std::string strideName, signalsmith::plot::Plot2D &comparisonPlot) {
		RunPlot plot("linear-" + opName + "-stride-" + strideName, benchmarkSeconds);
		auto &comparisonLine = comparisonPlot.line();

		auto linearDouble = plot.runner<Op<double, strideIn, strideOut, false, true>>("linear (double)");
		auto linearFloat = plot.runner<Op<float, strideIn, strideOut, false, true>>("linear (float)");
		auto forDouble = plot.runner<Op<double, strideIn, strideOut, false, false>>("for-dynamic (double)");
		auto forFloat = plot.runner<Op<float, strideIn, strideOut, false, false>>("for-dynamic (float)");
		auto forCDouble = plot.runner<Op<double, strideIn, strideOut, true, false>>("for-constexpr (double)");
		auto forCFloat = plot.runner<Op<float, strideIn, strideOut, true, false>>("for-constexpr (float)");

		auto runSize = [&](int n){
			int strideN = n/std::abs(linearDouble.wrapper.strideMax);
			double refTime = 1e-8*strideN;

			RunData<double> dataDouble(n);
			RunData<float> dataFloat(n);

			RunData<double> refDataDouble = dataDouble;
			RunData<float> refDataFloat = dataFloat;
			forDouble.wrapper.run(refDataDouble);
			forFloat.wrapper.run(refDataFloat);

			linearDouble.run(dataDouble, refTime, &refDataDouble);
			double comparisonRate = linearFloat.run(dataFloat, refTime, &refDataFloat);
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
};
void testLinear(int maxSize, double benchmarkSeconds) {
	TestLinear test(maxSize, benchmarkSeconds);
	test.opStrides<Copy>("copy");
	test.opStrides<RealNorm2>("norm2r");
	test.opStrides<ComplexNorm2>("norm2c");
	test.opStrides<ComplexNorm2Real>("norm2cr");
}
