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


template<typename V, bool useLinear=true>
struct BaseOp {
	using Linear = signalsmith::linear::Linear<V, !useLinear>;
	Linear linear;

	void prepare(int /*size*/, int maxSize) {
		linear.reserve(maxSize);
	}
};

template<typename V, bool useLinear=true>
struct Copy : public BaseOp<V, useLinear> {
	static constexpr int strideMax = maxAbs(strideIn, strideOut);

	void run(RunData<V> &data, int strideIn, int strideOut) {
		this->linear.copy(data.size/strideMax, data.rpA, strideIn, data.rpB, strideOut);
	}
};

template<typename V, bool useLinear=true>
struct RealNorm2 : public BaseOp<V, useLinear> {
	static constexpr int strideMax = maxAbs(strideIn);

	void run(RunData<V> &data, int strideIn, int) {
		data.rA = this->linear.norm2(data.size/strideMax, data.rpA, strideIn);
	}
};
template<typename V, bool useLinear=true>
struct ComplexNorm2 : public BaseOp<V, useLinear> {
	static constexpr int strideMax = maxAbs(strideIn);

	void run(RunData<V> &data, int strideIn, int) {
		data.rA = this->linear.norm2(data.size/strideMax, data.cpA, strideIn);
	}
};

template<typename V, bool useLinear=true>
struct ComplexNorm2Real : public BaseOp<V, strideIn, strideOut, constexprStride, useLinear> {
	static constexpr int strideMax = maxAbs(strideIn, strideOut);

	void run(RunData<V> &data, int strideIn, int strideOut) {
		this->linear.norm2(data.size/strideMax, data.cpA, strideIn, data.rpA, strideOut);
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

	template<template<typename, bool> class Op>
	void opStrides(std::string name) {
		auto &plot = figure(plotIndex%plotColumns, plotIndex/plotColumns).plot(200, 100);
		++plotIndex;
		
		plot.y.major(0);
		plot.x.major(0).label(name);

		testLinearOp<Op>(name, "1", plot, 1, 1);
		testLinearOp<Op>(name, "2", plot, 2, 2);
		testLinearOp<Op>(name, "m3-5", plot, -3, 5);
		testLinearOp<Op>(name, "13-7", plot, 13, 7);
		testLinearOp<Op>(name, "4-m4", plot, 4, -4);
		testLinearOp<Op>(name, "16-1", plot, 16, 1);
	}

	template<template<typename, bool> class Op>
	void testLinearOp(std::string opName, std::string strideName, signalsmith::plot::Plot2D &comparisonPlot, int strideIn=1, int strideOut=1) {
		RunPlot plot("linear-" + opName + "-stride-" + strideName, benchmarkSeconds);
		auto &comparisonLine = comparisonPlot.line();

		volatile int si = strideIn;
		volatile int so = strideOut;

		auto linearDouble = plot.runner<Op<double, true>>("linear (double)");
		auto linearFloat = plot.runner<Op<float, true>>("linear (float)");
		auto forDouble = plot.runner<Op<double, false>>("for-dynamic (double)");
		auto forFloat = plot.runner<Op<float, false>>("for-dynamic (float)");
		auto forCDouble = plot.runner<Op<double, false>>("for-constexpr (double)");
		auto forCFloat = plot.runner<Op<float, false>>("for-constexpr (float)");

		auto runSize = [&](int n){
			int strideN = n/std::abs(linearDouble.wrapper.strideMax);
			double refTime = 1e-8*strideN;

			RunData<double> dataDouble(n);
			RunData<float> dataFloat(n);

			RunData<double> refDataDouble = dataDouble;
			RunData<float> refDataFloat = dataFloat;
			forDouble.wrapper.run(refDataDouble, 0, nullptr, si, so);
			forFloat.wrapper.run(refDataFloat, 0, nullptr, si, so);

			linearDouble.run(dataDouble, refTime, &refDataDouble, si, so);
			double comparisonRate = linearFloat.run(dataFloat, refTime, &refDataFloat, si, so);
			forDouble.run(dataDouble, refTime, &refDataDouble, si, so);
			forFloat.run(dataFloat, refTime, &refDataFloat, si, so);
			forCDouble.run(dataDouble, refTime, &refDataDouble, si, so);
			forCFloat.run(dataFloat, refTime, &refDataFloat, si, so);
			
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
