#include "../linear.h"

#include "./test-runner.h"

#include <type_traits>

int strideMax(int a, int b) {
	return std::max(std::abs(a), std::abs(b));
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
	void run(RunData<V> &data, int strideIn, int strideOut) {
		this->linear.copy(data.size/strideMax(strideIn, strideOut), data.rpA, strideIn, data.rpB, strideOut);
	}
};

template<typename V, bool useLinear=true>
struct RealNorm2 : public BaseOp<V, useLinear> {
	void run(RunData<V> &data, int strideIn, int strideOut) {
		data.rA = this->linear.norm2(data.size/strideMax(strideIn, strideOut), data.rpA, strideIn);
	}
};
template<typename V, bool useLinear=true>
struct ComplexNorm2 : public BaseOp<V, useLinear> {
	void run(RunData<V> &data, int strideIn, int strideOut) {
		data.rA = this->linear.norm2(data.size/strideMax(strideIn, strideOut), data.cpA, strideIn);
	}
};

template<typename V, bool useLinear=true>
struct ComplexNorm2Real : public BaseOp<V, useLinear> {
	void run(RunData<V> &data, int strideIn, int strideOut) {
		this->linear.norm2(data.size/strideMax(strideIn, strideOut), data.cpA, strideIn, data.rpA, strideOut);
	}
};

template<typename V, bool useLinear=true>
struct SplitNorm2 : public BaseOp<V, useLinear> {
	void run(RunData<V> &data, int strideIn, int strideOut) {
		data.rA = this->linear.norm2(data.size/strideMax(strideIn, strideOut), data.sA, strideIn);
	}
};

template<typename V, bool useLinear=true>
struct SplitNorm2Real : public BaseOp<V, useLinear> {
	void run(RunData<V> &data, int strideIn, int strideOut) {
		this->linear.norm2(data.size/strideMax(strideIn, strideOut), data.sA, strideIn, data.rpA, strideOut);
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
		auto forDouble = plot.runner<Op<double, false>>("for-loop (double)");
		auto forFloat = plot.runner<Op<float, false>>("for-loop (float)");

		auto runSize = [&](int n){
			int strideN = n/strideMax(strideIn, strideOut);
			double refTime = 1e-8*strideN;

			RunData<double> dataDouble(n);
			RunData<float> dataFloat(n);

			RunData<double> refDataDouble = dataDouble;
			RunData<float> refDataFloat = dataFloat;
			forDouble.wrapper.run(refDataDouble, si, so);
			forFloat.wrapper.run(refDataFloat, si, so);

			linearDouble.run(dataDouble, refTime, &refDataDouble, si, so);
			double comparisonRate = linearFloat.run(dataFloat, refTime, &refDataFloat, si, so);
			forDouble.run(dataDouble, refTime, &refDataDouble, si, so);
			forFloat.run(dataFloat, refTime, &refDataFloat, si, so);
			
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
	test.opStrides<SplitNorm2>("norm2s");
	test.opStrides<ComplexNorm2Real>("norm2cr");
	test.opStrides<SplitNorm2Real>("norm2sr");
}
