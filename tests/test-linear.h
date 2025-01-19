#include "../linear.h"

#include "./test-runner.h"

#include <type_traits>

template<class Op>
struct OpVoidArBr {
	Op op;

	template<typename V>
	void reference(RunData<V> &data) const {
		auto a = data.real(0), b = data.real(1);
		for (size_t i = 0; i < data.size; ++i) {
			op.op(a[i], b[i]);
		}
	}

	template<typename L, typename V>
	void linear(L &linear, RunData<V> &data) const {
		op.op(linear.wrap(data.real(0), data.size), linear.wrap(data.real(1), data.size));
	}
};

struct TestLinear {
	static const int bigPlotWidth = 350, bigPlotHeight = 250;
	static const int tinyPlotWidth = 80, tinyPlotHeight = 80;
	int tinyPlotIndex = 0;
	static const int tinyPlotColumns = 8;
	int maxSize;
	double benchmarkSeconds;

	static double nToX(int n) {
		return (n >= 1) ? std::log(n) + 1 : n;
	}

	void addTicks(signalsmith::plot::Plot2D &plot) {
		plot.y.blank().major(0, "").label("relative speed");
		plot.x.major(nToX(1), "1");
		for (int n = 4; n <= maxSize; n *= 4) {
			int log = std::round(std::log2(n));
			std::string direct = std::to_string(n);
			std::string sci = "2^" + std::to_string(log);
			std::string &label = (direct.size() <= sci.size()) ? direct : sci;
			plot.x.tick(nToX(n/2), "");
			plot.x.tick(nToX(n), label);
		}
	}

	signalsmith::plot::Figure figure;
	signalsmith::plot::Figure tinyFigure;
	struct Config {
		std::string name;
		signalsmith::plot::Plot2D &plot;
	};
	std::vector<Config> configs;
	void addConfig(int x, int y, const char *name) {
		auto &newPlot = figure(x, y).plot(bigPlotWidth, bigPlotHeight);
		configs.push_back({name, newPlot});
		auto &firstPlot = configs[0].plot;
		if (configs.size() == 1) {
			addTicks(firstPlot);
		} else {
			newPlot.x.linkFrom(firstPlot.x);
			newPlot.y.linkFrom(firstPlot.y);
		}
		newPlot.title(name);
	}
	signalsmith::plot::Legend *legend;
	
	TestLinear(int maxSize, double benchmarkSeconds) : maxSize(maxSize), benchmarkSeconds(benchmarkSeconds) {
		addConfig(0, 0, "float (ref)");
		addConfig(1, 0, "float (Linear)");
		addConfig(0, 1, "double (ref)");
		addConfig(1, 1, "double (Linear)");
		legend = &configs[1].plot.legend(2, 1);
	}
	TestLinear(const TestLinear &other) = delete;
	~TestLinear() {
		if (benchmarkSeconds) {
			figure.write("linear.svg");
			tinyFigure.write("linear-comparison.svg");
		}
	}

	template<class OpWithArgs>
	void addOp(std::string plotName) {
		OpWithArgs opWithArgs;

		auto &tinyPlot = tinyFigure(tinyPlotIndex%tinyPlotColumns, tinyPlotIndex/tinyPlotColumns).plot(tinyPlotWidth, tinyPlotHeight);
		tinyPlot.x.blank().major(0, "").label(plotName);
		tinyPlot.y.blank().major(0, "");
		++tinyPlotIndex;
		signalsmith::plot::Plot2D opPlot(bigPlotWidth*1.5, bigPlotHeight);
		auto &opLegend = opPlot.legend(2, 1);
		addTicks(opPlot);
		opPlot.title(opWithArgs.op.name);
		struct OpLine {
			signalsmith::plot::Line2D &main;
			signalsmith::plot::Line2D &op;
			signalsmith::plot::Line2D &tiny;
			
			void add(int n, double y) {
				double x = nToX(n);
				main.add(x, y);
				op.add(x, y);
				tiny.add(x, y);
			}
		};
		std::vector<OpLine> opLines;
		while (opLines.size() < configs.size()) {
			auto &config = configs[opLines.size()];
			auto &opPlotLine = opPlot.line();
			opLegend.add(opPlotLine, config.name);
			auto &mainLine = config.plot.line();
			auto &tinyLine = tinyPlot.line();
			opLines.push_back(OpLine{mainLine, opPlotLine, tinyLine});
		}
		legend->add(opLines[0].main, opWithArgs.op.name);
		
		signalsmith::linear::Linear<float> linearFloat;
		signalsmith::linear::Linear<double> linearDouble;
		auto runSize = [&](int n){
			linearFloat.reserve(n);
			linearDouble.reserve(n);

			RunData<double> dataDouble(n);
			RunData<float> dataFloat(n);

			RunData<double> refDataDouble = dataDouble;
			RunData<float> refDataFloat = dataFloat;
			opWithArgs.reference(refDataDouble);
			opWithArgs.reference(refDataFloat);
			opWithArgs.linear(linearDouble, dataDouble);
			opWithArgs.linear(linearFloat, dataFloat);

			double refTime = n;
			{
				RunData<float> copy = dataFloat;
				opLines[0].add(n, refTime*runBenchmark(benchmarkSeconds, [&](){
					opWithArgs.reference(copy);
				}));
			}
			{
				RunData<float> copy = dataFloat;
				opLines[1].add(n, refTime*runBenchmark(benchmarkSeconds, [&](){
					opWithArgs.linear(linearFloat, copy);
				}));
			}
			{
				RunData<double> copy = dataDouble;
				opLines[2].add(n, refTime*runBenchmark(benchmarkSeconds, [&](){
					opWithArgs.reference(copy);
				}));
			}
			{
				RunData<double> copy = dataDouble;
				opLines[3].add(n, refTime*runBenchmark(benchmarkSeconds, [&](){
					opWithArgs.linear(linearDouble, copy);
				}));
			}
		};
		for (int n = 1; n <= maxSize; n *= 2) {
			if (n >= 8) runSize(int(n*2/M_PI));
			if (n >= 4) runSize(int(n*M_PI/4));
			runSize(n);
		}
		
		if (benchmarkSeconds) {
			opPlot.write("op-" + plotName + ".svg");
		}
	}
};

struct AssignR {
	const char *name = "a = b";

	template<class A, class B>
	void op(A &&a, B &&b) const {
		a = b;
	}
};
struct MulR {
	const char *name = "a = b*c";

	template<class A, class B, class C>
	void op(A &&a, B &&b, C &&c) const {
		a = b*c;
	}
};

void testLinear(int maxSize, double benchmarkSeconds) {
	TestLinear test(maxSize, benchmarkSeconds);
	test.addOp<OpVoidArBr<AssignR>>("AssignR");
}
