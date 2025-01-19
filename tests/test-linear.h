#include "../linear.h"

#include "./test-runner.h"

#include <type_traits>

template<class Op>
struct OpVoidArBr {
	Op op;

	template<typename V>
	void reference(RunData<V> &data) const {
		for (size_t i = 0; i < data.size; ++i) {
			op.op(data.real(0)[i], data.real(1)[i]);
		}
	}

	template<typename L, typename V>
	void linear(L &linear, RunData<V> &data) const {
		op.op(linear.wrap(data.real(0), data.size), linear.wrap(data.real(1), data.size));
	}
};

struct TestLinear {
	static const int bigPlotWidth = 400, bigPlotHeight = 250;
	int maxSize;
	double benchmarkSeconds;

	static double nToX(int n) {
		return (n >= 1) ? std::log(n) + 1 : n;
	}

	void addTicks(signalsmith::plot::Plot2D &plot) {
		plot.x.major(nToX(1), "1");
		for (int n = 2; n <= maxSize; n *= 2) {
			int log = std::round(std::log2(n));
			std::string direct = std::to_string(n);
			std::string sci = "2^" + std::to_string(log);
			std::string &label = (direct.size() <= sci.size()) ? direct : sci;
			plot.x.tick(nToX(n), label);
		}
	}

	signalsmith::plot::Figure figure;
	struct Config {
		std::string name;
		signalsmith::plot::Plot2D &plot;
	};
	std::vector<Config> configs;
	void addConfig(int x, int y, const char *name) {
		configs.push_back({name, figure(x, y).plot(bigPlotWidth, bigPlotHeight)});
		auto &firstPlot = configs[0].plot;
		auto &newPlot = configs.back().plot;
		if (configs.size() == 1) {
			addTicks(firstPlot);
		} else {
			newPlot.x.linkFrom(firstPlot.x);
			newPlot.y.linkFrom(firstPlot.y);
		}
	}
	signalsmith::plot::Legend *legend;
	
	TestLinear(int maxSize, double benchmarkSeconds) : maxSize(maxSize), benchmarkSeconds(benchmarkSeconds) {
		addConfig(0, 0, "float (ref)");
		addConfig(1, 0, "float (Linear)");
		addConfig(0, 1, "double (ref)");
		addConfig(1, 1, "double (Linear)");
		legend = &configs[1].plot.legend(2, 1);
	}
	~TestLinear() {
		if (benchmarkSeconds) {
			figure.write("linear-comparison.svg");
		}
	}

	template<class OpWithArgs>
	void addOp(std::string plotName) {
		OpWithArgs opWithArgs;
	
		signalsmith::plot::Plot2D opPlot(bigPlotWidth, bigPlotHeight);
		auto &opLegend = opPlot.legend(2, 1);
		opPlot.y.major(0);
		opPlot.x.major(0).label(opWithArgs.op.name);
		struct OpLine {
			signalsmith::plot::Line2D &main;
			signalsmith::plot::Line2D &op;
			
			void add(int n, double y) {
				double x = nToX(n);
				main.add(x, y);
				op.add(x, y);
			}
		};
		std::vector<OpLine> opLines;
		while (opLines.size() < configs.size()) {
			auto &config = configs[opLines.size()];
			auto &opPlotLine = opPlot.line();
			opLegend.add(opPlotLine, config.name);
			opLines.push_back({opPlotLine, config.plot.line()});
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

			{
				RunData<float> copy = dataFloat;
				opLines[0].add(n, runBenchmark(benchmarkSeconds, [&](){
					opWithArgs.reference(copy);
				}));
			}
			{
				RunData<float> copy = dataFloat;
				opLines[1].add(n, runBenchmark(benchmarkSeconds, [&](){
					opWithArgs.linear(linearFloat, copy);
				}));
			}
			{
				RunData<double> copy = dataDouble;
				opLines[2].add(n, runBenchmark(benchmarkSeconds, [&](){
					opWithArgs.reference(copy);
				}));
			}
			{
				RunData<double> copy = dataDouble;
				opLines[3].add(n, runBenchmark(benchmarkSeconds, [&](){
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
