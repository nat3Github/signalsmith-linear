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
template<class Op>
struct OpVoidArBrp {
	Op op;

	template<typename V>
	void reference(RunData<V> &data) const {
		auto a = data.real(0), b = data.positive(0);
		for (size_t i = 0; i < data.size; ++i) {
			op.op(a[i], b[i]);
		}
	}

	template<typename L, typename V>
	void linear(L &linear, RunData<V> &data) const {
		op.op(linear.wrap(data.real(0), data.size), linear.wrap(data.positive(0), data.size));
	}
};
template<class Op>
struct OpVoidArBc {
	Op op;

	template<typename V>
	void reference(RunData<V> &data) const {
		auto a = data.real(0);
		auto b = data.complex(0);
		for (size_t i = 0; i < data.size; ++i) {
			op.op(a[i], b[i]);
		}
	}

	template<typename L, typename V>
	void linear(L &linear, RunData<V> &data) const {
		op.op(linear.wrap(data.real(0), data.size), linear.wrap(data.complex(0), data.size));
	}
};
template<class Op>
struct OpVoidArBrCr {
	Op op;

	template<typename V>
	void reference(RunData<V> &data) const {
		auto a = data.real(0), b = data.real(1), c = data.real(2);
		for (size_t i = 0; i < data.size; ++i) {
			op.op(a[i], b[i], c[i]);
		}
	}

	template<typename L, typename V>
	void linear(L &linear, RunData<V> &data) const {
		auto a = data.real(0), b = data.real(1), c = data.real(2);
		op.op(linear.wrap(a, data.size), linear.wrap(b, data.size), linear.wrap(c, data.size));
	}
};

struct TestLinear {
	static const int bigPlotWidth = 350, bigPlotHeight = 250;
	static const int tinyPlotWidth = 80, tinyPlotHeight = 80;
	int tinyPlotIndex = 0;
	static const int tinyPlotColumns = 8;
	signalsmith::plot::Plot2D *firstTinyPlot = nullptr;
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
		tinyFigure.style.titlePadding = 0;
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
		for (size_t i = plotName.size(); i < 12; ++i) std::cout << " ";
		std::cout << plotName;
		std::cout << " : " << opWithArgs.op.name << "\n";

		auto &tinyPlot = tinyFigure(tinyPlotIndex%tinyPlotColumns, tinyPlotIndex/tinyPlotColumns).plot(tinyPlotWidth, tinyPlotHeight);
		tinyPlot.title(plotName, 0.5, -1);
		if (!firstTinyPlot) {
			tinyPlot.x.blank();
			tinyPlot.y.blank().major(0, "");
			firstTinyPlot = &tinyPlot;
		} else {
			tinyPlot.x.linkFrom(firstTinyPlot->x);
			tinyPlot.y.linkFrom(firstTinyPlot->y);
		}
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
			std::cout << "\tn = " << n << "\r" << std::flush;
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
			if (dataDouble.distance(refDataDouble) > 1e-12*n) {
				std::cout << "double: n = " << n << ", error = " << dataDouble.distance(refDataDouble) << "\n";
				std::cout << "\nReference:\n";
				refDataDouble.log();
				std::cout << "\nLinear:\n";
				dataDouble.log();
				abort();
			}
			if (dataFloat.distance(refDataFloat) > 1e-6*n) {
				std::cout << "float: n = " << n << ", error = " << dataFloat.distance(refDataFloat) << "\n";
				std::cout << "\nReference:\n";
				refDataFloat.log();
				std::cout << "\nLinear:\n";
				dataFloat.log();
				abort();
			}

			double refTime = n*1e-8;
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
		std::cout << "\t                                \r" << std::flush;
		
		if (benchmarkSeconds) {
			opPlot.write("op-" + plotName + ".svg");
		}
	}
};

#define TEST_EXPR2(Name, expr) \
struct Name { \
	const char *name = #expr; \
	template<class A, class B> \
	void op(A &&a, B &&b) const { \
		expr; \
	} \
};
#define TEST_EXPR3(Name, expr) \
struct Name { \
	const char *name = #expr; \
	template<class A, class B, class C> \
	void op(A &&a, B &&b, C &&c) const { \
		expr; \
	} \
};
TEST_EXPR2(Assign, a = b);
TEST_EXPR3(Add, a = b + c);
TEST_EXPR3(Sub, a = b - c);
TEST_EXPR3(Mul, a = b*c);
TEST_EXPR3(Div, a = b/c);

#define TEST_REAL_METHOD0(Name, scalarExpr, linearExpr) \
struct Name { \
	const char *name = #linearExpr; \
	void op(float &a, float &b) const { \
		scalarExpr; \
	} \
	void op(double &a, double &b) const { \
		scalarExpr; \
	} \
	template<class A, class B> \
	void op(A &&a, B &&b) const { \
		linearExpr; \
	} \
};
#define TEST_COMPLEX_METHOD0(Name, scalarExpr, linearExpr) \
struct Name { \
	const char *name = #linearExpr; \
	void op(float &a, float &b) const { \
		scalarExpr; \
	} \
	void op(double &a, double &b) const { \
		scalarExpr; \
	} \
	void op(std::complex<float> &a, std::complex<float> &b) const { \
		scalarExpr; \
	} \
	void op(std::complex<double> &a, std::complex<double> &b) const { \
		scalarExpr; \
	} \
	void op(float &a, std::complex<float> &b) const { \
		scalarExpr; \
	} \
	void op(double &a, std::complex<double> &b) const { \
		scalarExpr; \
	} \
	void op(std::complex<float> &a, float &b) const { \
		scalarExpr; \
	} \
	void op(std::complex<double> &a, double &b) const { \
		scalarExpr; \
	} \
	template<class A, class B> \
	void op(A &&a, B &&b) const { \
		linearExpr; \
	} \
};
TEST_REAL_METHOD0(Mod1, a = b - std::floor(b), a = b.mod1());
TEST_COMPLEX_METHOD0(Abs, a = std::abs(b), a = b.abs());
TEST_COMPLEX_METHOD0(Norm, a = std::norm(b), a = b.norm());
TEST_REAL_METHOD0(Exp, a = std::exp(b), a = b.exp());
TEST_REAL_METHOD0(Log, a = std::log(b), a = b.log());
TEST_REAL_METHOD0(Log10, a = std::log10(b), a = b.log10());
TEST_REAL_METHOD0(Sqrt, a = std::sqrt(b), a = b.sqrt());
TEST_COMPLEX_METHOD0(SqrtNorm, a = std::sqrt(std::norm(b)), a = b.norm().sqrt());

void testLinear(int maxSize, double benchmarkSeconds) {
	std::cout << "\nExpressions\n-----------\n";
	TestLinear test(maxSize, benchmarkSeconds);
	test.addOp<OpVoidArBr<Assign>>("AssignR");
	test.addOp<OpVoidArBrCr<Add>>("AddR");
	test.addOp<OpVoidArBrCr<Sub>>("SubR");
	test.addOp<OpVoidArBrCr<Mul>>("MulR");
	test.addOp<OpVoidArBrCr<Div>>("DivR");
	test.addOp<OpVoidArBr<Mod1>>("Mod1");
	test.addOp<OpVoidArBr<Abs>>("AbsR");
	test.addOp<OpVoidArBc<Abs>>("AbsC");
	test.addOp<OpVoidArBc<Norm>>("NormC");
	test.addOp<OpVoidArBc<SqrtNorm>>("SqrtNorm");
	test.addOp<OpVoidArBr<Exp>>("ExpR");
	test.addOp<OpVoidArBrp<Log>>("LogR");
	test.addOp<OpVoidArBrp<Log10>>("Log10R");
	test.addOp<OpVoidArBrp<Sqrt>>("Sqrt");
}
