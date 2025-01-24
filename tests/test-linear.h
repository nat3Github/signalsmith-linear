#include "../linear.h"

#include "./test-runner.h"

#include <type_traits>

template<class Op>
struct OpVoidArBr {
	Op op;

	template<class Data>
	void reference(Data &data) const {
		auto a = data.real(0), b = data.real(1);
		for (size_t i = 0; i < data.size; ++i) {
			op.opRef(a[i], b[i]);
		}
	}

	template<class L, class Data>
	void linear(L &linear, Data &data) const {
		op.op(linear.wrap(data.real(0), data.size), linear.wrap(data.real(1), data.size));
	}
};
template<class Op>
struct OpVoidArBrp {
	Op op;

	template<class Data>
	void reference(Data &data) const {
		auto a = data.real(0), b = data.positive(0);
		for (size_t i = 0; i < data.size; ++i) {
			op.opRef(a[i], b[i]);
		}
	}

	template<class L, class Data>
	void linear(L &linear, Data &data) const {
		op.op(linear.wrap(data.real(0), data.size), linear.wrap(data.positive(0), data.size));
	}
};
template<class Op>
struct OpVoidArBc {
	Op op;

	template<class Data>
	void reference(Data &data) const {
		auto a = data.real(0);
		auto b = data.complex(0);
		for (size_t i = 0; i < data.size; ++i) {
			op.opRef(a[i], b[i]);
		}
	}

	template<class L, class Data>
	void linear(L &linear, Data &data) const {
		op.op(linear.wrap(data.real(0), data.size), linear.wrap(data.complex(0), data.size));
	}
};
template<class Op>
struct OpVoidAcBc {
	Op op;

	template<class Data>
	void reference(Data &data) const {
		auto a = data.complex(0);
		auto b = data.complex(1);
		for (size_t i = 0; i < data.size; ++i) {
			op.opRef(a[i], b[i]);
		}
	}

	template<class L, class Data>
	void linear(L &linear, Data &data) const {
		op.op(linear.wrap(data.complex(0), data.size), linear.wrap(data.complex(1), data.size));
	}
};
template<class Op>
struct OpVoidArBrCr {
	Op op;

	template<class Data>
	void reference(Data &data) const {
		auto a = data.real(0), b = data.real(1), c = data.real(2);
		for (size_t i = 0; i < data.size; ++i) {
			op.opRef(a[i], b[i], c[i]);
		}
	}

	template<class L, class Data>
	void linear(L &linear, Data &data) const {
		auto a = data.real(0), b = data.real(1), c = data.real(2);
		op.op(linear.wrap(a, data.size), linear.wrap(b, data.size), linear.wrap(c, data.size));
	}
};
template<class Op>
struct OpVoidArBrCrDr {
	Op op;

	template<class Data>
	void reference(Data &data) const {
		auto a = data.real(0), b = data.real(1), c = data.real(2), d = data.real(3);
		for (size_t i = 0; i < data.size; ++i) {
			op.opRef(a[i], b[i], c[i], d[i]);
		}
	}

	template<class L, class Data>
	void linear(L &linear, Data &data) const {
		auto a = data.real(0), b = data.real(1), c = data.real(2), d = data.real(3);
		op.op(linear.wrap(a, data.size), linear.wrap(b, data.size), linear.wrap(c, data.size), linear.wrap(d, data.size));
	}
};
template<class Op>
struct OpVoidArBrCrDrEr {
	Op op;

	template<class Data>
	void reference(Data &data) const {
		auto a = data.real(0), b = data.real(1), c = data.real(2), d = data.real(3), e = data.real(4);
		for (size_t i = 0; i < data.size; ++i) {
			op.opRef(a[i], b[i], c[i], d[i], e[i]);
		}
	}

	template<class L, class Data>
	void linear(L &linear, Data &data) const {
		auto a = data.real(0), b = data.real(1), c = data.real(2), d = data.real(3), e = data.real(4);
		op.op(linear.wrap(a, data.size), linear.wrap(b, data.size), linear.wrap(c, data.size), linear.wrap(d, data.size), linear.wrap(e, data.size));
	}
};

struct TestLinear {
	static const int bigPlotWidth = 300, bigPlotHeight = 250;
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
		plot.y.blank().major(0, "");
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
		addConfig(2, 0, "float (Linear)");
		addConfig(1, 0, "float (fallback)");
		addConfig(0, 0, "float (ref)");
		addConfig(2, 1, "double (Linear)");
		addConfig(1, 1, "double (fallback)");
		addConfig(0, 1, "double (ref)");
		legend = &configs[0].plot.legend(2, 1);
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
	void addOp(std::string plotName, std::string nameSuffix="") {
		OpWithArgs opWithArgs;
		std::string opName = opWithArgs.op.name;
		if (nameSuffix.size() > 0) {
			opName += " [" + nameSuffix + "]";
		}
		for (size_t i = plotName.size(); i < 12; ++i) std::cout << " ";
		std::cout << plotName;
		std::cout << " : " << opName << "\n";

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
		opPlot.title(opName);
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
		legend->add(opLines[0].main, opName);
		
		signalsmith::linear::Linear linear;
		signalsmith::linear::LinearImpl<false> linearFallback;
		auto runSize = [&](int n){
			std::cout << "\tn = " << n << "\r" << std::flush;
			linear.reserve<float>(n);
			linear.reserve<double>(n);
			linearFallback.reserve<float>(n);
			linearFallback.reserve<double>(n);

			RunData<double> dataDouble(n);
			RunData<float> dataFloat(n);

			RunData<double> refDataDouble = dataDouble;
			RunData<float> refDataFloat = dataFloat;
			opWithArgs.reference(refDataDouble);
			opWithArgs.reference(refDataFloat);
			{
				RunData<double> copyDouble(n);
				RunData<float> copyFloat(n);
				opWithArgs.linear(linear, copyDouble);
				opWithArgs.linear(linear, copyFloat);
				if (copyDouble.distance(refDataDouble) > 1e-12*n) {
					std::cout << "Linear double: n = " << n << ", error = " << copyDouble.distance(refDataDouble) << "\n";
					std::cout << "\nReference:\n";
					refDataDouble.log();
					std::cout << "\nLinear:\n";
					copyDouble.log();
					abort();
				}
				if (copyFloat.distance(refDataFloat) > 1e-6*n) {
					std::cout << "Linear float: n = " << n << ", error = " << copyFloat.distance(refDataFloat) << "\n";
					std::cout << "\nReference:\n";
					refDataFloat.log();
					std::cout << "\nLinear:\n";
					copyFloat.log();
					abort();
				}
			}
			{
				RunData<double> copyDouble(n);
				RunData<float> copyFloat(n);
				opWithArgs.linear(linearFallback, copyDouble);
				opWithArgs.linear(linearFallback, copyFloat);
				if (copyDouble.distance(refDataDouble) > 1e-12*n) {
					std::cout << "fallback double: n = " << n << ", error = " << copyDouble.distance(refDataDouble) << "\n";
					std::cout << "\nReference:\n";
					refDataDouble.log();
					std::cout << "\nLinear:\n";
					copyDouble.log();
					abort();
				}
				if (copyFloat.distance(refDataFloat) > 1e-6*n) {
					std::cout << "fallback float: n = " << n << ", error = " << copyFloat.distance(refDataFloat) << "\n";
					std::cout << "\nReference:\n";
					refDataFloat.log();
					std::cout << "\nLinear:\n";
					copyFloat.log();
					abort();
				}
			}
			
			double refTime = n*1e-8;
			{
				RunData<float> copy = dataFloat;
				opLines[0].add(n, refTime*runBenchmark(benchmarkSeconds, [&](){
					opWithArgs.linear(linear, copy);
				}));
			}
			{
				RunData<float> copy = dataFloat;
				opLines[1].add(n, refTime*runBenchmark(benchmarkSeconds, [&](){
					opWithArgs.linear(linearFallback, copy);
				}));
			}
			{
				RunData<float> copy = dataFloat;
				opLines[2].add(n, refTime*runBenchmark(benchmarkSeconds, [&](){
					opWithArgs.reference(copy);
				}));
			}
			{
				RunData<double> copy = dataDouble;
				opLines[3].add(n, refTime*runBenchmark(benchmarkSeconds, [&](){
					opWithArgs.linear(linear, copy);
				}));
			}
			{
				RunData<double> copy = dataDouble;
				opLines[4].add(n, refTime*runBenchmark(benchmarkSeconds, [&](){
					opWithArgs.linear(linearFallback, copy);
				}));
			}
			{
				RunData<double> copy = dataDouble;
				opLines[5].add(n, refTime*runBenchmark(benchmarkSeconds, [&](){
					opWithArgs.reference(copy);
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

#define TEST_EXPR2(Name, refExpr, expr) \
struct Name { \
	const char *name = #expr; \
	template<class A, class B> \
	void opRef(A &a, B &b) const { \
		refExpr; \
	} \
	template<class A, class B> \
	void op(A &&a, B &&b) const { \
		expr; \
	} \
};
#define TEST_EXPR3(Name, refExpr, expr) \
struct Name { \
	const char *name = #expr; \
	template<class A, class B, class C> \
	void opRef(A &a, B &b, C &c) const { \
		refExpr; \
	} \
	template<class A, class B, class C> \
	void op(A &&a, B &&b, C &&c) const { \
		expr; \
	} \
};
#define TEST_EXPR4(Name, refExpr, expr) \
struct Name { \
	const char *name = #expr; \
	template<class A, class B, class C, class D> \
	void opRef(A &a, B &b, C &c, D &d) const { \
		refExpr; \
	} \
	template<class A, class B, class C, class D> \
	void op(A &&a, B &&b, C &&c, D &&d) const { \
		expr; \
	} \
};
#define TEST_EXPR5(Name, refExpr, expr) \
struct Name { \
	const char *name = #expr; \
	template<class A, class B, class C, class D, class E> \
	void opRef(A &a, B &b, C &c, D &d, E &e) const { \
		refExpr; \
	} \
	template<class A, class B, class C, class D, class E> \
	void op(A &&a, B &&b, C &&c, D &&d, E &&e) const { \
		expr; \
	} \
};
TEST_EXPR2(Assign, a = b, a = b);
TEST_EXPR3(Add, a = b + c, a = b + c);
TEST_EXPR3(Sub, a = b - c, a = b - c);
TEST_EXPR3(Mul, a = b*c, a = b*c);
TEST_EXPR3(Div, a = b/c, a = b/c);

TEST_EXPR2(Abs, a = std::abs(b), a = b.abs());
TEST_EXPR2(Norm, a = std::norm(b), a = b.norm());
TEST_EXPR2(Exp, a = std::exp(b), a = b.exp());
TEST_EXPR2(Exp2, a = std::exp2(b), a = b.exp2());
TEST_EXPR2(Log, a = std::log(b), a = b.log());
TEST_EXPR2(Log2, a = std::log2(b), a = b.log2());
TEST_EXPR2(Log10, a = std::log10(b), a = b.log10());
TEST_EXPR2(Sqrt, a = std::sqrt(b), a = b.sqrt());
TEST_EXPR2(Cbrt, a = std::cbrt(b), a = b.cbrt());
TEST_EXPR2(SqrtNorm, a = std::sqrt(std::norm(b)), a = b.norm().sqrt());
TEST_EXPR2(Conj, a = std::conj(b), a = b.conj());
TEST_EXPR2(Real, a = std::real(b), a = b.real());
TEST_EXPR2(Imag, a = std::imag(b), a = b.imag());
TEST_EXPR2(Arg, a = std::arg(b), a = b.arg());
TEST_EXPR2(Floor, a = std::floor(b), a = b.floor());
TEST_EXPR2(MinusFloor, a = b - std::floor(b), a = b - b.floor());

TEST_EXPR4(MulAdd, a = b*c + d, a = b*c + d);
TEST_EXPR4(MulAdd2, a = b + c*d, a = b + c*d);
TEST_EXPR4(AddMul, a = b*(c + d), a = b*(c + d));
TEST_EXPR4(AddMul2, a = (b + c)*d, a = (b + c)*d);
TEST_EXPR4(MulSub, a = b*c - d, a = b*c - d);
TEST_EXPR4(MulSub2, a = b - c*d, a = b - c*d);
TEST_EXPR4(SubMul, a = b*(c - d), a = b*(c - d));
TEST_EXPR4(SubMul2, a = (b - c)*d, a = (b - c)*d);

void testLinear(int maxSize, double benchmarkSeconds) {
	std::cout << "\nExpressions\n-----------\n";
	TestLinear test(maxSize, benchmarkSeconds);

	test.addOp<OpVoidArBrCrDr<MulAdd>>("MulAddR");
	test.addOp<OpVoidArBrCrDr<MulAdd2>>("MulAddR2");
	test.addOp<OpVoidArBrCrDr<AddMul>>("AddMulR");
	test.addOp<OpVoidArBrCrDr<AddMul2>>("AddMul2R");
	test.addOp<OpVoidArBrCrDr<MulSub>>("MulSubR");
	test.addOp<OpVoidArBrCrDr<MulSub2>>("MulSubR2");
	test.addOp<OpVoidArBrCrDr<SubMul>>("SubMulR");
	test.addOp<OpVoidArBrCrDr<SubMul2>>("SubMul2R");

	test.addOp<OpVoidArBr<Assign>>("AssignR");
	test.addOp<OpVoidArBrCr<Add>>("AddR");
	test.addOp<OpVoidArBrCr<Sub>>("SubR");
	test.addOp<OpVoidArBrCr<Mul>>("MulR");
	test.addOp<OpVoidArBrCr<Div>>("DivR");
	test.addOp<OpVoidArBr<Abs>>("AbsR");
	test.addOp<OpVoidArBc<Abs>>("AbsC", "for complex b");
	test.addOp<OpVoidArBc<Norm>>("NormC");
	test.addOp<OpVoidArBc<SqrtNorm>>("SqrtNormC");
	test.addOp<OpVoidAcBc<Conj>>("ConjC");
	test.addOp<OpVoidArBc<Real>>("Real");
	test.addOp<OpVoidArBc<Imag>>("Imag");
	test.addOp<OpVoidArBc<Arg>>("Arg");
	test.addOp<OpVoidArBr<Exp>>("ExpR");
	test.addOp<OpVoidArBr<Exp2>>("Exp2R");
	test.addOp<OpVoidArBrp<Log>>("LogR");
	test.addOp<OpVoidArBrp<Log2>>("Log2R");
	test.addOp<OpVoidArBrp<Log10>>("Log10R");
	test.addOp<OpVoidArBrp<Sqrt>>("Sqrt");
	test.addOp<OpVoidArBr<Cbrt>>("Cbrt");
	test.addOp<OpVoidArBr<Floor>>("Floor");
	test.addOp<OpVoidArBr<MinusFloor>>("MinusFloor");
}
