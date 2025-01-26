#include "../linear.h"

#include "./test-runner.h"

#include <type_traits>

template<class Op, int typeA, int typeB>
struct OpVoidAB {
	Op op;

	template<class Data>
	void reference(Data &data) const {
		auto a = data.template get<typeA>(0);
		auto b = data.template get<typeB>(1);
		for (size_t i = 0; i < data.size; ++i) {
			op.opRef(a[i], b[i]);
		}
	}

	template<class L, class Data>
	void linear(L &linear, Data &data) const {
		auto a = data.template get<typeA>(0);
		auto b = data.template get<typeB>(1);
		op.op(linear.wrap(a, data.size), linear.wrap(b, data.size));
	}
};
template<class Op>
using OpVoidArBr = OpVoidAB<Op, 0, 0>;
template<class Op>
using OpVoidArBp = OpVoidAB<Op, 0, 1>;
template<class Op>
using OpVoidArBc = OpVoidAB<Op, 0, 2>;
template<class Op>
using OpVoidAcBr = OpVoidAB<Op, 2, 0>;
template<class Op>
using OpVoidAcBp = OpVoidAB<Op, 2, 1>;
template<class Op>
using OpVoidAcBc = OpVoidAB<Op, 2, 2>;

template<class Op, int typeA, int typeB, int typeC>
struct OpVoidABC {
	Op op;

	template<class Data>
	void reference(Data &data) const {
		auto a = data.template get<typeA>(0);
		auto b = data.template get<typeB>(1);
		auto c = data.template get<typeC>(2);
		for (size_t i = 0; i < data.size; ++i) {
			op.opRef(a[i], b[i], c[i]);
		}
	}

	template<class L, class Data>
	void linear(L &linear, Data &data) const {
		auto a = data.template get<typeA>(0);
		auto b = data.template get<typeB>(1);
		auto c = data.template get<typeC>(2);
		size_t size = data.size;
		op.op(linear(a, size), linear(b, size), linear(c, size));
	}
};
template<class Op>
using OpVoidArBrCr = OpVoidABC<Op, 0, 0, 0>;
template<class Op>
using OpVoidArBrCc = OpVoidABC<Op, 0, 0, 2>;
template<class Op>
using OpVoidArBcCr = OpVoidABC<Op, 0, 2, 0>;
template<class Op>
using OpVoidArBcCc = OpVoidABC<Op, 0, 2, 2>;
template<class Op>
using OpVoidAcBrCr = OpVoidABC<Op, 2, 0, 0>;
template<class Op>
using OpVoidAcBrCc = OpVoidABC<Op, 2, 0, 2>;
template<class Op>
using OpVoidAcBcCr = OpVoidABC<Op, 2, 2, 0>;
template<class Op>
using OpVoidAcBcCc = OpVoidABC<Op, 2, 2, 2>;

template<class Op, int typeA, int typeB, int typeC, int typeD>
struct OpVoidABCD {
	Op op;

	template<class Data>
	void reference(Data &data) const {
		auto a = data.template get<typeA>(0);
		auto b = data.template get<typeB>(1);
		auto c = data.template get<typeC>(2);
		auto d = data.template get<typeD>(3);
		for (size_t i = 0; i < data.size; ++i) {
			op.opRef(a[i], b[i], c[i], d[i]);
		}
	}

	template<class L, class Data>
	void linear(L &linear, Data &data) const {
		auto a = data.template get<typeA>(0);
		auto b = data.template get<typeB>(1);
		auto c = data.template get<typeC>(2);
		auto d = data.template get<typeD>(3);
		size_t size = data.size;
		op.op(linear(a, size), linear(b, size), linear(c, size), linear(d, size));
	}
};

template<class Op>
using OpVoidArBrCrDr = OpVoidABCD<Op, 0, 0, 0, 0>;
template<class Op>
using OpVoidArBrCrDc = OpVoidABCD<Op, 0, 0, 0, 2>;
template<class Op>
using OpVoidArBrCcDr = OpVoidABCD<Op, 0, 0, 2, 0>;
template<class Op>
using OpVoidArBrCcDc = OpVoidABCD<Op, 0, 0, 2, 2>;
template<class Op>
using OpVoidArBcCrDr = OpVoidABCD<Op, 0, 2, 0, 0>;
template<class Op>
using OpVoidArBcCrDc = OpVoidABCD<Op, 0, 2, 0, 2>;
template<class Op>
using OpVoidArBcCcDr = OpVoidABCD<Op, 0, 2, 2, 0>;
template<class Op>
using OpVoidArBcCcDc = OpVoidABCD<Op, 0, 2, 2, 2>;
template<class Op>
using OpVoidAcBrCrDr = OpVoidABCD<Op, 2, 0, 0, 0>;
template<class Op>
using OpVoidAcBrCrDc = OpVoidABCD<Op, 2, 0, 0, 2>;
template<class Op>
using OpVoidAcBrCcDr = OpVoidABCD<Op, 2, 0, 2, 0>;
template<class Op>
using OpVoidAcBrCcDc = OpVoidABCD<Op, 2, 0, 2, 2>;
template<class Op>
using OpVoidAcBcCrDr = OpVoidABCD<Op, 2, 2, 0, 0>;
template<class Op>
using OpVoidAcBcCrDc = OpVoidABCD<Op, 2, 2, 0, 2>;
template<class Op>
using OpVoidAcBcCcDr = OpVoidABCD<Op, 2, 2, 2, 0>;
template<class Op>
using OpVoidAcBcCcDc = OpVoidABCD<Op, 2, 2, 2, 2>;

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

	std::string prefix;
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
	
	TestLinear(int maxSize, double benchmarkSeconds, std::string name="linear") : maxSize(maxSize), benchmarkSeconds(benchmarkSeconds), prefix(name) {
		std::cout << name << ":\n";
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
			figure.write(prefix + "-all.svg");
			tinyFigure.write(prefix + "-comparison.svg");
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
			opPlot.write(prefix + "-expr-" + plotName + ".svg");
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
	{
		TestLinear test(maxSize, benchmarkSeconds, "real");
		
		test.addOp<OpVoidArBr<Assign>>("Assign");
		test.addOp<OpVoidArBr<Abs>>("Abs");
		
		test.addOp<OpVoidArBrCr<Add>>("Add");
		test.addOp<OpVoidArBrCr<Sub>>("Sub");
		test.addOp<OpVoidArBrCr<Mul>>("Mul");
		test.addOp<OpVoidArBrCr<Div>>("Div");

		test.addOp<OpVoidArBr<Exp>>("Exp");
		test.addOp<OpVoidArBr<Exp2>>("Exp2");
		test.addOp<OpVoidArBp<Log>>("Log");
		test.addOp<OpVoidArBp<Log2>>("Log2");
		test.addOp<OpVoidArBp<Log10>>("Log10");
		test.addOp<OpVoidArBp<Sqrt>>("Sqrt");
		test.addOp<OpVoidArBr<Cbrt>>("Cbrt");
		test.addOp<OpVoidArBr<Floor>>("Floor");
		test.addOp<OpVoidArBr<MinusFloor>>("MinusFloor");
	}
	{
		TestLinear test(maxSize, benchmarkSeconds, "complex");

		test.addOp<OpVoidAcBcCc<Add>>("Add");
		test.addOp<OpVoidAcBcCc<Sub>>("Sub");
		test.addOp<OpVoidAcBcCc<Mul>>("Mul");
		test.addOp<OpVoidAcBcCc<Div>>("Div");

		test.addOp<OpVoidAcBc<Exp>>("Exp");
		test.addOp<OpVoidAcBc<Log>>("Log");
		test.addOp<OpVoidAcBc<Log10>>("Log10");
	}
	{
		TestLinear test(maxSize, benchmarkSeconds, "mixed");
		test.addOp<OpVoidAcBr<Assign>>("AssignCR");
		
		test.addOp<OpVoidAcBrCr<Add>>("AddCRR");
		test.addOp<OpVoidAcBcCr<Add>>("AddCCR");
		test.addOp<OpVoidAcBrCc<Add>>("AddCRC");
		
		test.addOp<OpVoidAcBrCr<Sub>>("SubCRR");
		test.addOp<OpVoidAcBcCr<Sub>>("SubCCR");
		test.addOp<OpVoidAcBrCc<Sub>>("SubCRC");

		test.addOp<OpVoidAcBrCr<Mul>>("MulCRR");
		test.addOp<OpVoidAcBcCr<Mul>>("MulCCR");
		test.addOp<OpVoidAcBrCc<Mul>>("MulCRC");

		test.addOp<OpVoidAcBrCr<Div>>("DivCRR");
		test.addOp<OpVoidAcBcCr<Div>>("DivCCR");
		test.addOp<OpVoidAcBrCc<Div>>("DivCRC");
	}

	{
		TestLinear test(maxSize, benchmarkSeconds, "complex2real");
		test.addOp<OpVoidArBc<Norm>>("NormC");
		test.addOp<OpVoidArBc<Abs>>("AbsC", "for complex b");
		test.addOp<OpVoidArBc<SqrtNorm>>("SqrtNormC");
		test.addOp<OpVoidAcBc<Conj>>("ConjC");
		test.addOp<OpVoidArBc<Real>>("Real");
		test.addOp<OpVoidArBc<Imag>>("Imag");
		test.addOp<OpVoidArBc<Arg>>("Arg");
	}
	{
		TestLinear test(maxSize, benchmarkSeconds, "compound");
		
		test.addOp<OpVoidArBrCrDr<MulAdd>>("MulAddR");
		test.addOp<OpVoidArBrCrDr<MulAdd2>>("MulAddR2");
		test.addOp<OpVoidArBrCrDr<AddMul>>("AddMulR");
		test.addOp<OpVoidArBrCrDr<AddMul2>>("AddMul2R");
		test.addOp<OpVoidArBrCrDr<MulSub>>("MulSubR");
		test.addOp<OpVoidArBrCrDr<MulSub2>>("MulSubR2");
		test.addOp<OpVoidArBrCrDr<SubMul>>("SubMulR");
		test.addOp<OpVoidArBrCrDr<SubMul2>>("SubMul2R");
	}
}
