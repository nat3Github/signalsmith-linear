#ifndef SIGNALSMITH_AUDIO_LINEAR_H
#define SIGNALSMITH_AUDIO_LINEAR_H

#include <cmath>
#include <complex>
#include <array>

#define SIGNALSMITH_AUDIO_LINEAR_CHUNK_SIZE 8
#define SIGNALSMITH_AUDIO_LINEAR_CHUNK_FOREACH_STEP(value, indexName, ...) \
	{ \
		size_t indexName = value; \
		__VA_ARGS__; \
	}
#define SIGNALSMITH_AUDIO_LINEAR_CHUNK_FOREACH(indexName, ...) \
	SIGNALSMITH_AUDIO_LINEAR_CHUNK_FOREACH_STEP(0, indexName, __VA_ARGS__) \
	SIGNALSMITH_AUDIO_LINEAR_CHUNK_FOREACH_STEP(1, indexName, __VA_ARGS__) \
	SIGNALSMITH_AUDIO_LINEAR_CHUNK_FOREACH_STEP(2, indexName, __VA_ARGS__) \
	SIGNALSMITH_AUDIO_LINEAR_CHUNK_FOREACH_STEP(3, indexName, __VA_ARGS__) \
	SIGNALSMITH_AUDIO_LINEAR_CHUNK_FOREACH_STEP(4, indexName, __VA_ARGS__) \
	SIGNALSMITH_AUDIO_LINEAR_CHUNK_FOREACH_STEP(5, indexName, __VA_ARGS__) \
	SIGNALSMITH_AUDIO_LINEAR_CHUNK_FOREACH_STEP(6, indexName, __VA_ARGS__) \
	SIGNALSMITH_AUDIO_LINEAR_CHUNK_FOREACH_STEP(7, indexName, __VA_ARGS__)

namespace signalsmith { namespace linear {

// Unsized pointers
template<typename V>
using ConstRealPointer = const V *;
template<typename V>
using RealPointer = V *;
template<typename V>
using RealChunk = std::array<V, SIGNALSMITH_AUDIO_LINEAR_CHUNK_SIZE>;

template<typename V>
using ConstComplexPointer = const std::complex<V> *;
template<typename V>
using ComplexPointer = std::complex<V> *;
template<typename V>
using ComplexChunk = std::array<std::complex<V>, SIGNALSMITH_AUDIO_LINEAR_CHUNK_SIZE>;

template<typename V>
struct ConstSplitPointer {
	ConstRealPointer<V> real, imag;
	ConstSplitPointer(ConstRealPointer<V> real, ConstRealPointer<V> imag) : real(real), imag(imag) {}
};
template<typename V>
struct SplitPointer {
	RealPointer<V> real, imag;
	SplitPointer(RealPointer<V> real, RealPointer<V> imag) : real(real), imag(imag) {}
	operator ConstSplitPointer<V>() {
		return {real, imag};
	}
};
template<typename V>
struct SplitChunk {
	RealChunk<V> real, imag;
};

#define SIGNALSMITH_LINEAR_SIZED_TYPE(Name) \
	template<typename V> \
	struct Const##Name { \
		Const##Name(Const##Name##Pointer<V> pointer, size_t size) : pointer(pointer), size(size) {} \
		Const##Name##Pointer<V> pointer; \
		const size_t size; \
	}; \
	template<typename V> \
	struct Name { \
		Name(Name##Pointer<V> pointer, size_t size) : pointer(pointer), size(size) {} \
		operator Const##Name<V>() const { \
			return {pointer, size}; \
		} \
		Name##Pointer<V> pointer; \
		const size_t size; \
	};

SIGNALSMITH_LINEAR_SIZED_TYPE(Real)
SIGNALSMITH_LINEAR_SIZED_TYPE(Complex)
SIGNALSMITH_LINEAR_SIZED_TYPE(Split)
#undef SIGNALSMITH_LINEAR_SIZED_TYPE

template<typename V>
struct Linear;

// Everything we deal with is actually one of these
template<class BaseExpr>
struct Expression;
// A subclass with = methods
template<class BaseExpr>
struct WritableExpression;

#define EXPRESSION_NAME(nameExpr) \
	static std::string name() {\
		return nameExpr; \
	}

// Expression templates, which always hold const pointers
namespace expression {
	size_t minSize(size_t a, size_t b) {
		return std::min<size_t>(a, b);
	}

	// Expressions that just read from a pointer
	template<typename V>
	struct ReadableReal {
		EXPRESSION_NAME("ReadableReal");
		ConstRealPointer<V> pointer;

		ReadableReal(ConstRealPointer<V> pointer) : pointer(pointer) {}
		
		V get(std::ptrdiff_t i) const {
			return pointer[i];
		}
		template<class L>
		ReadableReal maybeCache(L &&, size_t) const {
			return *this;
		}
	};
	template<typename V>
	struct ReadableComplex {
		EXPRESSION_NAME("ReadableComplex");
		ConstComplexPointer<V> pointer;

		ReadableComplex(ConstComplexPointer<V> pointer) : pointer(pointer) {}

		std::complex<V> get(std::ptrdiff_t i) const {
			return pointer[i];
		}
		template<class L>
		ReadableComplex maybeCache(L &&, size_t) const {
			return *this;
		}
	};
	template<typename V>
	struct ReadableSplit {
		EXPRESSION_NAME("ReadableSplit");
		ConstSplitPointer<V> pointer;

		ReadableSplit(ConstSplitPointer<V> pointer) : pointer(pointer) {}

		std::complex<V> get(std::ptrdiff_t i) const {
			return {pointer.real[i], pointer.imag[i]};
		}
		template<class L>
		ReadableSplit maybeCache(L &&, size_t) const {
			return *this;
		}
	};
	
	template<typename V>
	struct WritableReal {
		EXPRESSION_NAME("WritableReal");
		Linear<V> &linear;
		RealPointer<V> pointer;
		size_t size;
		WritableReal(Linear<V> &linear, RealPointer<V> pointer, size_t size) : linear(linear), pointer(pointer), size(size) {}
		
		template<class Expr>
		WritableReal & operator=(Expression<Expr> expr) {
			linear.fill(pointer, expr, size);
			return *this;
		}

		V get(std::ptrdiff_t i) const {
			return pointer[i];
		}
		template<class L>
		WritableReal maybeCache(L &&, size_t) const {
			return *this;
		}
	};
	template<typename V>
	struct WritableComplex {
		EXPRESSION_NAME("WritableComplex");
		Linear<V> &linear;
		ComplexPointer<V> pointer;
		size_t size;
		WritableComplex(Linear<V> &linear, ComplexPointer<V> pointer, size_t size) : linear(linear), pointer(pointer), size(size) {}
		
		template<class Expr>
		WritableComplex & operator=(Expression<Expr> expr) {
			linear.fill(pointer, expr, size);
			return *this;
		}

		std::complex<V> get(std::ptrdiff_t i) const {
			return pointer[i];
		}
		template<class L>
		WritableComplex maybeCache(L &&, size_t) const {
			return *this;
		}
	};
	template<typename V>
	struct WritableSplit {
		EXPRESSION_NAME("WritableSplit");
		Linear<V> &linear;
		SplitPointer<V> pointer;
		size_t size;
		WritableSplit(Linear<V> &linear, SplitPointer<V> pointer, size_t size) : linear(linear), pointer(pointer), size(size) {}
		
		template<class Expr>
		WritableSplit & operator=(Expression<Expr> expr) {
			linear.fill(pointer, expr, size);
			return *this;
		}

		std::complex<V> get(std::ptrdiff_t i) const {
			return {pointer.real[i], pointer.imag[i]};
		}
		template<class L>
		WritableSplit maybeCache(L &&, size_t) const {
			return *this;
		}
	};

	// + - * / % ^ & | ~ ! = < > += -= *= /= %= ^= &= |= << >> >>= <<= == != <= >= <=>(since C++20) && || ++ -- , ->* -> ( ) [ ]
/*
#define SIGNALSMITH_AUDIO_LINEAR_UNARY_PREFIX(Name, OP) \
	template<class Right> \
	struct Name { \
		const Right right; \
		Name(const Right &right) : right(right) {} \
		auto get(std::ptrdiff_t i) const -> decltype(OP right.get(i)) { \
			return OP right.get(i); \
		} \
	}; \
	template<class Right> \
	Expression<Name<Right>> operator OP(const Expression<Right> &right) { \
		return {right}; \
	}
	SIGNALSMITH_AUDIO_LINEAR_UNARY_PREFIX(Inc, ++)
	SIGNALSMITH_AUDIO_LINEAR_UNARY_PREFIX(Dec, --)
	SIGNALSMITH_AUDIO_LINEAR_UNARY_PREFIX(Not, !)
#undef SIGNALSMITH_AUDIO_LINEAR_UNARY_PREFIX
*/

#define SIGNALSMITH_AUDIO_LINEAR_BINARY_INFIX(Name, OP) \
	template<class A, class B> \
	struct Name { \
		EXPRESSION_NAME((#Name "<") + A::name() + "," + B::name() + ">"); \
		const A a; \
		const B b; \
		Name(const A &a, const B &b) : a(a), b(b) {} \
		auto get(std::ptrdiff_t i) const -> decltype(a.get(i) OP b.get(i)) { \
			return a.get(i) OP b.get(i); \
		} \
		template<class L> \
		auto maybeCache(L &&l, size_t size) const -> Name<decltype(l.maybeCache(a, size)),decltype(l.maybeCache(b, size))> { \
			return {l.maybeCache(a, size), l.maybeCache(b, size)}; \
		}\
	}; \
	template<class A, class B> \
	Expression<Name<A, B>> operator OP(const Expression<A> &a, const Expression<B> &b) { \
		return {a, b}; \
	}
	SIGNALSMITH_AUDIO_LINEAR_BINARY_INFIX(Add, +)
	SIGNALSMITH_AUDIO_LINEAR_BINARY_INFIX(Sub, -)
	SIGNALSMITH_AUDIO_LINEAR_BINARY_INFIX(Mul, *)
	SIGNALSMITH_AUDIO_LINEAR_BINARY_INFIX(Div, /)
#undef SIGNALSMITH_AUDIO_LINEAR_BINARY_INFIX

#define SIGNALSMITH_AUDIO_LINEAR_FUNC1(Name, func) \
	template<class A> \
	struct Name { \
		EXPRESSION_NAME((#Name "<") + A::name() + ">"); \
		const A a; \
		Name(const A &a) : a(a) {} \
		auto get(std::ptrdiff_t i) const -> decltype(func(a.get(i))) { \
			return func(a.get(i)); \
		} \
		template<class L> \
		auto maybeCache(L &&l, size_t size) const -> Name<decltype(l.maybeCache(a, size))> { \
			return {l.maybeCache(a, size)}; \
		}\
	};
	template<class A>
	A fastNorm(const A &a) {
		return a*a;
	}
	template<class A>
	A fastNorm(const std::complex<A> &a) {
		A real = a.real(), imag = a.imag();
		return real*real + imag*imag;
	}
	template<class A>
	A fastAbs(const A &a) {
		return std::abs(a);
	}
	template<class A>
	A fastAbs(const std::complex<A> &a) {
		return std::hypot(a.real(), a.imag());
	}
	SIGNALSMITH_AUDIO_LINEAR_FUNC1(Abs, fastAbs)
	SIGNALSMITH_AUDIO_LINEAR_FUNC1(Norm, std::norm)
	SIGNALSMITH_AUDIO_LINEAR_FUNC1(Exp, std::exp)
	SIGNALSMITH_AUDIO_LINEAR_FUNC1(Exp2, std::exp2)
	SIGNALSMITH_AUDIO_LINEAR_FUNC1(Log, std::log)
	SIGNALSMITH_AUDIO_LINEAR_FUNC1(Log2, std::log2)
	SIGNALSMITH_AUDIO_LINEAR_FUNC1(Log10, std::log10)
	SIGNALSMITH_AUDIO_LINEAR_FUNC1(Sqrt, std::sqrt)
	SIGNALSMITH_AUDIO_LINEAR_FUNC1(Cbrt, std::cbrt)
	SIGNALSMITH_AUDIO_LINEAR_FUNC1(Floor, std::floor)
	SIGNALSMITH_AUDIO_LINEAR_FUNC1(Conj, std::conj)
	SIGNALSMITH_AUDIO_LINEAR_FUNC1(Real, std::real)
	SIGNALSMITH_AUDIO_LINEAR_FUNC1(Imag, std::imag)
	SIGNALSMITH_AUDIO_LINEAR_FUNC1(Arg, std::arg)
#undef SIGNALSMITH_AUDIO_LINEAR_FUNC1
}

template<class BaseExpr>
struct Expression : public BaseExpr {
	template<class ...Args>
	Expression(Args &&...args) : BaseExpr(std::forward<Args>(args)...) {}

	auto operator[](std::ptrdiff_t i) -> decltype(BaseExpr::get(i)) const {
		return BaseExpr::get(i);
	}

	Expression<expression::Abs<BaseExpr>> abs() const {
		return {*this};
	}
	Expression<expression::Norm<BaseExpr>> norm() const {
		return {*this};
	}
	Expression<expression::Exp<BaseExpr>> exp() const {
		return {*this};
	}
	Expression<expression::Exp2<BaseExpr>> exp2() const {
		return {*this};
	}
	Expression<expression::Log<BaseExpr>> log() const {
		return {*this};
	}
	Expression<expression::Log2<BaseExpr>> log2() const {
		return {*this};
	}
	Expression<expression::Log10<BaseExpr>> log10() const {
		return {*this};
	}
	Expression<expression::Sqrt<BaseExpr>> sqrt() const {
		return {*this};
	}
	Expression<expression::Sqrt<BaseExpr>> cbrt() const {
		return {*this};
	}
	Expression<expression::Conj<BaseExpr>> conj() const {
		return {*this};
	}
	Expression<expression::Real<BaseExpr>> real() const {
		return {*this};
	}
	Expression<expression::Imag<BaseExpr>> imag() const {
		return {*this};
	}
	Expression<expression::Arg<BaseExpr>> arg() const {
		return {*this};
	}
	Expression<expression::Floor<BaseExpr>> floor() const {
		return {*this};
	}
};
template<class BaseExpr>
struct WritableExpression : public Expression<BaseExpr> {
	using Expression<BaseExpr>::Expression;
	
	template<class Expr>
	WritableExpression & operator=(Expr &&expr) {
		BaseExpr::operator=(expr);
		return *this;
	}
};

template<typename V>
struct LinearImplBase {
	void reserve(size_t) {}

	// Wrap a pointer as an expression
	Expression<expression::ReadableReal<V>> wrap(ConstRealPointer<V> pointer) {
		return {pointer};
	}
	Expression<expression::ReadableComplex<V>> wrap(ConstComplexPointer<V> pointer) {
		return {pointer};
	}
	Expression<expression::ReadableSplit<V>> wrap(ConstSplitPointer<V> pointer) {
		return {pointer};
	}

	// TODO: instead of assignment living in the Writable***, have it in WritableExpression only, so that it still looks like a Readable*** for fill/simplify
	// When a length is supplied, make it writable
	WritableExpression<expression::WritableReal<V>> wrap(RealPointer<V> pointer, size_t size) {
		return {self(), pointer, size};
	}
	WritableExpression<expression::WritableComplex<V>> wrap(ComplexPointer<V> pointer, size_t size) {
		return {self(), pointer, size};
	}
	WritableExpression<expression::WritableSplit<V>> wrap(SplitPointer<V> pointer, size_t size) {
		return {self(), pointer, size};
	}

	WritableExpression<expression::WritableReal<V>> wrap(std::vector<V> &vector) {
		return {self(), vector.data(), vector.size()};
	}
	WritableExpression<expression::WritableComplex<V>> wrap(std::vector<std::complex<V>> &vector) {
		return {self(), vector.data(), vector.size()};
	}
	WritableExpression<expression::WritableSplit<V>> wrap(std::vector<V> &real, std::vector<V> &imag) {
		SplitPointer<V> pointer{real.data(), imag.data()};
		size_t size = std::min<size_t>(real.size(), imag.size());
		return {self(), pointer, size};
	}

	Expression<expression::WritableReal<V>> wrap(const std::vector<V> &vector) {
		return {vector.data()};
	}
	Expression<expression::WritableComplex<V>> wrap(const std::vector<std::complex<V>> &vector) {
		return {vector.data()};
	}
	Expression<expression::WritableSplit<V>> wrap(const std::vector<V> &real, const std::vector<V> &imag) {
		ConstSplitPointer<V> pointer{real.data(), imag.data()};
		return {pointer};
	}

	template<class ...Args>
	auto operator()(Args &&...args) -> decltype(wrap(std::forward<Args>(args)...)) {
		return wrap(std::forward<Args>(args)...);
	}

	// If there are fast ways to compute specific expressions, this lets us store that result in temporary space, and then return a pointer expression
	template<class Expr>
	Expr maybeCache(const Expr &expr, size_t size) {
		return expr.maybeCache(self(), size);
	}

	template<class Pointer, class Expr>
	void fill(Pointer pointer, Expr expr, size_t size) {
		auto maybeCached = self().maybeCache(expr, size);
		for (size_t i = 0; i < size; ++i) {
			pointer[i] = maybeCached.get(i);
		}
	}

	// Remove the Expression<...> layer, so the simplification template-matching works
	template<class Pointer, class Expr>
	void fill(Pointer pointer, Expression<Expr> expr, size_t size) {
		return self().fill(pointer, (Expr &)expr, size);
	};

protected:
	LinearImplBase(Linear<V> *linearThis) {
		assert((LinearImplBase *)linearThis == this);
	}
	Linear<V> & self() {
		return *(Linear<V> *)this;
	}
};

template<typename V>
struct Linear : public LinearImplBase<V> {
	Linear() : LinearImplBase<V>(this) {}

	// An example of a simplification.  This is only called when being evaluated, so it could write to temporary storage and return a Readable??? expression.
	using LinearImplBase<V>::maybeCache;
	template<class Expr>
	expression::Abs<Expr> maybeCache(const expression::Sqrt<expression::Norm<Expr>> &expr, size_t) {
		return {expr.a.a};
	}
	
	// If the simplification is a better way to write values, then override .fill() for the specific pointer/expression
	using LinearImplBase<V>::fill;
	template<class Expr>
	void fill(RealPointer<V> pointer, expression::Sqrt<expression::Norm<Expr>> expr, size_t size) {
		auto replacedExpr = maybeCache(expr, size);
		checkSimplificationWorked(replacedExpr);
		for (size_t i = 0; i < size; ++i) {
			pointer[i] = replacedExpr.get(i);
		}
	}
private:
	template<class Expr>
	void checkSimplificationWorked(Expr) {}
	// This specific pattern should've been replaced
	template<class Expr>
	void checkSimplificationWorked(expression::Sqrt<expression::Norm<Expr>> expr) = delete;
};


}}; // namespace

//#if defined(SIGNALSMITH_USE_ACCELERATE)
//#	include "./platform/linear-accelerate.h"
//#elif 0//defined(SIGNALSMITH_USE_IPP)
//#	include "./platform/linear-ipp.h"
//#elif defined(SIGNALSMITH_USE_CBLAS)
//#	include "./platform/linear-cblas.h"
//#endif

#undef SIGNALSMITH_AUDIO_LINEAR_CHUNK_SIZE
#undef SIGNALSMITH_AUDIO_LINEAR_CHUNK_FOREACH_STEP
#undef SIGNALSMITH_AUDIO_LINEAR_CHUNK_FOREACH

#endif // include guard
