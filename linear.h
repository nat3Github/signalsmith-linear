#ifndef SIGNALSMITH_AUDIO_LINEAR_H
#define SIGNALSMITH_AUDIO_LINEAR_H

#include <cmath>
#include <complex>
#include <array>
#include <type_traits>
#include <cassert>

namespace signalsmith { namespace linear {

template<typename V>
using ConstRealPointer = const V *;
template<typename V>
using RealPointer = V *;

template<typename V>
using ConstComplexPointer = const std::complex<V> *;
template<typename V>
using ComplexPointer = std::complex<V> *;

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

template<bool=true>
struct LinearImpl;
using Linear = LinearImpl<true>;

// Everything we deal with is actually one of these
template<class BaseExpr>
struct Expression;
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
		ReadableReal maybeCache(L &, size_t) const {
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
		ReadableSplit maybeCache(L &, size_t) const {
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
		A a; \
		B b; \
		Name(const A &a, const B &b) : a(a), b(b) {} \
		auto get(std::ptrdiff_t i) const -> decltype(a.get(i) OP b.get(i)) const { \
			return a.get(i) OP b.get(i); \
		} \
		template<class L> \
		auto maybeCache(L &l, size_t size) const -> const Name<decltype(l.maybeCache(a, size)),decltype(l.maybeCache(b, size))> { \
			return {l.maybeCache(a, size), l.maybeCache(b, size)}; \
		}\
	}; \
} /*exit expression:: namespace */ \
template<class A, class B> \
const Expression<expression::Name<A, B>> operator OP(const Expression<A> &a, const Expression<B> &b) { \
	return {a, b}; \
} \
namespace expression {
	SIGNALSMITH_AUDIO_LINEAR_BINARY_INFIX(Add, +)
	SIGNALSMITH_AUDIO_LINEAR_BINARY_INFIX(Sub, -)
	SIGNALSMITH_AUDIO_LINEAR_BINARY_INFIX(Mul, *)
	SIGNALSMITH_AUDIO_LINEAR_BINARY_INFIX(Div, /)
#undef SIGNALSMITH_AUDIO_LINEAR_BINARY_INFIX

#define SIGNALSMITH_AUDIO_LINEAR_FUNC1(Name, func) \
	template<class A> \
	struct Name { \
		EXPRESSION_NAME((#Name "<") + A::name() + ">"); \
		A a; \
		Name(const A &a) : a(a) {} \
		auto get(std::ptrdiff_t i) const -> decltype(func(a.get(i))) { \
			return func(a.get(i)); \
		} \
		template<class L> \
		auto maybeCache(L &l, size_t size) const -> const Name<decltype(l.maybeCache(a, size))> { \
			return {l.maybeCache(a, size)}; \
		}\
	};
	
	template<class A>
	A fastAbs(const A &a) {
		return std::abs(a);
	}
	template<class A>
	A fastAbs(const std::complex<A> &a) {
		return std::hypot(a.real(), a.imag());
	}
	SIGNALSMITH_AUDIO_LINEAR_FUNC1(Abs, fastAbs)
	
	template<class A>
	A fastNorm(const A &a) {
		return a*a;
	}
	template<class A>
	A fastNorm(const std::complex<A> &a) {
		A real = a.real(), imag = a.imag();
		return real*real + imag*imag;
	}
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
	Expression(Args &&...args) : BaseExpr(std::forward<Args>(args)...) {
		static_assert(std::is_trivially_copyable<Expression>::value, "Expression<> must be trivially copyable");
		static_assert(std::is_trivially_copyable<BaseExpr>::value, "BaseExpr must be trivially copyable");
	}

	auto operator[](std::ptrdiff_t i) -> decltype(BaseExpr::get(i)) const {
		return BaseExpr::get(i);
	}

	const Expression<expression::Abs<BaseExpr>> abs() const {
		return {*this};
	}
	const Expression<expression::Norm<BaseExpr>> norm() const {
		return {*this};
	}
	const Expression<expression::Exp<BaseExpr>> exp() const {
		return {*this};
	}
	const Expression<expression::Exp2<BaseExpr>> exp2() const {
		return {*this};
	}
	const Expression<expression::Log<BaseExpr>> log() const {
		return {*this};
	}
	const Expression<expression::Log2<BaseExpr>> log2() const {
		return {*this};
	}
	const Expression<expression::Log10<BaseExpr>> log10() const {
		return {*this};
	}
	const Expression<expression::Sqrt<BaseExpr>> sqrt() const {
		return {*this};
	}
	const Expression<expression::Sqrt<BaseExpr>> cbrt() const {
		return {*this};
	}
	const Expression<expression::Conj<BaseExpr>> conj() const {
		return {*this};
	}
	const Expression<expression::Real<BaseExpr>> real() const {
		return {*this};
	}
	const Expression<expression::Imag<BaseExpr>> imag() const {
		return {*this};
	}
	const Expression<expression::Arg<BaseExpr>> arg() const {
		return {*this};
	}
	const Expression<expression::Floor<BaseExpr>> floor() const {
		return {*this};
	}
};
template<class BaseExpr>
struct WritableExpression : public Expression<BaseExpr> {
	using Expression<BaseExpr>::Expression;

	template<class Expr>
	WritableExpression & operator=(const Expr &expr) {
		BaseExpr::operator=(expr);
		return *this;
	}

	WritableExpression & operator=(const WritableExpression &expr) {
		BaseExpr::operator=(expr);
		return *this;
	}
};

template<bool useLinear=true>
struct LinearImplBase {
	using Linear = LinearImpl<useLinear>;

	template<class V>
	void reserve(size_t) {}

	template<typename V>
	struct WritableReal {
		EXPRESSION_NAME("WritableReal");
		Linear &linear;
		RealPointer<V> pointer;
		size_t size;
		WritableReal(Linear &linear, RealPointer<V> pointer, size_t size) : linear(linear), pointer(pointer), size(size) {
			static_assert(std::is_trivially_copyable<WritableReal>::value, "must be trivially copyable");
		}
		
		operator expression::ReadableReal<V>() const {
			return {pointer};
		}
		
		template<class Expr>
		WritableReal & operator=(Expression<Expr> expr) {
			linear.fill(pointer, expr, size);
			return *this;
		}

		V get(std::ptrdiff_t i) const {
			return pointer[i];
		}
		template<class L>
		expression::ReadableReal<V> maybeCache(L &&, size_t) const {
			return {pointer};
		}
	};
	template<typename V>
	struct WritableComplex {
		EXPRESSION_NAME("WritableComplex");
		Linear &linear;
		ComplexPointer<V> pointer;
		size_t size;
		WritableComplex(Linear &linear, ComplexPointer<V> pointer, size_t size) : linear(linear), pointer(pointer), size(size) {}

		operator expression::ReadableComplex<V>() const {
			return {pointer};
		}

		template<class Expr>
		WritableComplex & operator=(Expression<Expr> expr) {
			linear.fill(pointer, expr, size);
			return *this;
		}

		std::complex<V> get(std::ptrdiff_t i) const {
			return pointer[i];
		}
		template<class L>
		expression::ReadableComplex<V> maybeCache(L &, size_t) const {
			return {pointer};
		}
	};
	template<typename V>
	struct WritableSplit {
		EXPRESSION_NAME("WritableSplit");
		Linear &linear;
		SplitPointer<V> pointer;
		size_t size;
		WritableSplit(Linear &linear, SplitPointer<V> pointer, size_t size) : linear(linear), pointer(pointer), size(size) {}

		operator expression::ReadableSplit<V>() const {
			return {pointer};
		}

		template<class Expr>
		WritableSplit & operator=(Expression<Expr> expr) {
			linear.fill(pointer, expr, size);
			return *this;
		}

		std::complex<V> get(std::ptrdiff_t i) const {
			return {pointer.real[i], pointer.imag[i]};
		}
		template<class L>
		expression::ReadableSplit<V> maybeCache(L &&, size_t) const {
			return {pointer};
		}
	};
	
	// Wrap a pointer as an expression
	template<typename V>
	Expression<expression::ReadableReal<V>> wrap(ConstRealPointer<V> pointer) {
		return {pointer};
	}
	template<typename V>
	Expression<expression::ReadableComplex<V>> wrap(ConstComplexPointer<V> pointer) {
		return {pointer};
	}
	template<typename V>
	Expression<expression::ReadableSplit<V>> wrap(ConstSplitPointer<V> pointer) {
		return {pointer};
	}

	// TODO: instead of assignment living in the Writable***, have it in WritableExpression only, so that it still looks like a Readable*** for fill/simplify
	// When a length is supplied, make it writable
	template<typename V>
	WritableExpression<WritableReal<V>> wrap(RealPointer<V> pointer, size_t size) {
		return {self(), pointer, size};
	}
	template<typename V>
	WritableExpression<WritableComplex<V>> wrap(ComplexPointer<V> pointer, size_t size) {
		return {self(), pointer, size};
	}
	template<typename V>
	WritableExpression<WritableSplit<V>> wrap(SplitPointer<V> pointer, size_t size) {
		return {self(), pointer, size};
	}

	template<typename V>
	WritableExpression<WritableReal<V>> wrap(std::vector<V> &vector) {
		return {self(), vector.data(), vector.size()};
	}
	template<typename V>
	WritableExpression<WritableComplex<V>> wrap(std::vector<std::complex<V>> &vector) {
		return {self(), vector.data(), vector.size()};
	}
	template<typename V>
	WritableExpression<WritableSplit<V>> wrap(std::vector<V> &real, std::vector<V> &imag) {
		SplitPointer<V> pointer{real.data(), imag.data()};
		size_t size = std::min<size_t>(real.size(), imag.size());
		return {self(), pointer, size};
	}

	template<typename V>
	Expression<expression::ReadableReal<V>> wrap(const std::vector<V> &vector) {
		return {vector.data()};
	}
	template<typename V>
	Expression<expression::ReadableComplex<V>> wrap(const std::vector<std::complex<V>> &vector) {
		return {vector.data()};
	}
	template<typename V>
	Expression<expression::ReadableSplit<V>> wrap(const std::vector<V> &real, const std::vector<V> &imag) {
		ConstSplitPointer<V> pointer{real.data(), imag.data()};
		return {pointer};
	}

	template<class ...Args>
	auto operator()(Args &&...args) -> decltype(wrap(std::forward<Args>(args)...)) {
		return wrap(std::forward<Args>(args)...);
	}

	// If there are fast ways to compute specific expressions, this lets us store that result in temporary space, and then return a pointer expression
	template<class Expr>
	auto maybeCache(const Expr &expr, size_t size) -> decltype(expr.maybeCache(*this, size)) {
		return expr.maybeCache(*this, size);
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
	LinearImplBase(Linear *linearThis) {
		assert((LinearImplBase *)linearThis == this);
	}

	Linear & self() {
		return *(Linear *)this;
	}

	// Nothing in the fallback uses this, but specialisations might use/return internal temporary storage (possibly returning a struct which tracks lifetime/etc.) and then use .fill()
	template<class Expr>
	void toPointer(const Expr &expr, size_t size) = delete;
	
	template<typename V>
	ConstRealPointer<V> toPointer(const expression::ReadableReal<V> &expr, size_t) {
		return expr.pointer;
	}
	template<typename V>
	ConstComplexPointer<V> toPointer(const expression::ReadableComplex<V> &expr, size_t) {
		return expr.pointer;
	}
	template<typename V>
	ConstSplitPointer<V> toPointer(const expression::ReadableSplit<V> &expr, size_t) {
		return expr.pointer;
	}

	template<typename V>
	RealPointer<V> toPointer(const WritableReal<V> &expr, size_t) {
		return expr.pointer;
	}
	template<typename V>
	ComplexPointer<V> toPointer(const WritableComplex<V> &expr, size_t) {
		return expr.pointer;
	}
	template<typename V>
	SplitPointer<V> toPointer(const WritableSplit<V> &expr, size_t) {
		return expr.pointer;
	}
};

// SFINAE template for checking that an expression naturally returns a particular item type
template<class InputExpr, typename Item, class OutputExpr>
using ItemType = typename std::enable_if<
	std::is_same<
		typename std::decay<decltype(std::declval<InputExpr>().get(0))>::type,
		Item
	>::value,
	OutputExpr
>::type;

// Fallback implementation - this should be specialised (with useLinear=true) with faster methods where available
template<bool useLinear>
struct LinearImpl : public LinearImplBase<useLinear> {
	LinearImpl() : LinearImplBase<useLinear>(this) {}

	// An example of a simplification.  This is only called when being evaluated, so it could write to temporary storage and return a Readable??? expression.  Commenting this example out should make the `.fill()` below fail to compile.
	using LinearImplBase<useLinear>::maybeCache;

	template<class Expr>
	ItemType<Expr, std::complex<float>, expression::Abs<Expr>> maybeCache(const expression::Sqrt<expression::Norm<Expr>> &expr, size_t) {
		return {expr.a.a};
	}

	template<class Expr>
	ItemType<Expr, std::complex<double>, expression::Abs<Expr>> maybeCache(const expression::Sqrt<expression::Norm<Expr>> &expr, size_t) {
		return {expr.a.a};
	}

	// If the simplification is a better way to write values, then override .fill() for the specific pointer/expression.  Since `.maybeCache()` should only be called from `.fill()`, this still works even if `.maybeCache()` replaces the same pattern.
	using LinearImplBase<useLinear>::fill;
	template<typename V, class Expr>
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

/// Helper class for temporary storage, reserved up-front and tracked with values on the stack.
template<typename V>
struct Temporary {
	// This is called if we don't have enough reserved space and end up allocating
	std::function<void(size_t)> allocationWarning;
	
	void reserve(size_t size) {
		assert(start == buffer);
		
		if (buffer) delete[] buffer;
		start = buffer = new V[size];
		end = buffer + size;
	}
	
	// A chunk of temporary storage, valid as long as it's in scope on the stack.
	struct StackScoped {
		V *pointer;
		
		StackScoped(Temporary &temporary, size_t size) : pointer(temporary.start), temporary(temporary) {
			temporary.start += size;
			if (temporary.start > temporary.end) {
				// OK, actually we ran out of temporary space, so allocate
				pointer = new V[size];
				// but we're not happy about it. >:(
				temporary.allocationWarning(temporary.start - temporary.buffer);
			}
		}
		StackScoped(const StackScoped &other) = delete; // no copy/move/etc.
		~StackScoped() {
			if (pointer >= temporary.buffer && pointer < temporary.end) {
				assert(pointer <= temporary.start); // checks (although doesn't guarantee) that the storage wasn't re-allocated under our feet somehow
				temporary.start = pointer;
			} else {
				// We ran out of space, so it was allocated just for this
				delete[] pointer;
			}
		}
		
		operator V*() const {
			return pointer;
		}
	private:
		Temporary &temporary;
	};

	StackScoped scoped(size_t size) {
		return {*this, size};
	}
private:
	V *start = nullptr, *end = nullptr;
	V *buffer = nullptr;
};

}}; // namespace

#if defined(SIGNALSMITH_USE_ACCELERATE)
#	include "./platform/linear-accelerate.h"
#elif 0//defined(SIGNALSMITH_USE_IPP)
#	include "./platform/linear-ipp.h"
#elif 0//defined(SIGNALSMITH_USE_CBLAS)
#	include "./platform/linear-cblas.h"
#endif

#undef SIGNALSMITH_AUDIO_LINEAR_CHUNK_SIZE
#undef SIGNALSMITH_AUDIO_LINEAR_CHUNK_FOREACH_STEP
#undef SIGNALSMITH_AUDIO_LINEAR_CHUNK_FOREACH

#endif // include guard
