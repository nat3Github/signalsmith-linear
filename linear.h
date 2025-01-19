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

// Expression templates, which always hold const pointers
namespace expression {
	size_t minSize(size_t a, size_t b) {
		return std::min<size_t>(a, b);
	}

	// Expressions that just read from a pointer
	template<typename V>
	struct ReadableReal {
		ConstRealPointer<V> pointer;

		ReadableReal(ConstRealPointer<V> pointer) : pointer(pointer) {}
		
		V get(std::ptrdiff_t i) const {
			return pointer[i];
		}
		RealChunk<V> getChunk(std::ptrdiff_t i) const {
			RealChunk<V> result;
			auto *offsetPointer = pointer + i;
			SIGNALSMITH_AUDIO_LINEAR_CHUNK_FOREACH(o, result[o] = offsetPointer[o]);
			return result;
		}
	};
	template<typename V>
	struct ReadableComplex {
		ConstComplexPointer<V> pointer;

		ReadableComplex(ConstComplexPointer<V> pointer) : pointer(pointer) {}

		std::complex<V> get(std::ptrdiff_t i) const {
			return pointer[i];
		}
		ComplexChunk<V> getChunk(std::ptrdiff_t i) const {
			ComplexChunk<V> result;
			auto *offsetPointer = pointer + i;
			SIGNALSMITH_AUDIO_LINEAR_CHUNK_FOREACH(o, result[o] = offsetPointer[o]);
			return result;
		}
	};
	template<typename V>
	struct ReadableSplit {
		ConstSplitPointer<V> pointer;

		ReadableSplit(ConstSplitPointer<V> pointer) : pointer(pointer) {}

		std::complex<V> get(std::ptrdiff_t i) const {
			return {pointer.real[i], pointer.imag[i]};
		}
		SplitChunk<V> getChunk(std::ptrdiff_t i) const {
			SplitChunk<V> result;
			auto *offsetReal = pointer.real + i;
			SIGNALSMITH_AUDIO_LINEAR_CHUNK_FOREACH(o, result.real[o] = offsetReal[o]);
			auto *offsetImag = pointer.imag + i;
			SIGNALSMITH_AUDIO_LINEAR_CHUNK_FOREACH(o, result.imag[o] = offsetImag[o]);
			return result;
		}
	};
	
	template<typename V>
	struct WritableReal {
		Linear<V> &linear;
		RealPointer<V> pointer;
		size_t size;
		WritableReal(Linear<V> &linear, RealPointer<V> pointer, size_t size) : linear(linear), pointer(pointer), size(size) {}
		
		template<class Expr>
		WritableReal & operator=(Expr expr) {
			auto maybeCached = linear.cachedExpr(expr);
			for (size_t i = 0; i < size; ++i) {
				pointer[i] = maybeCached.get(i);
			}
			return *this;
		}

		V get(std::ptrdiff_t i) const {
			return pointer[i];
		}
		RealChunk<V> getChunk(std::ptrdiff_t i) const {
			RealChunk<V> result;
			auto *offsetPointer = pointer + i;
			SIGNALSMITH_AUDIO_LINEAR_CHUNK_FOREACH(o, result[o] = offsetPointer[o]);
			return result;
		}
	};
	template<typename V>
	struct WritableComplex {
		Linear<V> &linear;
		ComplexPointer<V> pointer;
		size_t size;
		WritableComplex(Linear<V> &linear, ComplexPointer<V> pointer, size_t size) : linear(linear), pointer(pointer), size(size) {}
		
		template<class Expr>
		WritableComplex & operator=(Expr expr) {
			auto maybeCached = linear.cachedExpr(expr);
			for (size_t i = 0; i < size; ++i) {
				pointer[i] = maybeCached.get(i);
			}
			return *this;
		}

		std::complex<V> get(std::ptrdiff_t i) const {
			return pointer[i];
		}
		ComplexChunk<V> getChunk(std::ptrdiff_t i) const {
			ComplexChunk<V> result;
			auto *offsetPointer = pointer + i;
			SIGNALSMITH_AUDIO_LINEAR_CHUNK_FOREACH(o, result[o] = offsetPointer[o]);
			return result;
		}
	};
	template<typename V>
	struct WritableSplit {
		Linear<V> &linear;
		SplitPointer<V> pointer;
		size_t size;
		WritableSplit(Linear<V> &linear, SplitPointer<V> pointer, size_t size) : linear(linear), pointer(pointer), size(size) {}
		
		template<class Expr>
		WritableSplit & operator=(Expr expr) {
			auto maybeCached = linear.cachedExpr(expr);
			for (size_t i = 0; i < size; ++i) {
				pointer[i] = maybeCached.get(i);
			}
			return *this;
		}

		std::complex<V> get(std::ptrdiff_t i) const {
			return {pointer.real[i], pointer.imag[i]};
		}
		SplitChunk<V> getChunk(std::ptrdiff_t i) const {
			SplitChunk<V> result;
			auto *offsetReal = pointer.real + i;
			SIGNALSMITH_AUDIO_LINEAR_CHUNK_FOREACH(o, result.real[o] = offsetReal[o]);
			auto *offsetImag = pointer.imag + i;
			SIGNALSMITH_AUDIO_LINEAR_CHUNK_FOREACH(o, result.imag[o] = offsetImag[o]);
			return result;
		}
	};
}

template<typename V>
struct Linear {
	void reserve(size_t) {}

	// This is where we'd add in convenience methods
	template<class BaseExpr>
	struct Expression : public BaseExpr {
		template<class ...Args>
		Expression(Args &&...args) : BaseExpr(std::forward<Args>(args)...) {}
		
		auto operator[](std::ptrdiff_t i) -> decltype(BaseExpr::get(i)) const {
			return BaseExpr::get(i);
		}
	};
	template<class BaseExpr>
	struct WritableExpression : public Expression<BaseExpr> {
		template<class ...Args>
		WritableExpression(Args &&...args) : Expression<BaseExpr>(std::forward<Args>(args)...) {}
		
		template<class Expr>
		WritableExpression & operator=(Expr &&expr) {
			BaseExpr::operator=(expr);
			return *this;
		}
	};

	// Wrap a pointer as an expression
	Expression<expression::ReadableReal<V>> wrap(RealPointer<V> pointer) {
		return {pointer};
	}
	Expression<expression::ReadableComplex<V>> wrap(ComplexPointer<V> pointer) {
		return {pointer};
	}
	Expression<expression::ReadableSplit<V>> wrap(SplitPointer<V> pointer) {
		return {pointer};
	}

	// When a length is supplied, make it writable
	WritableExpression<expression::WritableReal<V>> wrap(RealPointer<V> pointer, size_t size) {
		return {*this, pointer, size};
	}
	WritableExpression<expression::WritableComplex<V>> wrap(ComplexPointer<V> pointer, size_t size) {
		return {*this, pointer, size};
	}
	WritableExpression<expression::WritableSplit<V>> wrap(SplitPointer<V> pointer, size_t size) {
		return {*this, pointer, size};
	}

	// If there are fast ways to compute specific expressions, this lets us store that result in temporary space, and then return a pointer expression
	template<class Expr>
	Expr cachedExpr(const Expr &expr) {
		return expr;
	}
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
