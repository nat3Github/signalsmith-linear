#ifndef SIGNALSMITH_AUDIO_LINEAR_H
#define SIGNALSMITH_AUDIO_LINEAR_H

#include <cmath>
#include <complex>

namespace signalsmith { namespace linear {

/// Common linear operators (not thread-safe - use one per thread)
template<typename V, bool onlyGeneric=false>
struct Linear;

template<typename V, bool onlyGeneric=false>
struct LinearImplBase {
	using CV = std::complex<V>;

	/// Guarantees no operation will allocate if `N <= maxSize`
	/// This does nothing in the base implementation, but specialisations may working memory.
	void reserve(int /*maxSize*/) {}

#define LINEAR_CVX(ReturnType, fnName, XType, setupExpr, returnExpr, ...) \
	ReturnType fnName(const int N, const XType *x, const int xStride) { \
		if (xStride == 1) { \
			setupExpr; \
			for (int i = 0; i < N; ++i) { \
				auto xi = x[i]; \
				__VA_ARGS__; \
			} \
			returnExpr; \
		} else { \
			if (xStride < 0) x -= (N - 1)*xStride; \
			setupExpr; \
			for (int i = 0; i < N; ++i) { \
				auto xi = x[i*xStride]; \
				__VA_ARGS__; \
			} \
			returnExpr; \
		} \
	} \
	ReturnType fnName(const int N, const XType *x) { \
		return subclass().fnName(N, x, 1); \
	}

#define LINEAR_CVX_VY(ReturnType, fnName, XType, YType, setupExpr, returnExpr, ...) \
	ReturnType fnName(const int N, const XType *x, const int xStride, YType *y, const int yStride) { \
		if (xStride == 1 && yStride == 1) { \
			setupExpr; \
			for (int i = 0; i < N; ++i) { \
				auto xi = x[i]; \
				auto &yi = y[i]; \
				__VA_ARGS__; \
			} \
			returnExpr; \
		} else { \
			if (xStride < 0) x -= (N - 1)*xStride; \
			if (yStride < 0) y -= (N - 1)*yStride; \
			setupExpr; \
			for (int i = 0; i < N; ++i) { \
				auto xi = x[i*xStride]; \
				auto &yi = y[i*yStride]; \
				__VA_ARGS__; \
			} \
			returnExpr; \
		} \
	} \
	void fnName(const int N, const XType *x, YType *y) { \
		return subclass().fnName(N, x, 1, y, 1); \
	}

	LINEAR_CVX_VY(void, copy, V, V,
		/*setup*/,
		/*return*/,
		yi = xi;
	)
	LINEAR_CVX_VY(void, copy, CV, CV,
		/*setup*/,
		/*return*/,
		yi = xi;
	)

	LINEAR_CVX(V, norm2, V,
		V sum2 = 0,
		return sum2,
		sum2 += xi*xi;
	)
	LINEAR_CVX(V, norm2, CV,
		V sum2 = 0,
		return sum2,
		auto vr = xi.real(), vi = xi.imag();
		sum2 += vr*vr + vi*vi;
	)
	LINEAR_CVX_VY(void, norm2, CV, V,
		/*init*/,
		/*return*/,
		auto vr = xi.real(), vi = xi.imag();
		yi = vr*vr + vi*vi;
	);
#undef LINEAR_CVX
#undef LINEAR_CVX_VY

protected:
	/// Can only be constructed/copied when subclassed
	LinearImplBase(Linear<V, onlyGeneric> *subclassedThis) {
		// Tests for equality, and also (at compile-type) type inheritance
		if (this != subclassedThis) {
			abort();
		}
	}

private:
	Linear<V, onlyGeneric> subclass() {
		return *(Linear<V, onlyGeneric> *)this;
	}
}; // LinearImplBase

/// (Hopefully) faster version, using some acceleration library
template<typename V, bool onlyGeneric>
struct Linear : public LinearImplBase<V, onlyGeneric> {
	Linear() : LinearImplBase<V, onlyGeneric>(this) {}
};

}}; // namespace

#if 0//defined(SIGNALSMITH_USE_ACCELERATE)
#	include "./platform/linear-accelerate.h"
#elif 0//defined(SIGNALSMITH_USE_IPP)
#	include "./platform/linear-ipp.h"
#elif defined(SIGNALSMITH_USE_CBLAS)
#	include "./platform/linear-cblas.h"
#endif

#endif // include guard
