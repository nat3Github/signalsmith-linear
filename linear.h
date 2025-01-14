#ifndef SIGNALSMITH_AUDIO_LINEAR_H
#define SIGNALSMITH_AUDIO_LINEAR_H

#include <cmath>
#include <complex>

namespace signalsmith { namespace linear {

template<typename V>
struct ConstSplitComplex {
	const V *real, *imag;
	ConstSplitComplex(const V *real, const V *imag) : real(real), imag(imag) {}
	
	std::complex<V> get(std::ptrdiff_t i) const {
		return {real[i], imag[i]};
	}
	
	ConstSplitComplex & operator +=(std::ptrdiff_t i) {
		real += i;
		imag += i;
		return *this;
	}
	ConstSplitComplex & operator -=(std::ptrdiff_t i) {
		real -= i;
		imag -= i;
		return *this;
	}
};

template<typename V>
struct SplitComplex {
	V *real, *imag;
	SplitComplex(V *real, V *imag) : real(real), imag(imag) {}

	operator ConstSplitComplex<V>() const {
		return {real, imag};
	}

	std::complex<V> get(std::ptrdiff_t i) const {
		return {real[i], imag[i]};
	}

	void set(std::ptrdiff_t i, const std::complex<V> &v) {
		real[i] = v.real();
		imag[i] = v.imag();
	}

	SplitComplex & operator +=(std::ptrdiff_t i) {
		real += i;
		imag += i;
		return *this;
	}
	SplitComplex & operator -=(std::ptrdiff_t i) {
		real -= i;
		imag -= i;
		return *this;
	}
};

/// Common linear operators (may contain temporary storage - use one per thread)
template<typename V, bool onlyGeneric=false>
struct Linear;

template<typename V, bool onlyGeneric=false>
struct LinearImplBase {
	using CV = std::complex<V>;
	using SV = SplitComplex<V>;
	using SVc = ConstSplitComplex<V>;

	/// Guarantees no operation will allocate if `N <= maxSize`
	/// This does nothing in the base implementation, but specialisations may working memory.
	void reserve(size_t /*maxSize*/) {}

#define LINEAR_VX(ReturnType, fnName, XType, setupExpr, returnExpr, ...) \
	ReturnType fnName(const int N, XType x, const int xStride) { \
		if (xStride == 1) { \
			setupExpr; \
			for (int i = 0; i < N; ++i) { \
				const int xi = i; \
				__VA_ARGS__; \
			} \
			returnExpr; \
		} else { \
			if (xStride < 0) x -= (N - 1)*xStride; \
			setupExpr; \
			for (int i = 0; i < N; ++i) { \
				const int xi = i*xStride; \
				__VA_ARGS__; \
			} \
			returnExpr; \
		} \
	} \
	ReturnType fnName(const int N, XType x) { \
		return subclass().fnName(N, x, 1); \
	}


#define LINEAR_VX_VY(ReturnType, fnName, XType, YType, setupExpr, returnExpr, ...) \
	ReturnType fnName(const int N, XType x, const int xStride, YType y, const int yStride) { \
		if (xStride == 1 && yStride == 1) { \
			setupExpr; \
			for (int i = 0; i < N; ++i) { \
				const int xi = i; \
				const int yi = i; \
				__VA_ARGS__; \
			} \
			returnExpr; \
		} else { \
			if (xStride < 0) x -= (N - 1)*xStride; \
			if (yStride < 0) y -= (N - 1)*yStride; \
			setupExpr; \
			for (int i = 0; i < N; ++i) { \
				const int xi = i*xStride; \
				const int yi = i*yStride; \
				__VA_ARGS__; \
			} \
			returnExpr; \
		} \
	} \
	void fnName(const int N, XType x, YType y) { \
		return subclass().fnName(N, x, 1, y, 1); \
	}

	LINEAR_VX_VY(void, copy, const V *, V *,
		/*setup*/,
		/*return*/,
		y[yi] = x[xi];
	)
	LINEAR_VX_VY(void, copy, const CV *, CV *,
		/*setup*/,
		/*return*/,
		y[yi] = x[xi];
	)

	LINEAR_VX(V, norm2, const V *,
		V sum2 = 0,
		return sum2,
		V v = x[xi];
		sum2 += v*v;
	)
	LINEAR_VX(V, norm2, const CV *,
		V sum2 = 0,
		return sum2,
		CV v = x[xi];
		auto vr = v.real(), vi = v.imag();
		sum2 += vr*vr + vi*vi;
	)
	LINEAR_VX(V, norm2, SVc ,
		V sum2 = 0,
		return sum2,
		auto vr = x.real[xi], vi = x.imag[xi];
		sum2 += vr*vr + vi*vi;
	)
	LINEAR_VX_VY(void, norm2, const CV *, V *,
		/*init*/,
		/*return*/,
		CV v = x[xi];
		auto vr = v.real(), vi = v.imag();
		y[yi] = vr*vr + vi*vi;
	);
	LINEAR_VX_VY(void, norm2, SVc, V *,
		/*init*/,
		/*return*/,
		auto vr = x.real[xi], vi = x.imag[xi];
		y[yi] = vr*vr + vi*vi;
	);
#undef LINEAR_VX
#undef LINEAR_VX_VY

protected:
	/// Can only be constructed/copied when subclassed
	LinearImplBase(Linear<V, onlyGeneric> *subclassedThis) {
		// Tests for equality, and also (at compile-type) type inheritance
		if (this != subclassedThis) abort();
	}

private:
	Linear<V, onlyGeneric> subclass() {
		return *(Linear<V, onlyGeneric> *)this;
	}
}; // LinearImplBase

/// The main template, which may be specialised for float/double by faster libraries
template<typename V, bool onlyGeneric>
struct Linear : public LinearImplBase<V, onlyGeneric> {
	Linear() : LinearImplBase<V, onlyGeneric>(this) {}
};

}}; // namespace

#if defined(SIGNALSMITH_USE_ACCELERATE)
#	include "./platform/linear-accelerate.h"
#elif 0//defined(SIGNALSMITH_USE_IPP)
#	include "./platform/linear-ipp.h"
#elif defined(SIGNALSMITH_USE_CBLAS)
#	include "./platform/linear-cblas.h"
#endif

#endif // include guard
