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

#define LINEAR_VA(ReturnType, fnName, AType, setupExpr, returnExpr, ...) \
	ReturnType fnName(const int N, AType a, const int aStride) { \
		if (aStride == 1) { \
			setupExpr; \
			for (int i = 0; i < N; ++i) { \
				const int ai = i; \
				__VA_ARGS__; \
			} \
			returnExpr; \
		} else { \
			if (aStride < 0) a -= (N - 1)*aStride; \
			setupExpr; \
			for (int i = 0; i < N; ++i) { \
				const int ai = i*aStride; \
				__VA_ARGS__; \
			} \
			returnExpr; \
		} \
	} \
	ReturnType fnName(const int N, AType a) { \
		return subclass().fnName(N, a, 1); \
	}


#define LINEAR_VA_VB(ReturnType, fnName, AType, BType, setupExpr, returnExpr, ...) \
	ReturnType fnName(const int N, AType a, const int aStride, BType b, const int bStride) { \
		if (aStride == 1 && bStride == 1) { \
			setupExpr; \
			for (int i = 0; i < N; ++i) { \
				const int ai = i; \
				const int bi = i; \
				__VA_ARGS__; \
			} \
			returnExpr; \
		} else { \
			if (aStride < 0) a -= (N - 1)*aStride; \
			if (bStride < 0) b -= (N - 1)*bStride; \
			setupExpr; \
			for (int i = 0; i < N; ++i) { \
				const int ai = i*aStride; \
				const int bi = i*bStride; \
				__VA_ARGS__; \
			} \
			returnExpr; \
		} \
	} \
	void fnName(const int N, AType a, BType b) { \
		return subclass().fnName(N, a, 1, b, 1); \
	}

#define LINEAR_VA_VB_VC(ReturnType, fnName, AType, BType, CType, setupExpr, returnExpr, ...) \
	ReturnType fnName(const int N, AType a, const int aStride, BType b, const int bStride, CType c, const int cStride) { \
		if (aStride == 1 && bStride == 1 && cStride == 1) { \
			setupExpr; \
			for (int i = 0; i < N; ++i) { \
				const int ai = i; \
				const int bi = i; \
				const int ci = i; \
				__VA_ARGS__; \
			} \
			returnExpr; \
		} else { \
			if (aStride < 0) a -= (N - 1)*aStride; \
			if (bStride < 0) b -= (N - 1)*bStride; \
			if (cStride < 0) c -= (N - 1)*cStride; \
			setupExpr; \
			for (int i = 0; i < N; ++i) { \
				const int ai = i*aStride; \
				const int bi = i*bStride; \
				const int ci = i*cStride; \
				__VA_ARGS__; \
			} \
			returnExpr; \
		} \
	} \
	void fnName(const int N, AType a, BType b, CType c) { \
		return subclass().fnName(N, a, 1, b, 1, c, 1); \
	}

	LINEAR_VA_VB(void, copy, const V *, V *,
		/*setup*/,
		/*return*/,
		b[bi] = a[ai];
	)
	LINEAR_VA_VB(void, copy, const CV *, CV *,
		/*setup*/,
		/*return*/,
		b[bi] = a[ai];
	)

	LINEAR_VA(V, norm2, const V *,
		V sum2 = 0,
		return sum2,
		V v = a[ai];
		sum2 += v*v;
	)
	LINEAR_VA(V, norm2, const CV *,
		V sum2 = 0,
		return sum2,
		CV v = a[ai];
		auto vr = v.real(), vi = v.imag();
		sum2 += vr*vr + vi*vi;
	)
	LINEAR_VA(V, norm2, SVc ,
		V sum2 = 0,
		return sum2,
		auto vr = a.real[ai], vi = a.imag[ai];
		sum2 += vr*vr + vi*vi;
	)
	LINEAR_VA_VB(void, norm2, const CV *, V *,
		/*init*/,
		/*return*/,
		CV v = a[ai];
		auto vr = v.real(), vi = v.imag();
		b[bi] = vr*vr + vi*vi;
	);
	LINEAR_VA_VB(void, norm2, SVc, V *,
		/*init*/,
		/*return*/,
		auto vr = a.real[ai], vi = a.imag[ai];
		b[bi] = vr*vr + vi*vi;
	);

	LINEAR_VA_VB_VC(V, mul, const V *, const V *, V *,
		/*init*/,
		/*return*/,
		c[ci] = a[ai]*b[bi];
	)
	LINEAR_VA_VB_VC(V, mul, const CV *, const CV *, CV *,
		/*init*/,
		/*return*/,
		auto va = a[ai], vb = b[bi];
		auto var = va.real(), vai = va.imag(), vbr = vb.real(), vbi = vb.imag();
		c[ci] = {var*vbr - vai*vbi, var*vbi + vbr*vai};
	)
	LINEAR_VA_VB_VC(V, mul, SVc, SVc, SV,
		/*init*/,
		/*return*/,
		auto var = a.real[ai], vai = a.imag[ai], vbr = b.real[bi], vbi = b.imag[bi];
		c.real[ci] = var*vbr - vai*vbi;
		c.imag[ci] = var*vbi + vbr*vai;
	)
	LINEAR_VA_VB_VC(V, mul, const CV *, const V *, CV *,
		/*init*/,
		/*return*/,
		c[ci] = a[ai]*b[bi];
	)
	LINEAR_VA_VB_VC(V, mul, SVc, const V *, SV,
		/*init*/,
		/*return*/,
		auto vb = b[bi];
		c.real[ci] = a.real[ai]*vb;
		c.imag[ci] = a.imag[ai]*vb;
	)
	LINEAR_VA_VB_VC(V, mul, const V *, const CV *, CV *,
		/*init*/,
		/*return*/,
		c[ci] = a[ai]*b[bi];
	)
	LINEAR_VA_VB_VC(V, mul, const V *, SVc, SV,
		/*init*/,
		/*return*/,
		auto va = a[ai];
		c.real[ci] = va*b.real[bi];
		c.imag[ci] = va*b.imag[bi];
	)

	LINEAR_VA_VB(V, mul, const V *, V *,
		/*init*/,
		/*return*/,
		b[bi] = a[ai]*b[bi];
	)
	LINEAR_VA_VB(V, mul, const CV *, CV *,
		/*init*/,
		/*return*/,
		auto va = a[ai], vb = b[bi];
		auto var = va.real(), vai = va.imag(), vbr = vb.real(), vbi = vb.imag();
		b[bi] = {var*vbr - vai*vbi, var*vbi + vbr*vai};
	)
	LINEAR_VA_VB(V, mul, SVc, SV,
		/*init*/,
		/*return*/,
		auto var = a.real[ai], vai = a.imag[ai], vbr = b.real[bi], vbi = b.imag[bi];
		b.real[bi] = var*vbr - vai*vbi;
		b.imag[bi] = var*vbi + vbr*vai;
	)
	LINEAR_VA_VB(V, mul, const V *, CV *,
		/*init*/,
		/*return*/,
		b[bi] = a[ai]*b[bi];
	)
	LINEAR_VA_VB(V, mul, const V *, SV,
		/*init*/,
		/*return*/,
		auto va = a[ai];
		b.real[bi] *= va;
		b.imag[bi] *= va;
	)
#undef LINEAR_VA
#undef LINEAR_VA_VB
#undef LINEAR_VA_VB_VC

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
