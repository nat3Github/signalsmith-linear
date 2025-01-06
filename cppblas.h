#ifndef SIGNALSMITH_USE_CBLAS
#	error expected SIGNALSMITH_USE_CBLAS
#endif

#ifndef SIGNALSMITH_DSP_CPPBLAS_H
#define SIGNALSMITH_DSP_CPPBLAS_H

#include <cmath>
#include <complex>

#include <iostream>

namespace signalsmith { namespace blas {

template<typename V, bool useBlas=true>
void copy(const int N, const V *x, const int xStride, V *y, const int yStride) {
	if (xStride == 1 && yStride == 1) {
		for (int i = 0; i < N; ++i) {
			y[i] = x[i];
		}
	} else {
		if (xStride < 0) x -= (N - 1)*xStride;
		if (yStride < 0) y -= (N - 1)*yStride;
		for (int i = 0; i < N; ++i) {
			y[i*yStride] = x[i*xStride];
		}
	}
}

template<typename V, bool useBlas=true>
void copy(const int N, const V *x, V *y) {
	copy<V, useBlas>(N, x, 1, y, 1);
}

template<typename V, bool useBlas=true>
V norm2(const int N, const V *x, const int xStride) {
	if (xStride == 1) {
		V sum = 0;
		for (int i = 0; i < N; ++i) {
			auto v = x[i];
			sum += v*v;
		}
		return std::sqrt(sum);
	} else {
		const int stride = std::abs(xStride);
		V sum = 0;
		for (int i = 0; i < N; ++i) {
			auto v = x[i*stride];
			sum += v*v;
		}
		return std::sqrt(sum);
	}
}

template<typename V, bool useBlas=true>
V norm2(const int N, const V *x) {
	return norm2<V, useBlas>(N, x, 1);
}

template<typename V, bool useBlas=true>
V norm2(const int N, const std::complex<V> *x, const int xStride) {
	if (xStride == 1) {
		V sum = 0;
		for (int i = 0; i < N; ++i) {
			auto v = x[i];
			auto vr = v.real(), vi = v.imag();
			sum += vr*vr + vi*vi;
		}
		return std::sqrt(sum);
	} else {
		const int stride = std::abs(xStride);
		V sum = 0;
		for (int i = 0; i < N; ++i) {
			auto v = x[i*stride];
			auto vr = v.real(), vi = v.imag();
			sum += vr*vr + vi*vi;
		}
		return std::sqrt(sum);
	}
}

template<typename V, bool useBlas=true>
V norm2(const int N, const std::complex<V> *x) {
	return norm2<V, useBlas>(N, x, 1);
}

}}; // namespace

#ifdef SIGNALSMITH_USE_CBLAS
#define BLAS_SIZE int
#define BLAS_INDEX size_t
extern "C" {
	void cblas_scopy(const BLAS_SIZE N, const float *X, const BLAS_SIZE incX, float *Y, const BLAS_SIZE incY);
	void cblas_dcopy(const BLAS_SIZE N, const double *X, const BLAS_SIZE incX, double *Y, const BLAS_SIZE incY);
	void cblas_ccopy(const BLAS_SIZE N, const void *X, const BLAS_SIZE incX, void *Y, const BLAS_SIZE incY);
	void cblas_zcopy(const BLAS_SIZE N, const void *X, const BLAS_SIZE incX, void *Y, const BLAS_SIZE incY);

	float cblas_snrm2(const BLAS_SIZE N, const float *X, const BLAS_SIZE incX);
	double cblas_dnrm2(const BLAS_SIZE N, const double *X, const BLAS_SIZE incX);
	float cblas_scnrm2(const BLAS_SIZE N, const void *X, const BLAS_SIZE incX);
	double cblas_dznrm2(const BLAS_SIZE N, const void *X, const BLAS_SIZE incX);

	float cblas_sdot(const BLAS_SIZE N, const float *X, const BLAS_SIZE incX, const float *Y, const BLAS_SIZE incY);
	double cblas_ddot(const BLAS_SIZE N, const double *X, const BLAS_SIZE incX, const double *Y, const BLAS_SIZE incY);
	float cblas_cdotc(const BLAS_SIZE N, const void *X, const BLAS_SIZE incX, void *Y, const BLAS_SIZE incY);
	float cblas_cdotu(const BLAS_SIZE N, const void *X, const BLAS_SIZE incX, void *Y, const BLAS_SIZE incY);
	double cblas_zdotc(const BLAS_SIZE N, const void *X, const BLAS_SIZE incX, void *Y, const BLAS_SIZE incY);
	double cblas_zdotu(const BLAS_SIZE N, const void *X, const BLAS_SIZE incX, void *Y, const BLAS_SIZE incY);
};

namespace signalsmith { namespace blas {
	template<>
	void copy<float, true>(const int N, const float *x, const int xStride, float *y, const int yStride) {
		cblas_scopy(BLAS_SIZE(N), x, BLAS_SIZE(xStride), y, BLAS_SIZE(yStride));
	}
	template<>
	void copy<double, true>(const int N, const double *x, const int xStride, double *y, const int yStride) {
		cblas_dcopy(BLAS_SIZE(N), x, BLAS_SIZE(xStride), y, BLAS_SIZE(yStride));
	}
	template<>
	void copy<std::complex<float>, true>(const int N, const std::complex<float> *x, const int xStride, std::complex<float> *y, const int yStride) {
		cblas_ccopy(BLAS_SIZE(N), x, BLAS_SIZE(xStride), y, BLAS_SIZE(yStride));
	}
	template<>
	void copy<std::complex<double>, true>(const int N, const std::complex<double> *x, const int xStride, std::complex<double> *y, const int yStride) {
		cblas_zcopy(BLAS_SIZE(N), x, BLAS_SIZE(xStride), y, BLAS_SIZE(yStride));
	}

//#ifdef SIGNALSMITH_USE_ACCELERATE
//#	ifdef __FAST_MATH__
//#		define SMALL_N_REAL(...) if (N <= 64) return __VA_ARGS__;
//#	else
//#		define SMALL_N_REAL(...) if (N <= 16) return __VA_ARGS__;
//#	endif
//#endif

#ifndef SMALL_N_REAL
#	define SMALL_N_REAL(...)
#endif

// the BLAS nrm2 functions have extra protections against overflowing to Infinity, but they're slow
#ifdef SIGNALSMITH_EXACT_BLAS_NRM2
	template<>
	float norm2<float, true>(const int N, const float *x, const int xStride) {
		return cblas_snrm2(BLAS_SIZE(N), x, BLAS_SIZE(xStride));
	}
	template<>
	double norm2<double, true>(const int N, const double *x, const int xStride) {
		return cblas_dnrm2(BLAS_SIZE(N), x, BLAS_SIZE(xStride));
	}
	template<>
	float norm2<float, true>(const int N, const std::complex<float> *x, const int xStride) {
		return cblas_scnrm2(BLAS_SIZE(N), x, BLAS_SIZE(xStride));
	}
	template<>
	double norm2<double, true>(const int N, const std::complex<double> *x, const int xStride) {
		return cblas_dznrm2(BLAS_SIZE(N), x, BLAS_SIZE(xStride));
	}
#else
	template<>
	float norm2<float, true>(const int N, const float *x, const int xStride) {
		SMALL_N_REAL(norm2<float, false>(N, x, xStride))
		return std::sqrt(cblas_sdot(BLAS_SIZE(N), x, BLAS_SIZE(xStride), x, BLAS_SIZE(xStride)));
	}
	template<>
	double norm2<double, true>(const int N, const double *x, const int xStride) {
		SMALL_N_REAL(norm2<double, false>(N, x, xStride))
		return std::sqrt(cblas_ddot(BLAS_SIZE(N), x, BLAS_SIZE(xStride), x, BLAS_SIZE(xStride)));
	}
	template<>
	float norm2<float, true>(const int N, const std::complex<float> *x, const int xStride) {
		SMALL_N_REAL(norm2<float, false>(N, x, xStride))
		if (std::abs(xStride) == 1) {
			return std::sqrt(cblas_sdot(BLAS_SIZE(2*N), (const float *)x, 1, (const float *)x, 1));
		} else {
			// Aliasing crimes, but this specific one is as kosher as it gets
			float real2 = cblas_sdot(BLAS_SIZE(N), (const float *)x, BLAS_SIZE(xStride*2), (const float *)x, BLAS_SIZE(xStride*2));
			float imag2 = cblas_sdot(BLAS_SIZE(N), (const float *)x + 1, BLAS_SIZE(xStride*2), (const float *)x + 1, BLAS_SIZE(xStride*2));
			return std::sqrt(real2 + imag2);
		}
	}
	template<>
	double norm2<double, true>(const int N, const std::complex<double> *x, const int xStride) {
		SMALL_N_REAL(norm2<double, false>(N, x, xStride))
		if (std::abs(xStride) == 1) {
			return std::sqrt(cblas_ddot(BLAS_SIZE(2*N), (const double *)x, 1, (const double *)x, 1));
		} else {
			// Aliasing crimes, but this specific one is as kosher as it gets
			double real2 = cblas_ddot(BLAS_SIZE(N), (const double *)x, BLAS_SIZE(xStride*2), (const double *)x, BLAS_SIZE(xStride*2));
			double imag2 = cblas_ddot(BLAS_SIZE(N), (const double *)x + 1, BLAS_SIZE(xStride*2), (const double *)x + 1, BLAS_SIZE(xStride*2));
			return std::sqrt(real2 + imag2);
		}
	}
#endif
}}; // namespace
#undef SMALL_N_REAL
#undef BLAS_SIZE
#undef BLAS_INDEX
#endif

#endif // include guard
