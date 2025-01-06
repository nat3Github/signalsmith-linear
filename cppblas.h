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
void copy(const int N, const V *x, V *y) {
	for (int i = 0; i < N; ++i) {
		y[i] = x[i];
	}
}

template<typename V, bool useBlas=true>
void copy(const int N, const V *x, const int xStride, V *y, const int yStride) {
	if (xStride == 1 && yStride == 1) {
		copy<V, false>(N, x, y);
	} else {
		if (xStride < 0) x -= (N - 1)*xStride;
		if (yStride < 0) y -= (N - 1)*yStride;
		for (int i = 0; i < N; ++i) {
			y[i*yStride] = x[i*xStride];
		}
	}
}

//template<typename V>
//void copy<V, true>(const int N, const V *x, V *y) {
//	copy(N, x, 1, y, 1);
//}

template<typename V, bool useBlas=true>
V norm2(const int N, const V *x) {
	V sum = 0;
	for (int i = 0; i < N; ++i) {
		auto v = x[i];
		sum += v*v;
	}
	return std::sqrt(sum);
}

template<typename V, bool useBlas=true>
V norm2(const int N, const V *x, const int xStride) {
	if (xStride == 1) {
		return norm2<V, false>(N, x);
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

//template<typename V>
//V norm2<V, true>(const int N, const V *x) {
//	return norm2(N, x, 1);
//}

template<typename V, bool useBlas=true>
V norm2(const int N, const std::complex<V> *x) {
	V sum = 0;
	for (int i = 0; i < N; ++i) {
		auto v = x[i];
		sum += std::norm(v);
	}
	return std::sqrt(sum);
}

template<typename V, bool useBlas=true>
V norm2(const int N, const std::complex<V> *x, const int xStride) {
	if (xStride == 1) {
		return norm2<V, false>(N, x);
	} else {
		const int stride = std::abs(xStride);
		V sum = 0;
		for (int i = 0; i < N; ++i) {
			auto v = x[i*stride];
			sum += std::norm(v);
		}
		return std::sqrt(sum);
	}
}

//template<typename V>
//V norm2<V, true>(const int N, const std::complex<V> *x) {
//	return norm2(N, x, 1);
//}

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
}}; // namespace
#undef BLAS_SIZE
#undef BLAS_INDEX
#endif

#endif // include guard
