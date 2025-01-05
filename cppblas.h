#ifndef SIGNALSMITH_USE_CBLAS
#	error expected SIGNALSMITH_USE_CBLAS
#endif

#ifndef SIGNALSMITH_DSP_CPPBLAS_H
#define SIGNALSMITH_DSP_CPPBLAS_H

#include <cmath>
#include <complex>

namespace signalsmith { namespace blas {

template<typename V, bool useBlas=true>
void copy(const int N, const V *x, const int xStride, V *y, const int yStride) {
	for (int i = 0; i < N; ++i) {
		y[i*yStride] = x[i*xStride];
	}
}

}}; // namespace

#ifdef SIGNALSMITH_USE_CBLAS
#define BLAS_SIZE int
extern "C" {
	void cblas_scopy(const BLAS_SIZE N, const float *X, const BLAS_SIZE incX, float *Y, const BLAS_SIZE incY);
	void cblas_dcopy(const BLAS_SIZE N, const double *X, const BLAS_SIZE incX, double *Y, const BLAS_SIZE incY);
	void cblas_ccopy(const BLAS_SIZE N, const void *X, const BLAS_SIZE incX, void *Y, const BLAS_SIZE incY);
	void cblas_zcopy(const BLAS_SIZE N, const void *X, const BLAS_SIZE incX, void *Y, const BLAS_SIZE incY);
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
}}; // namespace
#undef BLAS_SIZE
#endif

#endif // include guard
