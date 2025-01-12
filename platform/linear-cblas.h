#if 1 // these make sense for Accelerate's CBLAS, but probably aren't universal
#	ifdef __FAST_MATH__
#		define SMALL_N_STRIDE(complexity, strideMaybe1, ...) if (N <= ((strideMaybe1) == 1 ? int(512/complexity) : int(64/complexity))) return __VA_ARGS__;
#	else
#		define SMALL_N_STRIDE(complexity, strideMaybe1, ...) if (N <= ((strideMaybe1) == 1 ? int(128/complexity) : int(64/complexity))) return __VA_ARGS__;
#	endif
#else
#	define SMALL_N_STRIDE(...)
#endif

#define BLAS_SIZE int
#define BLAS_INDEX size_t
extern "C" {
	void cblas_scopy(const BLAS_SIZE N, const float *X, const BLAS_SIZE incX, float *Y, const BLAS_SIZE incY);
	void cblas_dcopy(const BLAS_SIZE N, const double *X, const BLAS_SIZE incX, double *Y, const BLAS_SIZE incY);
	void cblas_ccopy(const BLAS_SIZE N, const void *X, const BLAS_SIZE incX, void *Y, const BLAS_SIZE incY);
	void cblas_zcopy(const BLAS_SIZE N, const void *X, const BLAS_SIZE incX, void *Y, const BLAS_SIZE incY);

	float cblas_sdot(const BLAS_SIZE N, const float *X, const BLAS_SIZE incX, const float *Y, const BLAS_SIZE incY);
	double cblas_ddot(const BLAS_SIZE N, const double *X, const BLAS_SIZE incX, const double *Y, const BLAS_SIZE incY);
	float cblas_cdotc(const BLAS_SIZE N, const void *X, const BLAS_SIZE incX, void *Y, const BLAS_SIZE incY);
	float cblas_cdotu(const BLAS_SIZE N, const void *X, const BLAS_SIZE incX, void *Y, const BLAS_SIZE incY);
	double cblas_zdotc(const BLAS_SIZE N, const void *X, const BLAS_SIZE incX, void *Y, const BLAS_SIZE incY);
	double cblas_zdotu(const BLAS_SIZE N, const void *X, const BLAS_SIZE incX, void *Y, const BLAS_SIZE incY);
};

namespace signalsmith { namespace linear {

template<>
struct Linear<float> : public LinearImplBase<float> {
	using Base = LinearImplBase<float>;

	Linear() : Base(this) {}

	using Base::copy;
	void copy(const int N, const float *x, const int xStride, float *y, const int yStride) {
		SMALL_N_STRIDE(1, xStride|yStride, Base::copy(N, x, xStride, y, yStride))
		cblas_scopy(BLAS_SIZE(N), x, BLAS_SIZE(xStride), y, BLAS_SIZE(yStride));
	}
	void copy(const int N, const std::complex<float> *x, const int xStride, std::complex<float> *y, const int yStride) {
		SMALL_N_STRIDE(2, xStride|yStride, Base::copy(N, x, xStride, y, yStride))
		cblas_ccopy(BLAS_SIZE(N), x, BLAS_SIZE(xStride), y, BLAS_SIZE(yStride));
	}

	using Base::norm2;
	float norm2(const int N, const float *x, const int xStride) {
		SMALL_N_STRIDE(2, xStride, Base::norm2(N, x, xStride))
		return cblas_sdot(BLAS_SIZE(N), x, BLAS_SIZE(xStride), x, BLAS_SIZE(xStride));
	}
	float norm2(const int N, const std::complex<float> *x, const int xStride) {
		SMALL_N_STRIDE(4, xStride, Base::norm2(N, x, xStride))
		if (std::abs(xStride) == 1) {
			return cblas_sdot(BLAS_SIZE(2*N), (const float *)x, 1, (const float *)x, 1);
		} else {
			// Aliasing crimes, but this specific one is as kosher as it gets
			float real2 = cblas_sdot(BLAS_SIZE(N), (const float *)x, BLAS_SIZE(xStride*2), (const float *)x, BLAS_SIZE(xStride*2));
			float imag2 = cblas_sdot(BLAS_SIZE(N), (const float *)x + 1, BLAS_SIZE(xStride*2), (const float *)x + 1, BLAS_SIZE(xStride*2));
			return real2 + imag2;
		}
	}
};

template<>
struct Linear<double> : public LinearImplBase<double> {
	using Base = LinearImplBase<double>;

	Linear() : Base(this) {}

	using Base::copy;
	void copy(const int N, const double *x, const int xStride, double *y, const int yStride) {
		SMALL_N_STRIDE(1, xStride|yStride, Base::copy(N, x, xStride, y, yStride))
		cblas_dcopy(BLAS_SIZE(N), x, BLAS_SIZE(xStride), y, BLAS_SIZE(yStride));
	}
	void copy(const int N, const std::complex<double> *x, const int xStride, std::complex<double> *y, const int yStride) {
		SMALL_N_STRIDE(2, xStride|yStride, Base::copy(N, x, xStride, y, yStride))
		cblas_zcopy(BLAS_SIZE(N), x, BLAS_SIZE(xStride), y, BLAS_SIZE(yStride));
	}

	using Base::norm2;
	double norm2(const int N, const double *x, const int xStride) {
		SMALL_N_STRIDE(2, xStride, Base::norm2(N, x, xStride))
		return cblas_ddot(BLAS_SIZE(N), x, BLAS_SIZE(xStride), x, BLAS_SIZE(xStride));
	}
	double norm2(const int N, const std::complex<double> *x, const int xStride) {
		SMALL_N_STRIDE(4, xStride, Base::norm2(N, x, xStride))
		if (std::abs(xStride) == 1) {
			return cblas_ddot(BLAS_SIZE(2*N), (const double *)x, 1, (const double *)x, 1);
		} else {
			// Aliasing crimes, but this specific one is as kosher as it gets
			double real2 = cblas_ddot(BLAS_SIZE(N), (const double *)x, BLAS_SIZE(xStride*2), (const double *)x, BLAS_SIZE(xStride*2));
			double imag2 = cblas_ddot(BLAS_SIZE(N), (const double *)x + 1, BLAS_SIZE(xStride*2), (const double *)x + 1, BLAS_SIZE(xStride*2));
			return real2 + imag2;
		}
	}
};

}}; // namespace

#undef SMALL_N_STRIDE
#undef BLAS_SIZE
#undef BLAS_INDEX
