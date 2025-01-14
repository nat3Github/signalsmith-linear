#define ACCELERATE_NEW_LAPACK

#include <Accelerate/Accelerate.h>

#ifndef SIGNALSMITH_USE_CBLAS
#	warning Linked to Accelerate, but not its CBLAS interface (`-lcblas`)
#endif

#ifdef SIGNALSMITH_USE_CBLAS
#	define CBLAS_INT __LAPACK_int
#	define CBLAS_INDEX size_t
#	define CBLAS_S __LAPACK_float_complex
#	define CBLAS_Z __LAPACK_double_complex
extern "C" {
	void cblas_scopy(const CBLAS_INT N, const float *X, const CBLAS_INT incX, float *Y, const CBLAS_INT incY);
	void cblas_dcopy(const CBLAS_INT N, const double *X, const CBLAS_INT incX, double *Y, const CBLAS_INT incY);
	void cblas_ccopy(const CBLAS_INT N, const CBLAS_S *X, const CBLAS_INT incX, CBLAS_S *Y, const CBLAS_INT incY);
	void cblas_zcopy(const CBLAS_INT N, const CBLAS_Z *X, const CBLAS_INT incX, CBLAS_Z *Y, const CBLAS_INT incY);
};
#endif

#if 1 // roughly tuned against Accelerate's CBLAS
#	ifdef __FAST_MATH__
#		define SMALL_N_STRIDE(complexity, strideMaybe1, ...) if (N <= ((strideMaybe1) == 1 ? int(512/complexity) : int(64/complexity))) return __VA_ARGS__;
#	else
#		define SMALL_N_STRIDE(complexity, strideMaybe1, ...) if (N <= ((strideMaybe1) == 1 ? int(128/complexity) : int(64/complexity))) return __VA_ARGS__;
#	endif
#else
#	define SMALL_N_STRIDE(...)
#endif

namespace signalsmith { namespace linear {

template<>
struct Linear<float> : public LinearImplBase<float> {
	using Base = LinearImplBase<float>;

	Linear() : Base(this) {}

	using Base::copy;
#ifdef SIGNALSMITH_USE_CBLAS
	void copy(const int N, const float *x, const int xStride, float *y, const int yStride) {
		SMALL_N_STRIDE(1, xStride|yStride, Base::copy(N, x, xStride, y, yStride))
		cblas_scopy(CBLAS_INT(N), x, CBLAS_INT(xStride), y, CBLAS_INT(yStride));
	}
	void copy(const int N, const std::complex<float> *x, const int xStride, std::complex<float> *y, const int yStride) {
		SMALL_N_STRIDE(2, xStride|yStride, Base::copy(N, x, xStride, y, yStride))
		cblas_ccopy(CBLAS_INT(N), x, CBLAS_INT(xStride), y, CBLAS_INT(yStride));
	}
#endif

	using Base::norm2;
	float norm2(const int N, const float *x, const int xStride) {
//		SMALL_N_STRIDE(2, xStride, Base::norm2(N, x, xStride))
		float sum2 = 0;
		vDSP_svesq(x, std::abs(xStride), &sum2, N);
		return sum2;
	}
	float norm2(const int N, const std::complex<float> *x, const int xStride) {
//		SMALL_N_STRIDE(4, xStride, Base::norm2(N, x, xStride))
		if (std::abs(xStride) == 1) {
			float sum2 = 0;
			// Aliasing crimes, but this specific one is as kosher as it gets
			vDSP_svesq((const float *)x, std::abs(xStride), &sum2, N*2);
			return sum2;
		} else {
			// criiiimes
			float real2 = 0, imag2 = 0;
			vDSP_svesq((const float *)x, std::abs(xStride)*2, &real2, N);
			vDSP_svesq((const float *)x + 1, std::abs(xStride)*2, &imag2, N);
			return real2 + imag2;
		}
	}

	float norm2(const int N, ConstSplitComplex<float> x, const int xStride) {
//		SMALL_N_STRIDE(2, xStride|yStride, Base::norm2(N, x, xStride, y, yStride))
		float real2 = 0, imag2 = 0;
		vDSP_svesq(x.real, std::abs(xStride), &real2, N);
		vDSP_svesq(x.imag, std::abs(xStride), &imag2, N);
		return real2 + imag2;
	}

	void norm2(const int N, const std::complex<float> *x, const int xStride, float *y, const int yStride) {
//		SMALL_N_STRIDE(2, xStride|yStride, Base::norm2(N, x, xStride, y, yStride))
		if (xStride < 0) x += (1 - N)*xStride;
		if (yStride < 0) y += (1 - N)*yStride;
		
		auto splitX = dspSplit(x);
		vDSP_zvmags(&splitX, 2*xStride, y, yStride, N);
	}

	void norm2(const int N, ConstSplitComplex<float> x, const int xStride, float *y, const int yStride) {
		if (xStride < 0) x += (1 - N)*xStride;
		if (yStride < 0) y += (1 - N)*yStride;

		auto splitX = dspSplit(x);
		vDSP_zvmags(&splitX, xStride, y, yStride, N);
	}

	void reserve(size_t max) {}
private:
	static DSPSplitComplex dspSplit(const ConstSplitComplex<float> &x) {
		DSPSplitComplex dsp;
		dsp.realp = (float *)x.real;
		dsp.imagp = (float *)x.imag;
		return dsp;
	}
	static DSPSplitComplex dspSplit(const SplitComplex<float> &x) {
		DSPSplitComplex dsp;
		dsp.realp = x.real;
		dsp.imagp = x.imag;
		return dsp;
	}
	static DSPSplitComplex dspSplit(const std::complex<float> *x) {
		DSPSplitComplex dsp;
		dsp.realp = (float *)x;
		dsp.imagp = (float *)x + 1;
		return dsp;
	}
	static DSPSplitComplex dspSplit(std::complex<float> *x) {
		DSPSplitComplex dsp;
		dsp.realp = (float *)x;
		dsp.imagp = (float *)x + 1;
		return dsp;
	}
};

template<>
struct Linear<double> : public LinearImplBase<double> {
	using Base = LinearImplBase<double>;

	Linear() : Base(this) {}

	using Base::copy;
#ifdef SIGNALSMITH_USE_CBLAS
	void copy(const int N, const double *x, const int xStride, double *y, const int yStride) {
		SMALL_N_STRIDE(1, xStride|yStride, Base::copy(N, x, xStride, y, yStride))
		cblas_dcopy(CBLAS_INT(N), x, CBLAS_INT(xStride), y, CBLAS_INT(yStride));
	}
	void copy(const int N, const std::complex<double> *x, const int xStride, std::complex<double> *y, const int yStride) {
		SMALL_N_STRIDE(2, xStride|yStride, Base::copy(N, x, xStride, y, yStride))
		cblas_zcopy(CBLAS_INT(N), x, CBLAS_INT(xStride), y, CBLAS_INT(yStride));
	}
#endif

	using Base::norm2;
	double norm2(const int N, const double *x, const int xStride) {
//		SMALL_N_STRIDE(2, xStride, Base::norm2(N, x, xStride))
		double sum2 = 0;
		vDSP_svesqD(x, std::abs(xStride), &sum2, N);
		return sum2;
	}
	double norm2(const int N, const std::complex<double> *x, const int xStride) {
		if (std::abs(xStride) == 1) {
			double sum2 = 0;
			// Aliasing crimes, but this specific one is as kosher as it gets
			vDSP_svesqD((const double *)x, std::abs(xStride), &sum2, N*2);
			return sum2;
		} else {
			// criiiimes
			double real2 = 0, imag2 = 0;
			vDSP_svesqD((const double *)x, std::abs(xStride)*2, &real2, N);
			vDSP_svesqD((const double *)x + 1, std::abs(xStride)*2, &imag2, N);
			return real2 + imag2;
		}
	}
	double norm2(const int N, ConstSplitComplex<double> x, const int xStride) {
		double real2 = 0, imag2 = 0;
		vDSP_svesqD(x.real, std::abs(xStride), &real2, N);
		vDSP_svesqD(x.imag, std::abs(xStride), &imag2, N);
		return real2 + imag2;
	}

	void norm2(const int N, const std::complex<double> *x, const int xStride, double *y, const int yStride) {
		if (xStride < 0) x += (1 - N)*xStride;
		if (yStride < 0) y += (1 - N)*yStride;
		
		auto splitX = dspSplit(x);
		vDSP_zvmagsD(&splitX, 2*xStride, y, yStride, N);
	}

	void norm2(const int N, ConstSplitComplex<double> x, const int xStride, double *y, const int yStride) {
		if (xStride < 0) x += (1 - N)*xStride;
		if (yStride < 0) y += (1 - N)*yStride;

		auto splitX = dspSplit(x);
		vDSP_zvmagsD(&splitX, xStride, y, yStride, N);
	}
	
	void reserve(size_t max) {}
private:
	static DSPDoubleSplitComplex dspSplit(const ConstSplitComplex<double> &x) {
		DSPDoubleSplitComplex dsp;
		dsp.realp = (double *)x.real;
		dsp.imagp = (double *)x.imag;
		return dsp;
	}
	static DSPDoubleSplitComplex dspSplit(const SplitComplex<double> &x) {
		DSPDoubleSplitComplex dsp;
		dsp.realp = x.real;
		dsp.imagp = x.imag;
		return dsp;
	}
	static DSPDoubleSplitComplex dspSplit(const std::complex<double> *x) {
		DSPDoubleSplitComplex dsp;
		dsp.realp = (double *)x;
		dsp.imagp = (double *)x + 1;
		return dsp;
	}
	static DSPDoubleSplitComplex dspSplit(std::complex<double> *x) {
		DSPDoubleSplitComplex dsp;
		dsp.realp = (double *)x;
		dsp.imagp = (double *)x + 1;
		return dsp;
	}
};

}}; // namespace

#undef SMALL_N_STRIDE
#undef CBLAS_INT
#undef BLAS_INDEX
