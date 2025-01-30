#include <Accelerate/Accelerate.h>

namespace signalsmith { namespace linear {

namespace _impl {
	template<>
	void complexMul<float>(std::complex<float> *a, const std::complex<float> *b, const std::complex<float> *c, size_t size) {
		DSPSplitComplex aSplit = {(float *)a, (float *)a + 1};
		DSPSplitComplex bSplit = {(float *)b, (float *)b + 1};
		DSPSplitComplex cSplit = {(float *)c, (float *)c + 1};
		vDSP_zvmul(&cSplit, 2, &bSplit, 2, &aSplit, 2, size, 1);
	}
	template<>
	void complexMulConj<float>(std::complex<float> *a, const std::complex<float> *b, const std::complex<float> *c, size_t size) {
		DSPSplitComplex aSplit = {(float *)a, (float *)a + 1};
		DSPSplitComplex bSplit = {(float *)b, (float *)b + 1};
		DSPSplitComplex cSplit = {(float *)c, (float *)c + 1};
		vDSP_zvmul(&cSplit, 2, &bSplit, 2, &aSplit, 2, size, -1);
	}
	template<>
	void complexMul<float>(float *ar, float *ai, const float *br, const float *bi, const float *cr, const float *ci, size_t size) {
		DSPSplitComplex aSplit = {ar, ai};
		DSPSplitComplex bSplit = {(float *)br, (float *)bi};
		DSPSplitComplex cSplit = {(float *)cr, (float *)ci};
		vDSP_zvmul(&cSplit, 1, &bSplit, 1, &aSplit, 1, size, 1);
	}
	template<>
	void complexMulConj<float>(float *ar, float *ai, const float *br, const float *bi, const float *cr, const float *ci, size_t size) {
		DSPSplitComplex aSplit = {ar, ai};
		DSPSplitComplex bSplit = {(float *)br, (float *)bi};
		DSPSplitComplex cSplit = {(float *)cr, (float *)ci};
		vDSP_zvmul(&cSplit, 1, &bSplit, 1, &aSplit, 1, size, -1);
	}

	/* not faster, because of the strides
	template<>
	void strideCopy<float>(const std::complex<float> *a, size_t aStride, std::complex<float> *b, size_t size) {
		DSPSplitComplex aSplit = {(float *)a, (float *)a + 1};
		DSPSplitComplex bSplit = {(float *)b, (float *)b + 1};
		vDSP_zvmov(&aSplit, 2*aStride, &bSplit, 2, size);
	}
	*/
	template<>
	void strideCopy<float>(const float *ar, const float *ai, size_t aStride, float *br, float *bi, size_t size) {
		DSPSplitComplex aSplit = {(float *)ar, (float *)ai};
		DSPSplitComplex bSplit = {br, bi};
		vDSP_zvmov(&aSplit, aStride, &bSplit, 1, size);
	}
	
	// doubles
	template<>
	void complexMul<double>(std::complex<double> *a, const std::complex<double> *b, const std::complex<double> *c, size_t size) {
		DSPDoubleSplitComplex aSplit = {(double *)a, (double *)a + 1};
		DSPDoubleSplitComplex bSplit = {(double *)b, (double *)b + 1};
		DSPDoubleSplitComplex cSplit = {(double *)c, (double *)c + 1};
		vDSP_zvmulD(&cSplit, 2, &bSplit, 2, &aSplit, 2, size, 1);
	}
	template<>
	void complexMulConj<double>(std::complex<double> *a, const std::complex<double> *b, const std::complex<double> *c, size_t size) {
		DSPDoubleSplitComplex aSplit = {(double *)a, (double *)a + 1};
		DSPDoubleSplitComplex bSplit = {(double *)b, (double *)b + 1};
		DSPDoubleSplitComplex cSplit = {(double *)c, (double *)c + 1};
		vDSP_zvmulD(&cSplit, 2, &bSplit, 2, &aSplit, 2, size, -1);
	}
	template<>
	void complexMul<double>(double *ar, double *ai, const double *br, const double *bi, const double *cr, const double *ci, size_t size) {
		DSPDoubleSplitComplex aSplit = {ar, ai};
		DSPDoubleSplitComplex bSplit = {(double *)br, (double *)bi};
		DSPDoubleSplitComplex cSplit = {(double *)cr, (double *)ci};
		vDSP_zvmulD(&cSplit, 1, &bSplit, 1, &aSplit, 1, size, 1);
	}
	template<>
	void complexMulConj<double>(double *ar, double *ai, const double *br, const double *bi, const double *cr, const double *ci, size_t size) {
		DSPDoubleSplitComplex aSplit = {ar, ai};
		DSPDoubleSplitComplex bSplit = {(double *)br, (double *)bi};
		DSPDoubleSplitComplex cSplit = {(double *)cr, (double *)ci};
		vDSP_zvmulD(&cSplit, 1, &bSplit, 1, &aSplit, 1, size, -1);
	}

	/* not faster, because of the strides
	template<>
	void strideCopy<double>(const std::complex<double> *a, size_t aStride, std::complex<double> *b, size_t size) {
		DSPDoubleSplitComplex aSplit = {(double *)a, (double *)a + 1};
		DSPDoubleSplitComplex bSplit = {(double *)b, (double *)b + 1};
		vDSP_zvmovD(&aSplit, 2*aStride, &bSplit, 2, size);
	}
	*/
	template<>
	void strideCopy<double>(const double *ar, const double *ai, size_t aStride, double *br, double *bi, size_t size) {
		DSPDoubleSplitComplex aSplit = {(double *)ar, (double *)ai};
		DSPDoubleSplitComplex bSplit = {br, bi};
		vDSP_zvmovD(&aSplit, aStride, &bSplit, 1, size);
	}}

template<>
struct Pow2FFT<float> {
	static constexpr bool prefersSplit = true;

	using Complex = std::complex<float>;

	Pow2FFT(size_t size = 0) {
		if (size > 0) resize(size);
	}
	~Pow2FFT() {
		if (hasSetup) vDSP_destroy_fftsetup(fftSetup);
	}

	void resize(size_t size) {
		_size = size;
		if (hasSetup) vDSP_destroy_fftsetup(fftSetup);
		log2 = std::round(std::log2(size));
		fftSetup = vDSP_create_fftsetup(log2, FFT_RADIX2);
		hasSetup = true;

		splitReal.resize(size);
		splitImag.resize(size);
	}

	void fft(const Complex* input, Complex* output) {
		DSPSplitComplex splitComplex{ splitReal.data(), splitImag.data() };
		vDSP_ctoz((DSPComplex*)input, 2, &splitComplex, 1, _size);
		vDSP_fft_zip(fftSetup, &splitComplex, 1, log2, kFFTDirection_Forward);
		vDSP_ztoc(&splitComplex, 1, (DSPComplex*)output, 2, _size);
	}
	void fft(const float *inR, const float *inI, float *outR, float *outI) {
		DSPSplitComplex inSplit{(float *)inR, (float *)inI};
		DSPSplitComplex outSplit{outR, outI};
		vDSP_fft_zop(fftSetup, &inSplit, 1, &outSplit, 1, log2, kFFTDirection_Forward);
	}

	void ifft(const Complex* input, Complex* output) {
		DSPSplitComplex splitComplex{ splitReal.data(), splitImag.data() };
		vDSP_ctoz((DSPComplex*)input, 2, &splitComplex, 1, _size);
		vDSP_fft_zip(fftSetup, &splitComplex, 1, log2, kFFTDirection_Inverse);
		vDSP_ztoc(&splitComplex, 1, (DSPComplex*)output, 2, _size);
	}
	void ifft(const float *inR, const float *inI, float *outR, float *outI) {
		DSPSplitComplex inSplit{(float *)inR, (float *)inI};
		DSPSplitComplex outSplit{outR, outI};
		vDSP_fft_zop(fftSetup, &inSplit, 1, &outSplit, 1, log2, kFFTDirection_Inverse);
	}

private:
	size_t _size = 0;
	bool hasSetup = false;
	FFTSetup fftSetup;
	int log2 = 0;
	std::vector<float> splitReal, splitImag;
};
template<>
struct Pow2FFT<double> {
	using Complex = std::complex<double>;

	Pow2FFT(size_t size = 0) {
		if (size > 0) resize(size);
	}
	~Pow2FFT() {
		if (hasSetup) vDSP_destroy_fftsetupD(fftSetup);
	}

	void resize(size_t size) {
		_size = size;
		if (hasSetup) vDSP_destroy_fftsetupD(fftSetup);
		log2 = std::round(std::log2(size));
		fftSetup = vDSP_create_fftsetupD(log2, FFT_RADIX2);
		hasSetup = true;

		splitReal.resize(size);
		splitImag.resize(size);
	}

	void fft(const Complex* input, Complex* output) {
		DSPDoubleSplitComplex splitComplex{ splitReal.data(), splitImag.data() };
		vDSP_ctozD((DSPDoubleComplex*)input, 2, &splitComplex, 1, _size);
		vDSP_fft_zipD(fftSetup, &splitComplex, 1, log2, kFFTDirection_Forward);
		vDSP_ztocD(&splitComplex, 1, (DSPDoubleComplex*)output, 2, _size);
	}
	void fft(const double *inR, const double *inI, double *outR, double *outI) {
		DSPDoubleSplitComplex inSplit{(double *)inR, (double *)inI};
		DSPDoubleSplitComplex outSplit{outR, outI};
		vDSP_fft_zopD(fftSetup, &inSplit, 1, &outSplit, 1, log2, kFFTDirection_Forward);
	}

	void ifft(const Complex* input, Complex* output) {
		DSPDoubleSplitComplex splitComplex{ splitReal.data(), splitImag.data() };
		vDSP_ctozD((DSPDoubleComplex*)input, 2, &splitComplex, 1, _size);
		vDSP_fft_zipD(fftSetup, &splitComplex, 1, log2, kFFTDirection_Inverse);
		vDSP_ztocD(&splitComplex, 1, (DSPDoubleComplex*)output, 2, _size);
	}
	void ifft(const double *inR, const double *inI, double *outR, double *outI) {
		DSPDoubleSplitComplex inSplit{(double *)inR, (double *)inI};
		DSPDoubleSplitComplex outSplit{outR, outI};
		vDSP_fft_zopD(fftSetup, &inSplit, 1, &outSplit, 1, log2, kFFTDirection_Inverse);
	}

private:
	size_t _size = 0;
	bool hasSetup = false;
	FFTSetupD fftSetup;
	int log2 = 0;
	std::vector<double> splitReal, splitImag;
};

}} // namespace
