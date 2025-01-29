#include <Accelerate/Accelerate.h>

namespace signalsmith { namespace linear {

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
