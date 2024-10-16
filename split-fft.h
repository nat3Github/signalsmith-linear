#ifndef SIGNALSMITH_DSP_FFT2_SPLIT_H
#define SIGNALSMITH_DSP_FFT2_SPLIT_H

#include "./simple-fft.h"

#include <complex>
#include <vector>

#ifdef SIGNALSMITH_USE_ACCELERATE
#include <Accelerate/Accelerate.h>
#endif

namespace signalsmith { namespace fft2 {

template<typename Sample>
struct Pow2FFT {
	using Complex = std::complex<Sample>;
	
	Pow2FFT(size_t size=0) {
		resize(size);
	}
	
	void resize(size_t size) {
		_size = size;
		innerFFT.resize(size);
		tmpTime.resize(size);
	}
	
	void fft(const Complex *time, Complex *freq) {
		innerFFT.fft(_size, time, freq);
	}
	
	void fftStrideTime(size_t stride, const Complex *time, Complex *freq) {
		const Complex *input = time;
		if (stride != 1) {
			input = tmpTime.data();
			for (size_t i = 0; i < _size; ++i) {
				tmpTime[i] = time[i*stride];
			}
		}
		innerFFT.fft(_size, input, freq);
	}

private:
	size_t _size;
	std::vector<Complex> tmpTime;
	SimpleFFT<Sample> innerFFT;
};

/// An FFT which can be computed in chunks
template<typename Sample>
struct SplitFFT {
	using Complex = std::complex<Sample>;
	static constexpr size_t maxSplit = 4;
	static constexpr size_t minInnerSize = 32;
	
	static size_t fastSizeAbove(size_t size) {
		size_t pow2 = 1;
		while (pow2 < 16 && pow2 < size) pow2 *= 2;
		while (pow2*8 < size) pow2 *= 2;
		size_t multiple = (size + pow2 - 1)/pow2; // will be 1-8
		if (multiple == 7) ++multiple;
		return multiple*pow2;
	}
	
	SplitFFT(size_t size=0) {
		resize(size);
	}
	
	void resize(size_t size) {
		innerSize = 1;
		outerSize = size;
		// Inner size = largest power of 2 such that either the inner size >= maxSize, or we have the target number of splits
		while (!(outerSize&1) && (outerSize > maxSplit || innerSize < minInnerSize)) {
			innerSize *= 2;
			outerSize /= 2;
		}
		tmpFreq.resize(innerSize);
		innerFFT.resize(innerSize);
		
		outerTwiddles.resize(innerSize*outerSize);
		for (size_t i = 0; i < innerSize; ++i) {
			for (size_t s = 0; s < outerSize; ++s) {
				Sample twiddlePhase = Sample(-2*M_PI*i/innerSize*s/outerSize);
				outerTwiddles[i + s*innerSize] = std::polar(Sample(1), twiddlePhase);
			}
		}

		dftTwists.resize(outerSize);
		for (size_t s = 0; s < outerSize; ++s) {
			Sample dftPhase = Sample(-2*M_PI*s/outerSize);
			dftTwists[s] = std::polar(Sample(1), dftPhase);
		}

		StepType finalStep = (StepType)0; // invalid final step
		if (outerSize == 2) finalStep = StepType::finalOrder2;
		if (outerSize == 3) finalStep = StepType::finalOrder3;
		if (outerSize == 4) finalStep = StepType::finalOrder4;
		if (outerSize == 5) finalStep = StepType::finalOrder5;
		
		if (size <= 1) {
			stepTypes.clear();
			if (size > 0) stepTypes.push_back(StepType::firstWithFinal); // This should just copy, but it's a rare enough case to not need an enum value
		} else if (finalStep == (StepType)0) {
			stepTypes.assign(outerSize, StepType::middleWithoutFinal);
			stepTypes[0] = StepType::firstWithoutFinal;
		} else {
			stepTypes.assign(outerSize, StepType::middleWithFinal);
			stepTypes[0] = StepType::firstWithFinal;
			stepTypes.push_back(finalStep);
		}
	}
	
	size_t size() const {
		return innerSize*outerSize;
	}
	size_t steps() const {
		return stepTypes.size();
	}
	
	void fft(const Complex *time, Complex *freq) {
		for (size_t s = 0; s < stepTypes.size(); ++s) {
			fftStep(stepTypes[s], s, time, freq);
		}
	}

	void ifft(const Complex *freq, Complex *time) {
		abort(); // not implemented
		innerFFT.ifft(innerSize, time, freq);
	}
private:
	size_t innerSize, outerSize;
	std::vector<Complex> tmpFreq;
	std::vector<Complex> outerTwiddles;
	std::vector<Complex> dftTwists;

	Pow2FFT<Sample> innerFFT;
	
	enum class StepType{firstWithFinal, firstWithoutFinal, middleWithFinal, middleWithoutFinal, finalOrder2, finalOrder3, finalOrder4, finalOrder5};
	std::vector<StepType> stepTypes;
	
	void fftStep(StepType stepType, size_t s, const Complex *time, Complex *freq) {
		switch (stepType) {
			case (StepType::firstWithFinal): {
				innerFFT.fftStrideTime(outerSize, time, freq);
				break;
			}
			case (StepType::firstWithoutFinal): {
				innerFFT.fftStrideTime(outerSize, time, freq);
				// We're doing the DFT as part of these passes, so duplicate this one
				for (size_t s = 1; s < outerSize; ++s) {
					for (size_t i = 0; i < innerSize; ++i) {
						freq[i + s*innerSize] = freq[i];
					}
				}
				break;
			}
			case (StepType::middleWithFinal): {
				innerFFT.fftStrideTime(outerSize, time + s, tmpFreq.data());
				
				auto *twiddles = outerTwiddles.data() + s*innerSize;
				// We'll do the final DFT in-place, as extra passes
				for (size_t i = 0; i < innerSize; ++i) {
					freq[i + s*innerSize] = tmpFreq[i]*twiddles[i];
				}
				break;
			}
			case (StepType::middleWithoutFinal): {
				innerFFT.fftStrideTime(outerSize, time + s, tmpFreq.data());
				
				auto *twiddles = outerTwiddles.data() + s*innerSize;
				// We have to do the final DFT right here
				for (size_t i = 0; i < innerSize; ++i) {
					Complex v = tmpFreq[i]*twiddles[i];
					tmpFreq[i] = v;
					freq[i] += v;
				}

				for (size_t f = 1; f < outerSize; ++f) {
					Complex dftTwist = dftTwists[(f*s)%outerSize];
					for (size_t i = 0; i < innerSize; ++i) {
						freq[i + f*innerSize] += tmpFreq[i]*dftTwist;
					}
				}
				break;
			}
			case StepType::finalOrder2:
				finalPass2(freq);
				break;
			case StepType::finalOrder3:
				finalPass3(freq);
				break;
			case StepType::finalOrder4:
				finalPass4(freq);
				break;
			case StepType::finalOrder5:
				finalPass5(freq);
				break;
		}
	}
	
	void finalPass2(Complex *f0) {
		auto *f1 = f0 + innerSize;
		for (size_t i = 0; i < innerSize; ++i) {
			Complex a = f0[i], b = f1[i];
			f0[i] = a + b;
			f1[i] = a - b;
		}
	}
	void finalPass3(Complex *f0) {
		auto *f1 = f0 + innerSize;
		auto *f2 = f0 + innerSize*2;
		const Complex tw1{Sample(-0.5), Sample(-std::sqrt(0.75))};
		for (size_t i = 0; i < innerSize; ++i) {
			Complex a = f0[i], b = f1[i], c = f2[i];
			f0[i] = a + b + c;
			f1[i] = a + b*tw1 + c*std::conj(tw1);
			f2[i] = a + b*std::conj(tw1) + c*tw1;
		}
	}
	void finalPass4(Complex *f0) {
		auto *f1 = f0 + innerSize;
		auto *f2 = f0 + innerSize*2;
		auto *f3 = f0 + innerSize*3;
		for (size_t i = 0; i < innerSize; ++i) {
			Complex a = f0[i], b = f1[i], c = f2[i], d = f3[i];
			
			Complex ac0 = a + c, ac1 = a - c;
			Complex bd0 = b + d, bd1 = b - d;
			f0[i] = ac0 + bd0;
			f1[i] = ac1 + bd1*Complex{0, -1};
			f2[i] = ac0 - bd0;
			f3[i] = ac1 - bd1*Complex{0, -1};
		}
	}
	void finalPass5(Complex *f0) {
		auto *f1 = f0 + innerSize;
		auto *f2 = f0 + innerSize*2;
		auto *f3 = f0 + innerSize*3;
		auto *f4 = f0 + innerSize*4;
		const Complex tw1{0.30901699437494745, -0.9510565162951535}, tw2{-0.8090169943749473, -0.5877852522924732};
		for (size_t i = 0; i < innerSize; ++i) {
			Complex a = f0[i], b = f1[i], c = f2[i], d = f3[i], e = f4[i];
			f0[i] = a + b + c + d + e;
			f1[i] = a + b*tw1 + c*tw2 + d*std::conj(tw2) + e*std::conj(tw1);
			f2[i] = a + b*tw2 + c*std::conj(tw1) + d*tw1 + e*std::conj(tw2);
			f3[i] = a + b*std::conj(tw2) + c*tw1 + d*std::conj(tw1) + e*tw2;
			f4[i] = a + b*std::conj(tw1) + c*std::conj(tw2) + d*tw2 + e*tw1;
		}
	}
};

// Accelerate
#ifdef SIGNALSMITH_USE_ACCELERATE
template<>
struct Pow2FFT<float> {
	using Complex = std::complex<float>;
	
	Pow2FFT(size_t size=0) {
		if (size > 0) resize(size);
	}
	~Pow2FFT() {
		if (hasSetup) vDSP_destroy_fftsetup(fftSetup);
	}
	
	void resize(int size) {
		_size = size;
		if (hasSetup) vDSP_destroy_fftsetup(fftSetup);
		log2 = std::round(std::log2(size));
		fftSetup =  vDSP_create_fftsetup(log2, FFT_RADIX2);
		hasSetup = true;
		
		splitReal.resize(size);
		splitImag.resize(size);
	}
	
	void fft(const Complex *input, Complex *output) {
		DSPSplitComplex splitComplex{splitReal.data(), splitImag.data()};
		vDSP_ctoz((DSPComplex *)input, 2, &splitComplex, 1, _size);
		vDSP_fft_zip(fftSetup, &splitComplex, 1, log2, kFFTDirection_Forward);
		vDSP_ztoc(&splitComplex, 1, (DSPComplex *)output, 2, _size);
	}
	
	void fftStrideTime(size_t stride, const Complex *input, Complex *output) {
		DSPSplitComplex splitComplex{splitReal.data(), splitImag.data()};
		vDSP_ctoz((DSPComplex *)input, 2*stride, &splitComplex, 1, _size);
		vDSP_fft_zip(fftSetup, &splitComplex, 1, log2, kFFTDirection_Forward);
		vDSP_ztoc(&splitComplex, 1, (DSPComplex *)output, 2, _size);
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

	Pow2FFT(size_t size=0) {
		if (size > 0) resize(size);
	}
	~Pow2FFT() {
		if (hasSetup) vDSP_destroy_fftsetupD(fftSetup);
	}
	
	void resize(int size) {
		_size = size;
		if (hasSetup) vDSP_destroy_fftsetupD(fftSetup);
		log2 = std::round(std::log2(size));
		fftSetup =  vDSP_create_fftsetupD(log2, FFT_RADIX2);
		hasSetup = true;
		
		splitReal.resize(size);
		splitImag.resize(size);
	}
	
	void fft(const Complex *input, Complex *output) {
		DSPDoubleSplitComplex splitComplex{splitReal.data(), splitImag.data()};
		vDSP_ctozD((DSPDoubleComplex *)input, 2, &splitComplex, 1, _size);
		vDSP_fft_zipD(fftSetup, &splitComplex, 1, log2, kFFTDirection_Forward);
		vDSP_ztocD(&splitComplex, 1, (DSPDoubleComplex *)output, 2, _size);
	}

	void fftStrideTime(size_t stride, const Complex *input, Complex *output) {
		DSPDoubleSplitComplex splitComplex{splitReal.data(), splitImag.data()};
		vDSP_ctozD((DSPDoubleComplex *)input, 2*stride, &splitComplex, 1, _size);
		vDSP_fft_zipD(fftSetup, &splitComplex, 1, log2, kFFTDirection_Forward);
		vDSP_ztocD(&splitComplex, 1, (DSPDoubleComplex *)output, 2, _size);
	}
private:
	size_t _size = 0;
	bool hasSetup = false;
	FFTSetupD fftSetup;
	int log2 = 0;
	std::vector<double> splitReal, splitImag;
};
#endif

}} // namespace
#endif // include guard
