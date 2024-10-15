#ifndef SIGNALSMITH_DSP_FFT2_SPLIT_H
#define SIGNALSMITH_DSP_FFT2_SPLIT_H

#include "./simple-fft.h"

#include <complex>
#include <vector>

namespace signalsmith { namespace fft2 {

/// An FFT which can be computed in chunks
template<typename Sample>
struct SplitFFT {
	using Complex = std::complex<Sample>;
	static constexpr size_t maxSplit = 8;
	static constexpr size_t minInnerSize = 8;
	
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
		tmpTime.resize(innerSize);
		tmpFreq.resize(innerSize);
		innerFFT.resize(innerSize);
		
		outerTwiddles.resize(innerSize*outerSize);
		for (size_t i = 0; i < innerSize; ++i) {
			for (size_t s = 0; s < outerSize; ++s) {
				Sample twiddlePhase = Sample(-2*M_PI*i/innerSize*s/outerSize);
				outerTwiddles[i + s*innerSize] = std::polar(Sample(1), twiddlePhase);
			}
		}
	}
	
	size_t size() const {
		return innerSize*outerSize;
	}
	size_t steps() const {
		return outerSize + 1;
	}
	
	void fft(const Complex *time, Complex *freq) {
		for (size_t s = 0; s < outerSize; ++s) {
			fftStep(s, time, freq);
		}
	}

	void ifft(const Complex *freq, Complex *time) {
		abort(); // not implemented
		innerFFT.ifft(innerSize, time, freq);
	}
private:
	size_t innerSize, outerSize;
	std::vector<Complex> tmpTime, tmpFreq;
	std::vector<Complex> outerTwiddles;
	SimpleFFT<Sample> innerFFT;
	
	void fftStep(size_t s, const Complex *time, Complex *freq) {
		if (s == 0) {
			for (size_t i = 0; i < innerSize; ++i) {
				tmpTime[i] = time[i*outerSize];
			}
			innerFFT.fft(innerSize, tmpTime.data(), freq);
			for (size_t s = 1; s < outerSize; ++s) {
				for (size_t i = 0; i < innerSize; ++i) {
					freq[i + s*innerSize] = freq[i];
				}
			}
		} else {
			for (size_t i = 0; i < innerSize; ++i) {
				tmpTime[i] = time[s + i*outerSize];
			}
			innerFFT.fft(innerSize, tmpTime.data(), tmpFreq.data());

			auto *twiddles = outerTwiddles.data() + s*innerSize;
			for (size_t i = 0; i < innerSize; ++i) {
				Complex v = tmpFreq[i]*twiddles[i];
				tmpFreq[i] = v;
				freq[i] += v;
			}

			for (size_t f = 1; f < outerSize; ++f) {
				Sample dftPhase = Sample(-2*M_PI*f/outerSize*s);
				Complex dftTwist = std::polar(Sample(1), dftPhase);

				for (size_t i = 0; i < innerSize; ++i) {
					freq[i + f*innerSize] += tmpFreq[i]*dftTwist;
				}
			}
		}
	}
};

}} // namespace
#endif // include guard
