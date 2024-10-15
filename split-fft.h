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
		tmpTime.resize(std::max(innerSize, outerSize));
		innerFFT.resize(innerSize);
	}
	
	size_t size() const {
		return innerSize*outerSize;
	}
	size_t steps() const {
		return outerSize + 1;
	}
	
	void fft(const Complex *time, Complex *freq) {
		for (size_t s = 0; s < outerSize; ++s) {
			for (size_t i = 0; i < innerSize; ++i) {
				tmpTime[i] = time[s + i*outerSize];
			}
			innerFFT.fft(innerSize, tmpTime.data(), freq + s*innerSize);
		}
		if (outerSize > 1) {
			for (size_t i = 0; i < innerSize; ++i) {
				Sample twiddlePhase = Sample(-2*M_PI*i/innerSize/outerSize);

				// Copy and apply twiddles
				for (size_t s = 0; s < outerSize; ++s) {
					Complex twiddle = std::polar(Sample(1), twiddlePhase*s);
					tmpTime[s] = freq[i + s*innerSize]*twiddle;
				}

				// Naive DFT, write results in-place
				for (size_t f = 0; f < outerSize; ++f) {
					Complex sum = 0;
					for (size_t s = 0; s < outerSize; ++s) {
						Sample dftPhase = Sample(-2*M_PI*f/outerSize*s);
						Complex dftTwist = std::polar(Sample(1), dftPhase);
						sum += tmpTime[s]*dftTwist;
					}
					freq[i + f*innerSize] = sum;
				}
			}
		}
	}

	void ifft(const Complex *freq, Complex *time) {
LOG_EXPR("ifft");
		abort(); // not implemented
		innerFFT.ifft(innerSize, time, freq);
	}
private:
	size_t innerSize, outerSize;
	std::vector<Complex> tmpTime;
	SimpleFFT<Sample> innerFFT;
};

}} // namespace
#endif // include guard
