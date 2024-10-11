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
	
	SplitFFT(size_t size=0) {
		resize(size);
	}
	
	size_t resize(size_t size) {
		innerFFT.resize(size);
		_size = size;
		_steps = 1;
		return _steps;
	}
	
	size_t size() const {
		return _size;
	}
	size_t steps() const {
		return _steps;
	}
	
	void fft(const Complex *time, Complex *freq) {
		innerFFT.fft(_size, time, freq);
	}

	void ifft(const Complex *freq, Complex *time) {
		innerFFT.ifft(_size, time, freq);
	}
private:
	size_t _size, _steps;
	SimpleFFT<Sample> innerFFT;
};

}} // namespace
#endif // include guard
