#ifndef SIGNALSMITH_AUDIO_LINEAR_FFT_H
#define SIGNALSMITH_AUDIO_LINEAR_FFT_H

#include <complex>
#include <vector>
#include <cmath>

#include "./linear.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace signalsmith { namespace linear {

/// Extremely simple and portable power-of-2 FFT
template<typename Sample>
struct SimpleFFT {
	using Complex = std::complex<Sample>;
	
	SimpleFFT(int maxSize=0) {
		resize(maxSize);
	}
	
	void resize(int maxSize) {
		twiddles.resize(maxSize/2);
		for (int i = 0; i < maxSize/2; ++i) {
			double twiddlePhase = -2*M_PI*i/maxSize;
			twiddles[i] = {
				Sample(std::cos(twiddlePhase)),
				Sample(std::sin(twiddlePhase))
			};
		}
		working.resize(maxSize);
	}
	
	void fft(int size, const Complex *time, Complex *freq) {
		if (size <= 1) {
			*freq = *time;
			return;
		}
		fftPass<false>(size, 1, time, freq, working.data());
	}

	void ifft(int size, const Complex *freq, Complex *time) {
		if (size <= 1) {
			*time = *freq;
			return;
		}
		fftPass<true>(size, 1, freq, time, working.data());
	}
private:
	std::vector<Complex> twiddles;
	std::vector<Complex> working;

	// Calculate a [size]-point FFT, where each element is a block of [stride] values
	template<bool inverse>
	void fftPass(int size, int stride, const Complex *input, Complex *output, Complex *working) const {
		if (size > 2) {
			// Calculate the two half-size FFTs (odd and even) by doubling the stride
			fftPass<inverse>(size/2, stride*2, input, working, output);
			combine2<inverse>(size, stride, working, output);
		} else {
			// The input can already be considered a 1-point FFT
			combine2<inverse>(size, stride, input, output);
		}
	}
	
	// Combine interleaved even/odd results into a single spectrum
	template<bool inverse>
	void combine2(int size, int stride, const Complex *input, Complex *output) const {
		auto twiddleStep = twiddles.size()*2/size;
		for (int i = 0; i < size/2; ++i) {
			Complex twiddle = twiddles[i*twiddleStep];
			
			const Complex *inputA = input + 2*i*stride;
			const Complex *inputB = input + (2*i + 1)*stride;
			Complex *outputA = output + i*stride;
			Complex *outputB = output + (i + size/2)*stride;
			for (int s = 0; s < stride; ++s) {
				Complex a = inputA[s];
				Complex b = inputB[s]*(inverse ? std::conj(twiddle) : twiddle);
				outputA[s] = a + b;
				outputB[s] = a - b;
			}
		}
	}
};

/// A power-of-2 only FFT, specialised with platform-specific fast implementations where available
template<typename Sample>
struct Pow2FFT : private ::signalsmith::linear::Linear {
	using Complex = std::complex<Sample>;
	using Linear = ::signalsmith::linear::Linear;
	
	Pow2FFT(size_t size=0) {
		resize(size);
	}
	
	void resize(size_t size) {
		_size = size;
		simpleFFT.resize(size);
		tmp.resize(size);
	}
	
	void fft(const Complex *time, Complex *freq) {
		simpleFFT.fft(_size, time, freq);
	}

	void fftStrideTime(size_t stride, const Complex *time, Complex *freq) {
		const Complex *input = time;
		if (stride != 1) {
			input = tmp.data();
			for (size_t i = 0; i < _size; ++i) {
				tmp[i] = time[i*stride];
			}
		}
		simpleFFT.fft(_size, input, freq);
	}

	void ifft(const Complex *freq, Complex *time) {
		simpleFFT.ifft(_size, freq, time);
	}

	void ifftStrideFreq(size_t stride, const Complex *freq, Complex *time) {
		const Complex *input = freq;
		if (stride != 1) {
			input = tmp.data();
			for (size_t i = 0; i < _size; ++i) {
				tmp[i] = freq[i*stride];
			}
		}
		simpleFFT.ifft(_size, input, time);
	}

private:
	size_t _size;
	std::vector<Complex> tmp;
	SimpleFFT<Sample> simpleFFT;
};

/// An FFT which can be computed in chunks
template<typename Sample>
struct SplitFFT : public Pow2FFT<Sample> {
	using Complex = std::complex<Sample>;
	using InnerFFT = Pow2FFT<Sample>;
	using Linear = typename InnerFFT::Linear;
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
		InnerFFT::resize(innerSize);
		
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
			fftStep<false>(stepTypes[s], s, time, freq);
		}
	}

	void ifft(const Complex *freq, Complex *time) {
		for (size_t s = 0; s < stepTypes.size(); ++s) {
			fftStep<true>(stepTypes[s], s, freq, time);
		}
	}
private:
	size_t innerSize, outerSize;
	std::vector<Complex> tmpFreq;
	std::vector<Complex> outerTwiddles;
	std::vector<Complex> dftTwists;

	enum class StepType{firstWithFinal, firstWithoutFinal, middleWithFinal, middleWithoutFinal, finalOrder2, finalOrder3, finalOrder4, finalOrder5};
	std::vector<StepType> stepTypes;
	
	template<bool inverse>
	void fftStep(StepType stepType, size_t s, const Complex *time, Complex *freq) {
		switch (stepType) {
			case (StepType::firstWithFinal): {
				if (inverse) {
					InnerFFT::ifftStrideFreq(outerSize, time, freq);
				} else {
					InnerFFT::fftStrideTime(outerSize, time, freq);
				}
				break;
			}
			case (StepType::firstWithoutFinal): {
				if (inverse) {
					InnerFFT::ifftStrideFreq(outerSize, time, freq);
				} else {
					InnerFFT::fftStrideTime(outerSize, time, freq);
				}
				// We're doing the DFT as part of these passes, so duplicate this one
				for (size_t s = 1; s < outerSize; ++s) {
					auto *offsetFreq = freq + s*innerSize;
					for (size_t i = 0; i < innerSize; ++i) {
						offsetFreq[i] = freq[i];
					}
				}
				break;
			}
			case (StepType::middleWithFinal): {
				if (inverse) {
					InnerFFT::ifftStrideFreq(outerSize, time + s, tmpFreq.data());
				} else {
					InnerFFT::fftStrideTime(outerSize, time + s, tmpFreq.data());
				}
				
				auto *twiddles = outerTwiddles.data() + s*innerSize;
				// We'll do the final DFT in-place, as extra passes
				if (inverse) {
					Linear::wrap(freq + s*innerSize, innerSize) = Linear::wrap(tmpFreq)*Linear::wrap(twiddles).conj();
				} else {
					Linear::wrap(freq + s*innerSize, innerSize) = Linear::wrap(tmpFreq)*Linear::wrap(twiddles);
				}
				break;
			}
			case (StepType::middleWithoutFinal): {
				if (inverse) {
					InnerFFT::ifftStrideFreq(outerSize, time + s, tmpFreq.data());
				} else {
					InnerFFT::fftStrideTime(outerSize, time + s, tmpFreq.data());
				}
				auto *twiddles = outerTwiddles.data() + s*innerSize;

				// We have to do the final DFT right here
				for (size_t i = 0; i < innerSize; ++i) {
					Complex v = tmpFreq[i]*(inverse ? std::conj(twiddles[i]) : twiddles[i]);
					tmpFreq[i] = v;
					freq[i] += v;
				}

				for (size_t f = 1; f < outerSize; ++f) {
					Complex dftTwist = dftTwists[(f*s)%outerSize];
					for (size_t i = 0; i < innerSize; ++i) {
						freq[i + f*innerSize] += tmpFreq[i]*(inverse ? std::conj(dftTwist) : dftTwist);
					}
				}
				break;
			}
			case StepType::finalOrder2:
				finalPass2(freq);
				break;
			case StepType::finalOrder3:
				finalPass3<inverse>(freq);
				break;
			case StepType::finalOrder4:
				finalPass4<inverse>(freq);
				break;
			case StepType::finalOrder5:
				finalPass5<inverse>(freq);
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
	template<bool inverse>
	void finalPass3(Complex *f0) {
		auto *f1 = f0 + innerSize;
		auto *f2 = f0 + innerSize*2;
		const Complex tw1{Sample(-0.5), Sample(-std::sqrt(0.75)*(inverse ? -1 : 1))};
		for (size_t i = 0; i < innerSize; ++i) {
			Complex a = f0[i], b = f1[i], c = f2[i];
			f0[i] = a + b + c;
			f1[i] = a + b*tw1 + c*std::conj(tw1);
			f2[i] = a + b*std::conj(tw1) + c*tw1;
		}
	}
	template<bool inverse>
	void finalPass4(Complex *f0) {
		auto *f1 = f0 + innerSize;
		auto *f2 = f0 + innerSize*2;
		auto *f3 = f0 + innerSize*3;
		for (size_t i = 0; i < innerSize; ++i) {
			Complex a = f0[i], b = f1[i], c = f2[i], d = f3[i];
			
			Complex ac0 = a + c, ac1 = a - c;
			Complex bd0 = b + d, bd1 = b - d;
			f0[i] = ac0 + bd0;
			f1[i] = ac1 + bd1*Complex{0, inverse ? 1 : -1};
			f2[i] = ac0 - bd0;
			f3[i] = ac1 - bd1*Complex{0, inverse ? 1 : -1};
		}
	}
	template<bool inverse>
	void finalPass5(Complex *f0) {
		auto *f1 = f0 + innerSize;
		auto *f2 = f0 + innerSize*2;
		auto *f3 = f0 + innerSize*3;
		auto *f4 = f0 + innerSize*4;
		const Complex tw1{0.30901699437494745, -0.9510565162951535*(inverse ? -1 : 1)}, tw2{-0.8090169943749473, -0.5877852522924732*(inverse ? -1 : 1)};
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

}} // namespace

// Platform-specific
#if defined(SIGNALSMITH_USE_ACCELERATE)
#	include "./platform/fft-accelerate.h"
#elif defined(SIGNALSMITH_USE_IPP)
#	include "./platform/fft-ipp.h"
#endif

#endif // include guard
