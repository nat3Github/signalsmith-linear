#ifndef SIGNALSMITH_AUDIO_LINEAR_FFT_H
#define SIGNALSMITH_AUDIO_LINEAR_FFT_H

#include <complex>
#include <vector>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace signalsmith { namespace linear {

namespace _impl {
	template<class V>
	void complexMul(std::complex<V> *a, const std::complex<V> *b, const std::complex<V> *c, size_t size) {
		for (size_t i = 0; i < size; ++i) {
			a[i] = b[i]*c[i];
		}
	}
	template<class V>
	void complexMulConj(std::complex<V> *a, const std::complex<V> *b, const std::complex<V> *c, size_t size) {
		for (size_t i = 0; i < size; ++i) {
			a[i] = b[i]*std::conj(c[i]);
		}
	}
	template<class V>
	void complexMul(V *ar, V *ai, const V *br, const V *bi, const V *cr, const V *ci, size_t size) {
		for (size_t i = 0; i < size; ++i) {
			V rr = br[i]*cr[i] - bi[i]*ci[i];
			V ri = br[i]*ci[i] + bi[i]*cr[i];
			ar[i] = rr;
			ai[i] = ri;
		}
	}
	template<class V>
	void complexMulConj(V *ar, V *ai, const V *br, const V *bi, const V *cr, const V *ci, size_t size) {
		for (size_t i = 0; i < size; ++i) {
			V rr = cr[i]*br[i] + ci[i]*bi[i];
			V ri = cr[i]*bi[i] - ci[i]*br[i];
			ar[i] = rr;
			ai[i] = ri;
		}
	}
}

/// Extremely simple and portable power-of-2 FFT
template<typename Sample>
struct SimpleFFT {
	using Complex = std::complex<Sample>;
	
	SimpleFFT(size_t maxSize=0) {
		resize(maxSize);
	}
	
	void resize(size_t maxSize) {
		twiddles.resize(maxSize/2);
		for (size_t i = 0; i < maxSize/2; ++i) {
			double twiddlePhase = -2*M_PI*i/maxSize;
			twiddles[i] = {
				Sample(std::cos(twiddlePhase)),
				Sample(std::sin(twiddlePhase))
			};
		}
		working.resize(maxSize);
	}
	
	void fft(size_t size, const Complex *time, Complex *freq) {
		if (size <= 1) {
			*freq = *time;
			return;
		}
		fftPass<false>(size, 1, time, freq, working.data());
	}

	void ifft(size_t size, const Complex *freq, Complex *time) {
		if (size <= 1) {
			*time = *freq;
			return;
		}
		fftPass<true>(size, 1, freq, time, working.data());
	}

	void fft(size_t size, const Sample *inR, const Sample *inI, Sample *outR, Sample *outI) {
		if (size <= 1) {
			*outR = *inR;
			*outI = *inI;
			return;
		}
		Sample *workingR = (Sample *)working.data(), *workingI = workingR + size;
		fftPass<false>(size, 1, inR, inI, outR, outI, workingR, workingI);
	}
	void ifft(size_t size, const Sample *inR, const Sample *inI, Sample *outR, Sample *outI) {
		if (size <= 1) {
			*outR = *inR;
			*outI = *inI;
			return;
		}
		Sample *workingR = (Sample *)working.data(), *workingI = workingR + size;
		fftPass<true>(size, 1, inR, inI, outR, outI, workingR, workingI);
	}
private:
	std::vector<Complex> twiddles;
	std::vector<Complex> working;

	// Calculate a [size]-point FFT, where each element is a block of [stride] values
	template<bool inverse>
	void fftPass(size_t size, size_t stride, const Complex *input, Complex *output, Complex *working) const {
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
	void combine2(size_t size, size_t stride, const Complex *input, Complex *output) const {
		auto twiddleStep = twiddles.size()*2/size;
		for (size_t i = 0; i < size/2; ++i) {
			Complex twiddle = twiddles[i*twiddleStep];
			
			const Complex *inputA = input + 2*i*stride;
			const Complex *inputB = input + (2*i + 1)*stride;
			Complex *outputA = output + i*stride;
			Complex *outputB = output + (i + size/2)*stride;
			for (size_t s = 0; s < stride; ++s) {
				Complex a = inputA[s];
				Complex b = inputB[s];
				b = inverse ? Complex{b.real()*twiddle.real() + b.imag()*twiddle.imag(), b.imag()*twiddle.real() - b.real()*twiddle.imag()} : Complex{b.real()*twiddle.real() - b.imag()*twiddle.imag(), b.imag()*twiddle.real() + b.real()*twiddle.imag()};
				outputA[s] = a + b;
				outputB[s] = a - b;
			}
		}
	}

	// Calculate a [size]-point FFT, where each element is a block of [stride] values
	template<bool inverse>
	void fftPass(size_t size, size_t stride, const Sample *inputR, const Sample *inputI, Sample *outputR, Sample *outputI, Sample *workingR, Sample *workingI) const {
		if (size > 2) {
			// Calculate the two half-size FFTs (odd and even) by doubling the stride
			fftPass<inverse>(size/2, stride*2, inputR, inputI, workingR, workingI, outputR, outputI);
			combine2<inverse>(size, stride, workingR, workingI, outputR, outputI);
		} else {
			// The input can already be considered a 1-point FFT
			combine2<inverse>(size, stride, inputR, inputI, outputR, outputI);
		}
	}

	// Combine interleaved even/odd results into a single spectrum
	template<bool inverse>
	void combine2(size_t size, size_t stride, const Sample *inputR, const Sample *inputI, Sample *outputR, Sample *outputI) const {
		auto twiddleStep = twiddles.size()*2/size;
		for (size_t i = 0; i < size/2; ++i) {
			Complex twiddle = twiddles[i*twiddleStep];
			
			const Sample *inputAR = inputR + 2*i*stride;
			const Sample *inputAI = inputI + 2*i*stride;
			const Sample *inputBR = inputR + (2*i + 1)*stride;
			const Sample *inputBI = inputI + (2*i + 1)*stride;
			Sample *outputAR = outputR + i*stride;
			Sample *outputAI = outputI + i*stride;
			Sample *outputBR = outputR + (i + size/2)*stride;
			Sample *outputBI = outputI + (i + size/2)*stride;
			for (size_t s = 0; s < stride; ++s) {
				Complex a = {inputAR[s], inputAI[s]};
				Complex b = {inputBR[s], inputBI[s]};
				b = inverse ? Complex{b.real()*twiddle.real() + b.imag()*twiddle.imag(), b.imag()*twiddle.real() - b.real()*twiddle.imag()} : Complex{b.real()*twiddle.real() - b.imag()*twiddle.imag(), b.imag()*twiddle.real() + b.real()*twiddle.imag()};
				Complex sum = a + b, diff = a - b;
				outputAR[s] = sum.real();
				outputAI[s] = sum.imag();
				outputBR[s] = diff.real();
				outputBI[s] = diff.imag();
			}
		}
	}
};

/// A power-of-2 only FFT, specialised with platform-specific fast implementations where available
template<typename Sample>
struct Pow2FFT {
	static constexpr bool prefersSplit = true; // whether this FFT implementation is faster when given split-complex inputs
	using Complex = std::complex<Sample>;
	
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
	void fft(const Sample *inR, const Sample *inI, Sample *outR, Sample *outI) {
		simpleFFT.fft(_size, inR, inI, outR, outI);
	}

	void ifft(const Complex *freq, Complex *time) {
		simpleFFT.ifft(_size, freq, time);
	}
	void ifft(const Sample *inR, const Sample *inI, Sample *outR, Sample *outI) {
		simpleFFT.ifft(_size, inR, inI, outR, outI);
	}

private:
	size_t _size;
	std::vector<Complex> tmp;
	SimpleFFT<Sample> simpleFFT;
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
		// Inner size = largest power of 2 such that either the inner size >= minInnerSize, or we have the target number of splits
		while (!(outerSize&1) && (outerSize > maxSplit || innerSize < minInnerSize)) {
			innerSize *= 2;
			outerSize /= 2;
		}
		tmpStride.resize(0);
		tmpFreq.resize(innerSize);
		innerFFT.resize(innerSize);
		
		outerTwiddles.resize(innerSize*outerSize);
		outerTwiddlesR.resize(innerSize*outerSize);
		outerTwiddlesI.resize(innerSize*outerSize);
		for (size_t i = 0; i < innerSize; ++i) {
			for (size_t s = 0; s < outerSize; ++s) {
				Sample twiddlePhase = Sample(-2*M_PI*i/innerSize*s/outerSize);
				outerTwiddles[i + s*innerSize] = std::polar(Sample(1), twiddlePhase);
			}
		}
		for (size_t i = 0; i < size; ++i) {
			outerTwiddlesR[i] = outerTwiddles[i].real();
			outerTwiddlesI[i] = outerTwiddles[i].imag();
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
			tmpStride.resize(innerSize); // We need more temporary data in this case
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
			fftStep<false>(s, time, freq);
		}
	}
	void fft(size_t step, const Complex *time, Complex *freq) {
		fftStep<false>(step, time, freq);
	}
	void fft(const Sample *inR, const Sample *inI, Sample *outR, Sample *outI) {
		for (size_t s = 0; s < stepTypes.size(); ++s) {
			fftStep<false>(s, inR, inI, outR, outI);
		}
	}
	void fft(size_t step, const Sample *inR, const Sample *inI, Sample *outR, Sample *outI) {
		fftStep<false>(step, inR, inI, outR, outI);
	}
	void ifft(const Complex *freq, Complex *time) {
		for (size_t s = 0; s < stepTypes.size(); ++s) {
			fftStep<true>(s, freq, time);
		}
	}
	void ifft(size_t step, const Complex *freq, Complex *time) {
		fftStep<true>(step, freq, time);
	}
	void ifft(const Sample *inR, const Sample *inI, Sample *outR, Sample *outI) {
		for (size_t s = 0; s < stepTypes.size(); ++s) {
			fftStep<true>(s, inR, inI, outR, outI);
		}
	}
	void ifft(size_t step, const Sample *inR, const Sample *inI, Sample *outR, Sample *outI) {
		fftStep<true>(step, inR, inI, outR, outI);
	}
private:
	using InnerFFT = Pow2FFT<Sample>;
	InnerFFT innerFFT;

	size_t innerSize, outerSize;
	std::vector<Complex> tmpFreq, tmpStride;
	std::vector<Complex> outerTwiddles;
	std::vector<Sample> outerTwiddlesR, outerTwiddlesI;
	std::vector<Complex> dftTwists;

	enum class StepType{firstWithFinal, firstWithoutFinal, middleWithFinal, middleWithoutFinal, finalOrder2, finalOrder3, finalOrder4, finalOrder5};
	std::vector<StepType> stepTypes;
	
	template<bool inverse>
	void fftStep(size_t s, const Complex *time, Complex *freq) {
		switch (stepTypes[s]) {
			case (StepType::firstWithFinal): {
				for (size_t i = 0; i < innerSize; ++i) {
					tmpFreq[i] = time[i*outerSize];
				}
				if (inverse) {
					innerFFT.ifft(tmpFreq.data(), freq);
				} else {
					innerFFT.fft(tmpFreq.data(), freq);
				}
				break;
			}
			case (StepType::firstWithoutFinal): {
				for (size_t i = 0; i < innerSize; ++i) {
					tmpFreq[i] = time[i*outerSize];
				}
				if (inverse) {
					innerFFT.ifft(tmpFreq.data(), freq);
				} else {
					innerFFT.fft(tmpFreq.data(), freq);
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
				const Complex *offsetTime = time + s;
				Complex *offsetOut = freq + s*innerSize;
				for (size_t i = 0; i < innerSize; ++i) {
					offsetOut[i] = offsetTime[i*outerSize];
				}
				if (inverse) {
					innerFFT.ifft(offsetOut, tmpFreq.data());
				} else {
					innerFFT.fft(offsetOut, tmpFreq.data());
				}
				
				auto *twiddles = outerTwiddles.data() + s*innerSize;
				// We'll do the final DFT in-place, as extra passes
				if (inverse) {
					_impl::complexMulConj(offsetOut, tmpFreq.data(), twiddles, innerSize);
				} else {
					_impl::complexMul(offsetOut, tmpFreq.data(), twiddles, innerSize);
				}
				break;
			}
			case (StepType::middleWithoutFinal): {
				const Complex *offsetTime = time + s;
				for (size_t i = 0; i < innerSize; ++i) {
					tmpStride[i] = offsetTime[i*outerSize];
				}
				if (inverse) {
					innerFFT.ifft(tmpStride.data(), tmpFreq.data());
				} else {
					innerFFT.fft(tmpStride.data(), tmpFreq.data());
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
	template<bool inverse>
	void fftStep(size_t s, const Sample *inR, const Sample *inI, Sample *outR, Sample *outI) {
		Sample *tmpR = (Sample *)tmpFreq.data(), *tmpI = tmpR + tmpFreq.size();
		switch (stepTypes[s]) {
			case (StepType::firstWithFinal): {
				for (size_t i = 0; i < innerSize; ++i) {
					tmpR[i] = inR[i*outerSize];
					tmpI[i] = inI[i*outerSize];
				}
				if (inverse) {
					innerFFT.ifft(tmpR, tmpI, outR, outI);
				} else {
					innerFFT.fft(tmpR, tmpI, outR, outI);
				}
				break;
			}
			case (StepType::firstWithoutFinal): {
				for (size_t i = 0; i < innerSize; ++i) {
					tmpR[i] = inR[i*outerSize];
					tmpI[i] = inI[i*outerSize];
				}
				if (inverse) {
					innerFFT.ifft(tmpR, tmpI, outR, outI);
				} else {
					innerFFT.fft(tmpR, tmpI, outR, outI);
				}
				// We're doing the DFT as part of these passes, so duplicate this one
				for (size_t s = 1; s < outerSize; ++s) {
					auto *offsetR = outR + s*innerSize;
					auto *offsetI = outI + s*innerSize;
					for (size_t i = 0; i < innerSize; ++i) {
						offsetR[i] = outR[i];
						offsetI[i] = outI[i];
					}
				}
				break;
			}
			case (StepType::middleWithFinal): {
				const Sample *offsetR = inR + s;
				const Sample *offsetI = inI + s;
				Sample *offsetOutR = outR + s*innerSize;
				Sample *offsetOutI = outI + s*innerSize;
				for (size_t i = 0; i < innerSize; ++i) {
					offsetOutR[i] = offsetR[i*outerSize];
					offsetOutI[i] = offsetI[i*outerSize];
				}
				if (inverse) {
					innerFFT.ifft(offsetOutR, offsetOutI, tmpR, tmpI);
				} else {
					innerFFT.fft(offsetOutR, offsetOutI, tmpR, tmpI);
				}
				
				auto *twiddlesR = outerTwiddlesR.data() + s*innerSize;
				auto *twiddlesI = outerTwiddlesI.data() + s*innerSize;
				// We'll do the final DFT in-place, as extra passes
				if (inverse) {
					_impl::complexMulConj(offsetOutR, offsetOutI, tmpR, tmpI, twiddlesR, twiddlesI, innerSize);
				} else {
					_impl::complexMul(offsetOutR, offsetOutI, tmpR, tmpI, twiddlesR, twiddlesI, innerSize);
				}
				break;
			}
			case (StepType::middleWithoutFinal): {
				const Sample *offsetR = inR + s;
				const Sample *offsetI = inI + s;
				Sample *tmpStrideR = (Sample *)tmpStride.data();
				Sample *tmpStrideI = tmpStrideR + innerSize;
				for (size_t i = 0; i < innerSize; ++i) {
					tmpStrideR[i] = offsetR[i*outerSize];
					tmpStrideI[i] = offsetI[i*outerSize];
				}
				if (inverse) {
					innerFFT.ifft(tmpStrideR, tmpStrideI, tmpR, tmpI);
				} else {
					innerFFT.fft(tmpStrideR, tmpStrideI, tmpR, tmpI);
				}

				auto *twiddlesR = outerTwiddlesR.data() + s*innerSize;
				auto *twiddlesI = outerTwiddlesI.data() + s*innerSize;
				// We have to do the final DFT right here
				for (size_t i = 0; i < innerSize; ++i) {
					Complex v = Complex{tmpR[i], tmpI[i]}*Complex{twiddlesR[i], inverse ? -twiddlesI[i] : twiddlesI[i]};
					tmpR[i] = v.real();
					tmpI[i] = v.imag();
					outR[i] += v.real();
					outI[i] += v.imag();
				}
				for (size_t f = 1; f < outerSize; ++f) {
					Complex dftTwist = dftTwists[(f*s)%outerSize];
					for (size_t i = 0; i < innerSize; ++i) {
						Complex v = Complex{tmpR[i], tmpI[i]}*(inverse ? std::conj(dftTwist) : dftTwist);
						outR[i + f*innerSize] += v.real();
						outI[i + f*innerSize] += v.imag();
					}
				}
				break;
			}
			case StepType::finalOrder2:
				finalPass2(outR, outI);
				break;
			case StepType::finalOrder3:
				finalPass3<inverse>(outR, outI);
				break;
			case StepType::finalOrder4:
				finalPass4<inverse>(outR, outI);
				break;
			case StepType::finalOrder5:
				finalPass5<inverse>(outR, outI);
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
	void finalPass2(Sample *f0r, Sample *f0i) {
		auto *f1r = f0r + innerSize;
		auto *f1i = f0i + innerSize;
		for (size_t i = 0; i < innerSize; ++i) {
			Sample ar = f0r[i], ai = f0i[i];
			Sample br = f1r[i], bi = f1i[i];
			f0r[i] = ar + br;
			f0i[i] = ai + bi;
			f1r[i] = ar - br;
			f1i[i] = ai - bi;
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
	void finalPass3(Sample *f0r, Sample *f0i) {
		auto *f1r = f0r + innerSize;
		auto *f1i = f0i + innerSize;
		auto *f2r = f0r + innerSize*2;
		auto *f2i = f0i + innerSize*2;
		const Sample tw1r = -0.5, tw1i = -std::sqrt(0.75)*(inverse ? -1 : 1);
		
		for (size_t i = 0; i < innerSize; ++i) {
			Sample ar = f0r[i], ai = f0i[i], br = f1r[i], bi = f1i[i], cr = f2r[i], ci = f2i[i];

			f0r[i] = ar + br + cr;
			f0i[i] = ai + bi + ci;
			f1r[i] = ar + br*tw1r - bi*tw1i + cr*tw1r + ci*tw1i;
			f1i[i] = ai + bi*tw1r + br*tw1i - cr*tw1i + ci*tw1r;
			f2r[i] = ar + br*tw1r + bi*tw1i + cr*tw1r - ci*tw1i;
			f2i[i] = ai + bi*tw1r - br*tw1i + cr*tw1i + ci*tw1r;
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
			Complex bd0 = b + d, bd1 = inverse ? (b - d) : (d - b);
			Complex bd1i = {-bd1.imag(), bd1.real()};
			f0[i] = ac0 + bd0;
			f1[i] = ac1 + bd1i;
			f2[i] = ac0 - bd0;
			f3[i] = ac1 - bd1i;
		}
	}
	template<bool inverse>
	void finalPass4(Sample *f0r, Sample *f0i) {
		auto *f1r = f0r + innerSize;
		auto *f1i = f0i + innerSize;
		auto *f2r = f0r + innerSize*2;
		auto *f2i = f0i + innerSize*2;
		auto *f3r = f0r + innerSize*3;
		auto *f3i = f0i + innerSize*3;
		for (size_t i = 0; i < innerSize; ++i) {
			Sample ar = f0r[i], ai = f0i[i], br = f1r[i], bi = f1i[i], cr = f2r[i], ci = f2i[i], dr = f3r[i], di = f3i[i];
			
			Sample ac0r = ar + cr, ac0i = ai + ci;
			Sample ac1r = ar - cr, ac1i = ai - ci;
			Sample bd0r = br + dr, bd0i = bi + di;
			Sample bd1r = br - dr, bd1i = bi - di;
			
			f0r[i] = ac0r + bd0r;
			f0i[i] = ac0i + bd0i;
			f1r[i] = inverse ? (ac1r - bd1i) : (ac1r + bd1i);
			f1i[i] = inverse ? (ac1i + bd1r) : (ac1i - bd1r);
			f2r[i] = ac0r - bd0r;
			f2i[i] = ac0i - bd0i;
			f3r[i] = inverse ? (ac1r + bd1i) : (ac1r - bd1i);
			f3i[i] = inverse ? (ac1i - bd1r) : (ac1i + bd1r);
		}
	}
	template<bool inverse>
	void finalPass5(Complex *f0) {
		auto *f1 = f0 + innerSize;
		auto *f2 = f0 + innerSize*2;
		auto *f3 = f0 + innerSize*3;
		auto *f4 = f0 + innerSize*4;
		const Sample tw1r = 0.30901699437494745;
		const Sample tw1i = -0.9510565162951535*(inverse ? -1 : 1);
		const Sample tw2r = -0.8090169943749473;
		const Sample tw2i = -0.5877852522924732*(inverse ? -1 : 1);
		for (size_t i = 0; i < innerSize; ++i) {
			Complex a = f0[i], b = f1[i], c = f2[i], d = f3[i], e = f4[i];

			Complex be0 = b + e, be1 = {e.imag() - b.imag(), b.real() - e.real()}; // (b - e)*i
			Complex cd0 = c + d, cd1 = {d.imag() - c.imag(), c.real() - d.real()};
			
			Complex bcde01 = be0*tw1r + cd0*tw2r;
			Complex bcde02 = be0*tw2r + cd0*tw1r;
			Complex bcde11 = be1*tw1i + cd1*tw2i;
			Complex bcde12 = be1*tw2i - cd1*tw1i;

			f0[i] = a + be0 + cd0;
			f1[i] = a + bcde01 + bcde11;
			f2[i] = a + bcde02 + bcde12;
			f3[i] = a + bcde02 - bcde12;
			f4[i] = a + bcde01 - bcde11;
		}
	}
	template<bool inverse>
	void finalPass5(Sample *f0r, Sample *f0i) {
		auto *f1r = f0r + innerSize;
		auto *f1i = f0i + innerSize;
		auto *f2r = f0r + innerSize*2;
		auto *f2i = f0i + innerSize*2;
		auto *f3r = f0r + innerSize*3;
		auto *f3i = f0i + innerSize*3;
		auto *f4r = f0r + innerSize*4;
		auto *f4i = f0i + innerSize*4;
		
		const Sample tw1r = 0.30901699437494745;
		const Sample tw1i = -0.9510565162951535*(inverse ? -1 : 1);
		const Sample tw2r = -0.8090169943749473;
		const Sample tw2i = -0.5877852522924732*(inverse ? -1 : 1);
		for (size_t i = 0; i < innerSize; ++i) {
			Sample ar = f0r[i], ai = f0i[i], br = f1r[i], bi = f1i[i], cr = f2r[i], ci = f2i[i], dr = f3r[i], di = f3i[i], er = f4r[i], ei = f4i[i];

			Sample be0r = br + er, be0i = bi + ei;
			Sample be1r = br - er, be1i = ei - bi;
			Sample cd0r = cr + dr, cd0i = ci + di;
			Sample cd1r = cr - dr, cd1i = di - ci;

			f0r[i] = ar + be0r + cd0r;
			f0i[i] = ai + be0i + cd0i;
			f1r[i] = ar + be0r*tw1r + be1i*tw1i + cd0r*tw2r + cd1i*tw2i;
			f1i[i] = ai + be0i*tw1r + be1r*tw1i + cd0i*tw2r + cd1r*tw2i;
			f2r[i] = ar + be0r*tw2r + be1i*tw2i + cd0r*tw1r - cd1i*tw1i;
			f2i[i] = ai + be0i*tw2r + be1r*tw2i + cd0i*tw1r - cd1r*tw1i;
			f3r[i] = ar + be0r*tw2r - be1i*tw2i + cd0r*tw1r + cd1i*tw1i;
			f3i[i] = ai + be0i*tw2r - be1r*tw2i + cd0i*tw1r + cd1r*tw1i;
			f4r[i] = ar + be0r*tw1r - be1i*tw1i + cd0r*tw2r - cd1i*tw2i;
			f4i[i] = ai + be0i*tw1r - be1r*tw1i + cd0i*tw2r - cd1r*tw2i;
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
