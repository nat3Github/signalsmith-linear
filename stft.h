#ifndef SIGNALSMITH_AUDIO_LINEAR_STFT_H
#define SIGNALSMITH_AUDIO_LINEAR_STFT_H

#include "./fft.h"

namespace signalsmith { namespace linear {

/// A self-normalising STFT, with variable position/window for output blocks
template<typename Sample, bool splitComputation=true>
struct DynamicSTFT {
	RealFFT<Sample, splitComputation> fft;

	using Complex = std::complex<Sample>;

	void configure(size_t inChannels, size_t outChannels, size_t blockSamples, size_t intervalSamples=0, size_t extraInputHistory=0) {
		_analysisChannels = inChannels;
		_synthesisChannels = outChannels;
		_blockSamples = blockSamples;
		_fftSamples = fft.fastSizeAbove(blockSamples/2)*2;
		fft.resize(_fftSamples);
		_fftBins = _fftSamples/2;
		
		_inputLengthSamples = _blockSamples + extraInputHistory;
		inputBuffer.resize(_inputLengthSamples*_analysisChannels);

		sumBuffer.resize(_blockSamples*_synthesisChannels);
		sumWindowProducts.resize(_blockSamples);
		spectrumBuffer.resize(_fftBins*std::max(_analysisChannels, _synthesisChannels));
		timeBuffer.resize(_fftSamples);

		_analysisWindow.resize(_blockSamples);
		_synthesisWindow.resize(_blockSamples);
		setInterval(intervalSamples ? intervalSamples : blockSamples/4);

		reset();
	}

	void reset() {
		inputPos = _inputLengthSamples;
		outputPos = 0;
		for (auto &v : inputBuffer) v = 0;
		for (auto &v : sumBuffer) v = 0;
		for (auto &v : sumWindowProducts) v = 1;
		for (auto &v : spectrumBuffer) v = 0;
	}

	void writeInput(size_t channel, size_t offset, size_t length, Sample *input) {
		Sample *buffer = inputBuffer.data() + channel*_inputLengthSamples;
		for (size_t i = 0; i < length; ++i) {
			size_t i2 = (offset + inputPos + i)%_inputLengthSamples;
			buffer[i2] = input[i];
		}
	}
	void moveInput(std::ptrdiff_t samples) {
		inputPos = (inputPos + samples)%_blockSamples;
	}

	void readOutput(size_t channel, size_t offset, size_t length, Sample *output) {
		Sample *buffer = sumBuffer.data() + channel*_blockSamples;
		for (size_t i = 0; i < length; ++i) {
			size_t i2 = (offset + outputPos + i)%_blockSamples;
			output[i] = buffer[i2]/sumWindowProducts[i2];
		}
	}
	
	size_t bands() const {
		return _fftBins;
	}
	
	Complex * spectrum(size_t channel) {
		return spectrumBuffer.data() + channel*_fftBins;
	}
	const Complex * spectrum(size_t channel) const {
		return spectrumBuffer.data() + channel*_fftBins;
	}

	Sample * analysisWindow() {
		return _analysisWindow.data();
	}
	const Sample * analysisWindow() const {
		return _analysisWindow.data();
	}
	// Sets the centre index of the window
	void analysisOffset(size_t offset) {
		_analysisOffset = offset;
	}
	size_t analysisOffset() const {
		return _analysisOffset;
	}
	Sample * synthesisWindow() {
		return _synthesisWindow.data();
	}
	const Sample * synthesisWindow() const {
		return _synthesisWindow.data();
	}
	// Sets the centre index of the window
	void synthesisOffset(size_t offset) {
		_synthesisOffset = offset;
	}
	size_t synthesisOffset() const {
		return _synthesisOffset;
	}
	
	void setInterval(size_t defaultInterval, bool updateWindows=true) {
		_defaultInterval = defaultInterval;
		if (!updateWindows) return;
		
		_synthesisOffset = _blockSamples/2;
		_analysisOffset = _blockSamples - _synthesisOffset;
		
		// ACG window with heuristically-good bandwidth (copied from DSP library)
		double gaussianFactor = -0.695*_blockSamples/defaultInterval;
		auto gaussian = [&](double x){
			return std::exp(x*x*gaussianFactor);
		};
		
		double invSize = 1.0/_blockSamples;
		double offsetScale = gaussian(1)/(gaussian(3) + gaussian(1));
		double norm = 1/(gaussian(0) - 2*offsetScale*(gaussian(2)));
		for (auto &v : _synthesisWindow) v = 0;
		for (size_t i = 0; i < _fftSamples; ++i) {
			double r = (2*i + 1)*invSize - 1;
			_analysisWindow[i] = _synthesisWindow[i] = norm*(gaussian(r) - offsetScale*(gaussian(r - 2) + gaussian(r + 2)));
		}
		// Force perfect reconstruction
		for (size_t i = 0; i < defaultInterval; ++i) {
			Sample sum = 0;
			for (size_t i2 = i; i2 < _fftSamples; i2 += defaultInterval) {
				sum += _analysisWindow[i2]*_synthesisWindow[i2];
			}
			Sample factor = 1/std::sqrt(sum);
			for (size_t i2 = i; i2 < _fftSamples; i2 += defaultInterval) {
				_analysisWindow[i2] *= factor;
				_synthesisWindow[i2] *= factor;
			}
		}
	}
	
	void analyse(std::ptrdiff_t samplesInPast=0) {
		for (size_t s = 0; s < analyseSteps(); ++s) {
			analyseStep(s, samplesInPast);
		}
	}
	size_t analyseSteps() const {
		return _analysisChannels*(fft.steps() + 1);
	}
	void analyseStep(size_t step, std::ptrdiff_t samplesInPast=0) {
		size_t fftSteps = fft.steps();
		size_t channel = step/(fftSteps + 1);
		step -= channel*(fftSteps + 1);

		if (step == 0) { // extra step at start of each channel
			Sample *buffer = inputBuffer.data() + channel*_inputLengthSamples;
			for (auto &v : timeBuffer) v = 0;
			for (size_t i = 0; i < _blockSamples; ++i) {
				Sample w = _analysisWindow[i];
				size_t ti = (i + _fftSamples - _synthesisOffset)%_fftSamples;
				size_t bi = (inputPos + i + _inputLengthSamples - _blockSamples - samplesInPast)%_inputLengthSamples;
				timeBuffer[ti] = buffer[bi]*w;
			}
			return;
		}
		--step;
		fft.fft(step, timeBuffer.data(), spectrum(channel));
	}

	void synthesise() {
		synthesise(_defaultInterval);
	}
	void synthesise(size_t samplesSincePrev) {
		for (size_t s = 0; s < synthesiseSteps(); ++s) {
			synthesiseStep(s, samplesSincePrev);
		}
	}
	size_t synthesiseSteps() const {
		return _synthesisChannels*(fft.steps() + 1) + 1;
	}
	void synthesiseStep(size_t step) {
		synthesiseStep(step, _defaultInterval);
	}
	void synthesiseStep(size_t step, size_t samplesSincePrev) {
		if (step == 0) { // Extra first step which adds in the effective gain for a pure analysis-synthesis cycle
			// Before we move the output, clear the old data
			for (size_t i = 0; i < samplesSincePrev; ++i) {
				size_t i2 = (outputPos + i)%_blockSamples;
				for (size_t c = 0; c < _synthesisChannels; ++c) {
					Sample *buffer = sumBuffer.data() + c*_blockSamples;
					buffer[i2] = 0;
				}
				sumWindowProducts[i2] = almostZero;
			}

			outputPos = (outputPos + samplesSincePrev)%_blockSamples;
			
			int windowShift = int(_synthesisOffset) - int(_analysisOffset);
			int wMin = std::max<int>(0, windowShift);
			int wMax = std::min<int>(_blockSamples, int(_blockSamples) + windowShift);

			Sample *windowProduct = sumWindowProducts.data();
			for (int i = wMin; i < wMax; ++i) {
				Sample wa = _analysisWindow[i - windowShift];
				Sample ws = _synthesisWindow[i];
				size_t bi = (outputPos + i)%_blockSamples;
				windowProduct[bi] += wa*ws*fft.size();
			}
			return;
		}
		--step;

		size_t fftSteps = fft.steps();
		size_t channel = step/(fftSteps + 1);
		step -= channel*(fftSteps + 1);

		if (step == fftSteps) { // extra step after each channel's FFT
			Sample *buffer = sumBuffer.data() + channel*_blockSamples;
			for (size_t i = 0; i < _blockSamples; ++i) {
				Sample w = _synthesisWindow[i];
				size_t ti = (i + _fftSamples - _synthesisOffset)%_fftSamples;
				size_t bi = (outputPos + i)%_blockSamples;
				buffer[bi] += timeBuffer[ti]*w;
			}
		} else {
			fft.ifft(step, spectrum(channel), timeBuffer.data());
		}
	}

#define COMPAT_SPELLING(name, alt) \
	template<class ...Args> \
	void alt(Args &&...args) { \
		name(std::forward<Args>(args)...); \
	}
	COMPAT_SPELLING(analyse, analyze);
	COMPAT_SPELLING(analyseStep, analyseStep);
	COMPAT_SPELLING(analyseSteps, analyzeSteps);
	COMPAT_SPELLING(synthesise, synthesize);
	COMPAT_SPELLING(synthesiseStep, synthesizeStep);
	COMPAT_SPELLING(synthesiseSteps, synthesizeSteps);
private:
	static constexpr Sample almostZero = 1e-30;

	size_t _analysisChannels, _synthesisChannels, _inputLengthSamples, _blockSamples, _fftSamples, _fftBins;
	size_t _defaultInterval = 0;

	std::vector<Sample> _analysisWindow, _synthesisWindow;
	size_t _analysisOffset = 0, _synthesisOffset = 0;

	std::vector<Complex> spectrumBuffer;
	std::vector<Sample> timeBuffer;

	size_t inputPos = 0;
	std::vector<Sample> inputBuffer;

	size_t outputPos = 0;
	std::vector<Sample> sumBuffer;
	std::vector<Sample> sumWindowProducts;
};

}} // namespace

#endif // include guard
