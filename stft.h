#ifndef SIGNALSMITH_AUDIO_LINEAR_STFT_H
#define SIGNALSMITH_AUDIO_LINEAR_STFT_H

#include "./fft.h"

namespace signalsmith { namespace linear {

/// A self-normalising STFT, with variable position/window for output blocks
template<typename Sample, bool splitComputation=true, bool modified=false>
struct DynamicSTFT {
	RealFFT<Sample, splitComputation, modified> fft;

	using Complex = std::complex<Sample>;

	enum class WindowShape {ignore, acg, kaiser};
	static constexpr WindowShape acg = WindowShape::acg;
	static constexpr WindowShape kaiser = WindowShape::kaiser;

	void configure(size_t inChannels, size_t outChannels, size_t blockSamples, size_t extraInputHistory=0, size_t intervalSamples=0) {
		_analysisChannels = inChannels;
		_synthesisChannels = outChannels;
		_blockSamples = blockSamples;
		_fftSamples = fft.fastSizeAbove((blockSamples + 1)/2)*2;
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
		setInterval(intervalSamples ? intervalSamples : blockSamples/4, kaiser);

		reset();
	}
	
	size_t blockSamples() const {
		return _blockSamples;
	}
	size_t defaultInterval() const {
		return _defaultInterval;
	}
	size_t bands() const {
		return _fftBins;
	}
	size_t latency() const {
		return _synthesisOffset + _blockSamples - _analysisOffset;
	}

	void reset(Sample productWeight=1) {
		inputPos = _blockSamples;
		outputPos = 0;
		for (auto &v : inputBuffer) v = 0;
		for (auto &v : sumBuffer) v = 0;
		for (auto &v : spectrumBuffer) v = 0;
		for (auto &v : sumWindowProducts) v = 0;
		addWindowProduct();
		for (auto &v : sumWindowProducts) v = v*productWeight + almostZero;
	}

	void writeInput(size_t channel, size_t offset, size_t length, Sample *input) {
		Sample *buffer = inputBuffer.data() + channel*_inputLengthSamples;
		for (size_t i = 0; i < length; ++i) {
			size_t i2 = (offset + inputPos + i)%_inputLengthSamples;
			buffer[i2] = input[i];
		}
	}
	void writeInput(size_t channel, size_t length, Sample *input) {
		writeInput(channel, 0, length, input);
	}
	void moveInput(std::ptrdiff_t samples, bool clearInput=false) {
		if (clearInput) {
			for (size_t i = 0; i < samples; ++i) {
				size_t i2 = (inputPos + i)%_inputLengthSamples;
				for (size_t c = 0; c < _analysisChannels; ++c) {
					Sample *buffer = inputBuffer.data() + c*_inputLengthSamples;
					buffer[i2] = 0;
				}
			}
		}

		inputPos = (inputPos + samples)%_inputLengthSamples;
		_samplesSinceAnalysis += samples;
	}
	size_t samplesSinceAnalysis() const {
		return _samplesSinceAnalysis;
	}

	void readOutput(size_t channel, size_t offset, size_t length, Sample *output) {
		Sample *buffer = sumBuffer.data() + channel*_blockSamples;
		for (size_t i = 0; i < length; ++i) {
			size_t i2 = (offset + outputPos + i)%_blockSamples;
			output[i] = buffer[i2]/sumWindowProducts[i2];
		}
	}
	void readOutput(size_t channel, size_t length, Sample *output) {
		return readOutput(channel, 0, length, output);
	}
	void moveOutput(size_t samples) {
		// Zero the output buffer as we cross it
		for (size_t i = 0; i < samples; ++i) {
			size_t i2 = (outputPos + i)%_blockSamples;
			for (size_t c = 0; c < _synthesisChannels; ++c) {
				Sample *buffer = sumBuffer.data() + c*_blockSamples;
				buffer[i2] = 0;
			}
			sumWindowProducts[i2] = almostZero;
		}
		outputPos = (outputPos + samples)%_blockSamples;
		_samplesSinceSynthesis += samples;
	}
	size_t samplesSinceSynthesis() const {
		return _samplesSinceSynthesis;
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
	
	void setInterval(size_t defaultInterval, WindowShape windowShape=WindowShape::ignore) {
		_defaultInterval = defaultInterval;
		if (windowShape == WindowShape::ignore) return;
		
		_analysisOffset = _synthesisOffset = _blockSamples/2;

		if (windowShape == acg) {
			auto window = ApproximateConfinedGaussian::withBandwidth(double(_blockSamples)/defaultInterval);
			window.fill(_analysisWindow, _blockSamples);
		} else if (windowShape == kaiser) {
			auto window = Kaiser::withBandwidth(double(_blockSamples)/defaultInterval);
			window.fill(_analysisWindow,  _blockSamples);
		}

		forcePerfectReconstruction(_analysisWindow, _blockSamples, _defaultInterval);
		for (size_t i = 0; i < _blockSamples; ++i) {
			_synthesisWindow[i] = _analysisWindow[i];
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

		if (step == 0) { // extra step at start of each channel: copy windowed input into buffer
			_samplesSinceAnalysis = samplesInPast;
			Sample *buffer = inputBuffer.data() + channel*_inputLengthSamples;
			for (auto &v : timeBuffer) v = 0;
			for (size_t i = 0; i < _blockSamples; ++i) {
				Sample w = _analysisWindow[i];
				size_t ti = (i + _fftSamples - _analysisOffset)%_fftSamples;
				size_t bi = (inputPos + i - _blockSamples - samplesInPast + _inputLengthSamples)%_inputLengthSamples;
				timeBuffer[ti] = buffer[bi]*w;
			}
			return;
		}
		--step;
		fft.fft(step, timeBuffer.data(), spectrum(channel));
	}

	void synthesise() {
		for (size_t s = 0; s < synthesiseSteps(); ++s) {
			synthesiseStep(s);
		}
	}
	size_t synthesiseSteps() const {
		return _synthesisChannels*(fft.steps() + 1) + 1;
	}
	void synthesiseStep(size_t step) {
		if (step == 0) { // Extra first step which adds in the effective gain for a pure analysis-synthesis cycle
			addWindowProduct();
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
	size_t _samplesSinceSynthesis = 0, _samplesSinceAnalysis = 0;
	
	void addWindowProduct() {
		_samplesSinceSynthesis = 0;

		int windowShift = int(_synthesisOffset) - int(_analysisOffset);
		int wMin = std::max<int>(0, windowShift);
		int wMax = std::min<int>(_blockSamples, int(_blockSamples) + windowShift);

		Sample *windowProduct = sumWindowProducts.data();
		for (int i = wMin; i < wMax; ++i) {
			Sample wa = _analysisWindow[i - windowShift];
			Sample ws = _synthesisWindow[i];
			size_t bi = (outputPos + i)%_blockSamples;
			windowProduct[bi] += wa*ws*_fftSamples;
		}
	}
	
	// Copied from DSP library `windows.h`
	class Kaiser {
		inline static double bessel0(double x) {
			const double significanceLimit = 1e-4;
			double result = 0;
			double term = 1;
			double m = 0;
			while (term > significanceLimit) {
				result += term;
				++m;
				term *= (x*x)/(4*m*m);
			}

			return result;
		}
		double beta;
		double invB0;
		
		static double heuristicBandwidth(double bandwidth) {
			return bandwidth + 8/((bandwidth + 3)*(bandwidth + 3)) + 0.25*std::max(3 - bandwidth, 0.0);
		}
	public:
		Kaiser(double beta) : beta(beta), invB0(1/bessel0(beta)) {}

		static Kaiser withBandwidth(double bandwidth, bool heuristicOptimal=false) {
			return Kaiser(bandwidthToBeta(bandwidth, heuristicOptimal));
		}
		static double bandwidthToBeta(double bandwidth, bool heuristicOptimal=false) {
			if (heuristicOptimal) { // Heuristic based on numerical search
				bandwidth = heuristicBandwidth(bandwidth);
			}
			bandwidth = std::max(bandwidth, 2.0);
			double alpha = std::sqrt(bandwidth*bandwidth*0.25 - 1);
			return alpha*M_PI;
		}
		
		static double betaToBandwidth(double beta) {
			double alpha = beta*(1.0/M_PI);
			return 2*std::sqrt(alpha*alpha + 1);
		}
		static double bandwidthToEnergyDb(double bandwidth, bool heuristicOptimal=false) {
			// Horrible heuristic fits
			if (heuristicOptimal) {
				if (bandwidth < 3) bandwidth += (3 - bandwidth)*0.5;
				return 12.9 + -3/(bandwidth + 0.4) - 13.4*bandwidth + (bandwidth < 3)*-9.6*(bandwidth - 3);
			}
			return 10.5 + 15/(bandwidth + 0.4) - 13.25*bandwidth + (bandwidth < 2)*13*(bandwidth - 2);
		}
		static double energyDbToBandwidth(double energyDb, bool heuristicOptimal=false) {
			double bw = 1;
			while (bw < 20 && bandwidthToEnergyDb(bw, heuristicOptimal) > energyDb) {
				bw *= 2;
			}
			double step = bw/2;
			while (step > 0.0001) {
				if (bandwidthToEnergyDb(bw, heuristicOptimal) > energyDb) {
					bw += step;
				} else {
					bw -= step;
				}
				step *= 0.5;
			}
			return bw;
		}
		static double bandwidthToPeakDb(double bandwidth, bool heuristicOptimal=false) {
			// Horrible heuristic fits
			if (heuristicOptimal) {
				return 14.2 - 20/(bandwidth + 1) - 13*bandwidth + (bandwidth < 3)*-6*(bandwidth - 3) + (bandwidth < 2.25)*5.8*(bandwidth - 2.25);
			}
			return 10 + 8/(bandwidth + 2) - 12.75*bandwidth + (bandwidth < 2)*4*(bandwidth - 2);
		}
		static double peakDbToBandwidth(double peakDb, bool heuristicOptimal=false) {
			double bw = 1;
			while (bw < 20 && bandwidthToPeakDb(bw, heuristicOptimal) > peakDb) {
				bw *= 2;
			}
			double step = bw/2;
			while (step > 0.0001) {
				if (bandwidthToPeakDb(bw, heuristicOptimal) > peakDb) {
					bw += step;
				} else {
					bw -= step;
				}
				step *= 0.5;
			}
			return bw;
		}

		static double bandwidthToEnbw(double bandwidth, bool heuristicOptimal=false) {
			if (heuristicOptimal) bandwidth = heuristicBandwidth(bandwidth);
			double b2 = std::max<double>(bandwidth - 2, 0);
			return 1 + b2*(0.2 + b2*(-0.005 + b2*(-0.000005 + b2*0.0000022)));
		}

		double operator ()(double unit) {
			double r = 2*unit - 1;
			double arg = std::sqrt(1 - r*r);
			return bessel0(beta*arg)*invB0;
		}
	
		template<typename Data>
		void fill(Data &&data, int size) const {
			double invSize = 1.0/size;
			for (int i = 0; i < size; ++i) {
				double r = (2*i + 1)*invSize - 1;
				double arg = std::sqrt(1 - r*r);
				data[i] = bessel0(beta*arg)*invB0;
			}
		}
	};

	class ApproximateConfinedGaussian {
		double gaussianFactor;
		
		double gaussian(double x) const {
			return std::exp(-x*x*gaussianFactor);
		}
	public:
		static double bandwidthToSigma(double bandwidth) {
			return 0.3/std::sqrt(bandwidth);
		}
		static ApproximateConfinedGaussian withBandwidth(double bandwidth) {
			return ApproximateConfinedGaussian(bandwidthToSigma(bandwidth));
		}

		ApproximateConfinedGaussian(double sigma) : gaussianFactor(0.0625/(sigma*sigma)) {}
	
		/// Fills an arbitrary container
		template<typename Data>
		void fill(Data &&data, int size) const {
			double invSize = 1.0/size;
			double offsetScale = gaussian(1)/(gaussian(3) + gaussian(-1));
			double norm = 1/(gaussian(0) - 2*offsetScale*(gaussian(2)));
			for (int i = 0; i < size; ++i) {
				double r = (2*i + 1)*invSize - 1;
				data[i] = norm*(gaussian(r) - offsetScale*(gaussian(r - 2) + gaussian(r + 2)));
			}
		}
	};

	template<typename Data>
	void forcePerfectReconstruction(Data &&data, int windowLength, int interval) {
		for (int i = 0; i < interval; ++i) {
			double sum2 = 0;
			for (int index = i; index < windowLength; index += interval) {
				sum2 += data[index]*data[index];
			}
			double factor = 1/std::sqrt(sum2);
			for (int index = i; index < windowLength; index += interval) {
				data[index] *= factor;
			}
		}
	}
};

}} // namespace

#endif // include guard
