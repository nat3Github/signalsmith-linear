#include "./stopwatch.h"
#if defined(__has_include) && __has_include("plot/signalsmith.h")
#	include "plot/signalsmith.h"
#else
#	include "plot/plot.h"
#endif

#include <iostream>
#include <complex>
#include <vector>
#include <random>

#define LOG_EXPR(expr) std::cout << #expr << " = " << (expr) << std::endl;

template<typename Sample>
struct RunData {
	using Complex = std::complex<Sample>;

	const int size;
	const int maxSize;
	std::vector<Complex> input, output;
	
	RunData(int size, int maxSize, int seed=0) : size(size), maxSize(maxSize), input(size), output(size), randomEngine(seed) {
		randomise();
	}
	
	void randomise() {
		std::uniform_real_distribution<Sample> dist{-1, 1};
		for (auto &v : input) {
			v = {dist(randomEngine), dist(randomEngine)};
		}
	}

private:
	std::default_random_engine randomEngine;
};

template<class FftWrapper>
struct Runner {
	static constexpr double testSeconds = 0.5;
	static constexpr double testChunk = 0.1;

	const char *name;
	signalsmith::plot::Line2D &line;
	Stopwatch stopwatch{false};
	FftWrapper fft;
	
	Runner(const char *name, signalsmith::plot::Line2D &line, signalsmith::plot::Legend &legend) : name(name), line(line) {
		legend.add(line, name);
	}
	
	template<class Data>
	void run(double x, Data &data) {
		fft.prepare(data.size, data.maxSize);
		size_t rounds = 0, roundStep = 1;
		
		double dummySum = 1;
		double seconds = 0;
		while (seconds < testSeconds) {
			stopwatch.start();
			
			for (size_t r = 0; r < roundStep; ++r) {
				dummySum += fft.run(data);
			}
			
			double lap = stopwatch.seconds(stopwatch.lap());
			if (lap < testChunk) {
				roundStep *= 2;
			} else {
				seconds += lap;
				rounds += roundStep;
			}
		}
		double rps = rounds/seconds;
		double ref = data.size*(std::log2(data.size) + 1);
		double scaledRps = rps*ref;
		line.add(x, scaledRps);
		
		std::cout << name << "\t@ " << data.size << ": " << scaledRps << " (" << dummySum << ")\n";
	}
};

// ---------- wrappers

#include "../simple-fft.h"
template<class Sample>
struct SimpleWrapper {
	signalsmith::fft2::SimpleFFT<Sample> fft;

	void prepare(int size, int) {
		fft.resize(size);
	}
	
	template<class Data>
	double run(Data &data) {
		fft.fft(data.size, data.input.data(), data.output.data());
		return data.output[0].real();
	}
};

#include "./others/signalsmith-fft.h"
template<class Sample>
struct SignalsmithFFTWrapper {
	signalsmith::FFT<Sample> fft{1};

	void prepare(int size, int) {
		fft.setSize(size);
	}
	
	template<class Data>
	double run(Data &data) {
		fft.fft(data.input.data(), data.output.data());
		return data.output[0].real();
	}
};

#include "dsp/fft.h"
template<class Sample>
struct SignalsmithDSPWrapper {
	signalsmith::fft::FFT<Sample> fft{1};

	void prepare(int size, int) {
		fft.setSize(size);
	}
	
	template<class Data>
	double run(Data &data) {
		fft.fft(data.input.data(), data.output.data());
		return data.output[0].real();
	}
};

#include <Accelerate/Accelerate.h>
struct AccelerateFloatWrapper {
	bool hasSetup = false;
	FFTSetup fftSetup;
	int log2 = 0;
	
	std::vector<float> splitReal, splitImag;

	AccelerateFloatWrapper() {}
	~AccelerateFloatWrapper() {
		if (hasSetup) vDSP_destroy_fftsetup(fftSetup);
	}
	
	void prepare(int size, int) {
		if (hasSetup) vDSP_destroy_fftsetup(fftSetup);
		log2 = std::round(std::log2(size));
		fftSetup =  vDSP_create_fftsetup(log2, FFT_RADIX2);
		hasSetup = true;
		
		splitReal.resize(size);
		splitImag.resize(size);
	}
	
	template<class Data>
	double run(Data &data) {
		DSPSplitComplex splitComplex{splitReal.data(), splitImag.data()};
		vDSP_ctoz((DSPComplex *)data.input.data(), 2, &splitComplex, 1, data.size);
		
		vDSP_fft_zip(fftSetup, &splitComplex, 1, log2, kFFTDirection_Forward);

		vDSP_ztoc(&splitComplex, 1, (DSPComplex *)data.output.data(), 2, data.size);
		return data.output[0].real();
	}
};
struct AccelerateDoubleWrapper {
	bool hasSetup = false;
	FFTSetupD fftSetup;
	int log2 = 0;
	
	std::vector<double> splitReal, splitImag;

	AccelerateDoubleWrapper() {}
	~AccelerateDoubleWrapper() {
		if (hasSetup) vDSP_destroy_fftsetupD(fftSetup);
	}
	
	void prepare(int size, int) {
		if (hasSetup) vDSP_destroy_fftsetupD(fftSetup);
		log2 = std::round(std::log2(size));
		fftSetup =  vDSP_create_fftsetupD(log2, FFT_RADIX2);
		hasSetup = true;
		
		splitReal.resize(size);
		splitImag.resize(size);
	}
	
	template<class Data>
	double run(Data &data) {
		DSPDoubleSplitComplex splitComplex{splitReal.data(), splitImag.data()};
		vDSP_ctozD((DSPDoubleComplex *)data.input.data(), 2, &splitComplex, 1, data.size);
		
		vDSP_fft_zipD(fftSetup, &splitComplex, 1, log2, kFFTDirection_Forward);

		vDSP_ztocD(&splitComplex, 1, (DSPDoubleComplex *)data.output.data(), 2, data.size);
		return data.output[0].real();
	}
};

// ---------- main code

int main() {
	signalsmith::plot::Figure figure;
	auto &plot = figure.plot(800, 250);
	plot.x.label("FFT size");
	
	auto &legend = plot.legend(0, 1);
	Runner<SimpleWrapper<double>> simpleDouble("simple (double)", plot.line(), legend);
	Runner<SimpleWrapper<float>> simpleFloat("simple (float)", plot.line(), legend);
//	Runner<SignalsmithFFTWrapper<double>> signalsmithDouble("Signalsmith (double)", plot.line(), legend);
//	Runner<SignalsmithFFTWrapper<float>> signalsmithFloat("Signalsmith (float)", plot.line(), legend);
	Runner<SignalsmithDSPWrapper<double>> dspDouble("DSP library (double)", plot.line(), legend);
	Runner<SignalsmithDSPWrapper<float>> dspFloat("DSP library (float)", plot.line(), legend);
	Runner<AccelerateDoubleWrapper> accelerateDouble("Accelerate (double)", plot.line(), legend);
	Runner<AccelerateFloatWrapper> accelerateFloat("Accelerate (float)", plot.line(), legend);

	int maxSize = 65536*8;
	bool first = true;
	for (int n = 1; n <= maxSize; n *= 2) {
		double x = std::log2(n);
		RunData<double> dataDouble(n, maxSize, n);
		RunData<float> dataFloat(n, maxSize, n);
		
		simpleDouble.run(x, dataDouble);
		simpleFloat.run(x, dataFloat);
//		signalsmithDouble.run(x, dataDouble);
//		signalsmithFloat.run(x, dataFloat);
		dspDouble.run(x, dataDouble);
		dspFloat.run(x, dataFloat);
		accelerateDouble.run(x, dataDouble);
		accelerateFloat.run(x, dataFloat);

		if (first) {
			first = false;
			plot.x.major(x, std::to_string(n));
		} else {
			plot.x.tick(x, std::to_string(n));
		}
	}
	plot.y.major(0); // auto-scaled range includes 0
	plot.y.blankLabels().label("speed"); // values don't matter, only the comparison
	figure.write("comparison.svg");
}
