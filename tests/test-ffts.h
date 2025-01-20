#include "../fft.h"

#include "./test-runner.h"

// ---------- wrappers

template<class Sample>
struct SimpleWrapper {
	signalsmith::linear::SimpleFFT<Sample> fft;

	void prepare(int size, int) {
		fft.resize(size);
	}
	
	template<class Data>
	void run(Data &data) {
		auto *time = data.complex(0), *freq = data.complex(1);
		fft.fft(data.size, time, freq);
	}
};

template<class Sample>
struct Pow2Wrapper {
	signalsmith::linear::Pow2FFT<Sample> fft;

	void prepare(int size, int) {
		fft.resize(size);
	}

	template<class Data>
	void run(Data& data) {
		auto *time = data.complex(0), *freq = data.complex(1);
		fft.fft(time, freq);
	}
};

template<class Sample>
struct SplitWrapper {
	signalsmith::linear::SplitFFT<Sample> fft;

	void prepare(int size, int) {
		fft.resize(size);
	}

	template<class Data>
	void run(Data& data) {
		auto *time = data.complex(0), *freq = data.complex(1);
		fft.fft(time, freq);
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
	void run(Data& data) {
		auto *time = data.complex(0), *freq = data.complex(1);
		fft.fft(time, freq);
	}
};

#include "./others/dsp/fft.h"
template<class Sample>
struct SignalsmithDSPWrapper {
	signalsmith::fft::FFT<Sample> fft{1};

	void prepare(int size, int) {
		fft.setSize(size);
	}
	
	template<class Data>
	void run(Data& data) {
		auto *time = data.complex(0), *freq = data.complex(1);
		fft.fft(time, freq);
	}
};

// ---------- main code

void testFfts(int maxSize, double benchmarkSeconds) {
	signalsmith::plot::Plot2D fastSizePlot(200, 200);
	auto &fastSizeLine = fastSizePlot.line();
	for (int n = 1; n < 65536; ++n) {
		size_t fastN = signalsmith::linear::SplitFFT<double>::fastSizeAbove(n);
		fastSizeLine.add(std::log2(n), std::log2(fastN));
	}
	fastSizePlot.line().add(0, 0).add(16, 16);
	fastSizePlot.line().add(0, std::log2(1.25)).add(16 - std::log2(1.25), 16);
	fastSizePlot.write("fast-sizes.svg");
	
	RunPlot plot("ffts", benchmarkSeconds);

	auto simpleDouble = plot.runner<SimpleWrapper<double>>("Simple (double)");
	auto simpleFloat = plot.runner<SimpleWrapper<float>>("Simple (float)");
	auto pow2Double = plot.runner<Pow2Wrapper<double>>("Pow2 (double)");
	auto pow2Float = plot.runner<Pow2Wrapper<float>>("Pow2 (float)");
	auto splitDouble = plot.runner<SplitWrapper<double>>("Split (double)");
	auto splitFloat = plot.runner<SplitWrapper<float>>("Split (float)");
	auto dspDouble = plot.runner<SignalsmithDSPWrapper<double>>("DSP library (double)");
	auto dspFloat = plot.runner<SignalsmithDSPWrapper<float>>("DSP library (float)");

	auto runSize = [&](int n, int pow3=0, int pow5=0, int pow7=0){
		std::cout << "\tn = " << n << "\r" << std::flush;
		double refTime = 1e-8*n*(std::log2(n) + 1); // heuristic for expected computation time, just to compare different sizes

		RunData<double> dataDouble(n);
		RunData<float> dataFloat(n);

		RunData<double> *refPtrDouble = nullptr;
		RunData<float> *refPtrFloat = nullptr;

		auto refDataDouble = dataDouble;
		auto refDataFloat = dataFloat;
		if (n < 256) {
			refPtrDouble = &refDataDouble;
			refPtrFloat = &refDataFloat;
			auto *inputDouble = refDataDouble.complex(0);
			auto *outputDouble = refDataDouble.complex(1);
			auto *inputFloat = refDataFloat.complex(0);
			auto *outputFloat = refDataFloat.complex(1);
			for (int f = 0; f < n; ++f) {
				std::complex<double> sumD = 0;
				std::complex<double> sumF = 0;
				for (int i = 0; i < n; ++i) {
					auto coeff = std::polar(1.0, -2*M_PI*i*f/n);
					sumD += coeff*inputDouble[i];
					std::complex<float> vf = inputFloat[i];
					sumF += coeff*std::complex<double>{vf.real(), vf.imag()};
				}
				outputDouble[f] = sumD;
				outputFloat[f] = {float(sumF.real()), float(sumF.imag())};
			}
		}

		if (pow3 + pow5 + pow7 == 0) {
			simpleDouble.run(dataDouble, refTime, refPtrDouble);
			simpleFloat.run(dataFloat, refTime, refPtrFloat);
			pow2Double.run(dataDouble, refTime, refPtrDouble);
			pow2Float.run(dataFloat, refTime, refPtrFloat);
		}
		if (pow5 + pow7 == 0) {
			dspDouble.run(dataDouble, refTime, refPtrDouble);
			dspFloat.run(dataFloat, refTime, refPtrFloat);
		}
		if (signalsmith::linear::SplitFFT<double>::fastSizeAbove(n) == size_t(n)) {
			splitDouble.run(dataDouble, refTime, refPtrDouble);
			splitFloat.run(dataFloat, refTime, refPtrFloat);
		}
	};
	for (int n = 1; n <= maxSize; n *= 2) {
		if (n/16) runSize(n*9/16, 2);
		if (n/8) runSize(n*5/8, 0, 1);
		if (n/4) runSize(n*3/4, 1);
		if (n/8) runSize(n*7/8, 1, 1);
		if (n/16) runSize(n*15/16, 1, 1);
		runSize(n);

		if (int(std::round(std::log2(n)))%2 == 0) {
			plot.tick(n, 2);
		} else {
			plot.tick(n, "");
		}
	}
	std::cout << "\t                                \r" << std::flush;

	plot.plot.x.range(std::log, 1, maxSize).label("FFT size");
}
