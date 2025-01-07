#include "../fft2.h"

#include "./test-runner.h"

// ---------- wrappers

template<class Sample>
struct SimpleWrapper {
	signalsmith::fft2::SimpleFFT<Sample> fft;

	void prepare(int size, int) {
		fft.resize(size);
	}
	
	template<class Data>
	void run(Data &data) {
		fft.fft(data.size, data.cpA, data.cpB);
	}
};

template<class Sample>
struct Pow2Wrapper {
	signalsmith::fft2::Pow2FFT<Sample> fft;

	void prepare(int size, int) {
		fft.resize(size);
	}

	template<class Data>
	void run(Data& data) {
		fft.fft(data.cpA, data.cpB);
	}
};

template<class Sample>
struct SplitWrapper {
	signalsmith::fft2::SplitFFT<Sample> fft;

	void prepare(int size, int) {
		fft.resize(size);
	}

	template<class Data>
	void run(Data& data) {
		fft.fft(data.cpA, data.cpB);
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
		fft.fft(data.cpA, data.cpB);
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
		fft.fft(data.cpA, data.cpB);
	}
};

// ---------- main code

void testFfts(int maxSize, double benchmarkSeconds) {
	signalsmith::plot::Plot2D fastSizePlot(200, 200);
	auto &fastSizeLine = fastSizePlot.line();
	for (int n = 1; n < 65536; ++n) {
		int fastN = signalsmith::fft2::SplitFFT<double>::fastSizeAbove(n);
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
		double refTime = 1e-8*n*(std::log2(n) + 1); // heuristic for expected computation time, just to compare different sizes

		RunData<double> dataDouble(n, maxSize, n);
		RunData<float> dataFloat(n, maxSize, n);

		RunData<double> *refPtrDouble = nullptr;
		RunData<float> *refPtrFloat = nullptr;

		auto refDataDouble = dataDouble;
		auto refDataFloat = dataFloat;
		if (n < 256) {
			refPtrDouble = &refDataDouble;
			refPtrFloat = &refDataFloat;
			for (int f = 0; f < n; ++f) {
				std::complex<double> sumD = 0;
				std::complex<double> sumF = 0;
				for (int i = 0; i < n; ++i) {
					auto coeff = std::polar(1.0, -2*M_PI*i*f/n);
					sumD += coeff*refDataDouble.cvA[i];
					std::complex<float> vf = refDataFloat.cvA[i];
					sumF += coeff*std::complex<double>{vf.real(), vf.imag()};
				}
				refDataDouble.cvB[f] = sumD;
				refDataFloat.cvB[f] = {float(sumF.real()), float(sumF.imag())};
			}
		}

		if (pow3 + pow5 + pow7 == 0) {
			simpleDouble.run(dataDouble, refTime, refPtrDouble);
			simpleFloat.run(dataFloat, refTime, refPtrFloat);
			pow2Double.run(dataDouble, refTime, refPtrDouble);
			pow2Float.run(dataFloat, refTime, refPtrFloat);

			plot.tick(n);
		}
		if (pow5 + pow7 == 0) {
			dspDouble.run(dataDouble, refTime, refPtrDouble);
			dspFloat.run(dataFloat, refTime, refPtrFloat);
		}
		if (signalsmith::fft2::SplitFFT<double>::fastSizeAbove(n) == size_t(n)) {
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
	}
}
