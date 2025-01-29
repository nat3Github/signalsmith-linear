#include "../fft.h"

#include "./test-runner.h"

// ---------- wrappers

template<class Sample>
struct SimpleWrapper {
	signalsmith::linear::SimpleFFT<Sample> fft;

	void prepare(size_t size, size_t) {
		fft.resize(size);
	}
	
	template<class Data>
	void run(Data &data) {
		auto *time = data.complex(0), *freq = data.complex(1), *time2 = data.complex(2);
		fft.fft(data.size, time, freq);
		fft.ifft(data.size, freq, time2);
	}
	
	template<class Data>
	void normalise(Data &) {}
};

template<class Sample>
struct SimpleWrapperSplit {
	signalsmith::linear::SimpleFFT<Sample> fft;

	void prepare(size_t size, size_t) {
		fft.resize(size);
	}
	
	template<class Data>
	void run(Data &data) {
		auto time = data.split(0), freq = data.split(1), time2 = data.split(2);
		fft.fft(data.size, time.real, time.imag, freq.real, freq.imag);
		fft.ifft(data.size, freq.real, freq.imag, time2.real, time2.imag);
	}
	
	template<class Data>
	void normalise(Data &data) {
		data.splitToComplex(1);
		data.splitToComplex(2);
	}
};

template<class Sample>
struct Pow2Wrapper {
	signalsmith::linear::Pow2FFT<Sample> fft;

	void prepare(size_t size, size_t) {
		fft.resize(size);
	}

	template<class Data>
	void run(Data& data) {
		auto *time = data.complex(0), *freq = data.complex(1), *time2 = data.complex(2);
		fft.fft(time, freq);
		fft.ifft(freq, time2);
	}

	template<class Data>
	void normalise(Data &) {}
};
template<class Sample>
struct Pow2WrapperSplit {
	signalsmith::linear::Pow2FFT<Sample> fft;

	void prepare(size_t size, size_t) {
		fft.resize(size);
	}

	template<class Data>
	void run(Data &data) {
		auto time = data.split(0), freq = data.split(1), time2 = data.split(2);
		fft.fft(time.real, time.imag, freq.real, freq.imag);
		fft.ifft(freq.real, freq.imag, time2.real, time2.imag);
	}
	
	template<class Data>
	void normalise(Data &data) {
		data.splitToComplex(1);
		data.splitToComplex(2);
	}
};

template<class Sample>
struct SplitWrapper {
	signalsmith::linear::SplitFFT<Sample> fft;

	void prepare(size_t size, size_t) {
		fft.resize(size);
	}

	template<class Data>
	void run(Data& data) {
		auto *time = data.complex(0), *freq = data.complex(1), *time2 = data.complex(2);
		fft.fft(time, freq);
		fft.ifft(freq, time2);
	}

	template<class Data>
	void normalise(Data &) {}
};
template<class Sample>
struct SplitWrapperSplit {
	signalsmith::linear::SplitFFT<Sample> fft;

	void prepare(size_t size, size_t) {
		fft.resize(size);
	}

	template<class Data>
	void run(Data &data) {
		auto time = data.split(0), freq = data.split(1), time2 = data.split(2);
		fft.fft(time.real, time.imag, freq.real, freq.imag);
		fft.ifft(freq.real, freq.imag, time2.real, time2.imag);
	}
	
	template<class Data>
	void normalise(Data &data) {
		data.splitToComplex(1);
		data.splitToComplex(2);
	}
};

#include "./others/signalsmith-fft.h"
template<class Sample>
struct SignalsmithFFTWrapper {
	signalsmith::FFT<Sample> fft{1};

	void prepare(size_t size, size_t) {
		fft.setSize(size);
	}
	
	template<class Data>
	void run(Data& data) {
		auto *time = data.complex(0), *freq = data.complex(1), *time2 = data.complex(2);
		fft.fft(time, freq);
		fft.ifft(freq, time2);
	}

	template<class Data>
	void normalise(Data &) {}
};

#include "./others/dsp/fft.h"
template<class Sample>
struct SignalsmithDSPWrapper {
	signalsmith::fft::FFT<Sample> fft{1};

	void prepare(size_t size, size_t) {
		fft.setSize(size);
	}
	
	template<class Data>
	void run(Data& data) {
		auto *time = data.complex(0), *freq = data.complex(1), *time2 = data.complex(2);
		fft.fft(time, freq);
		fft.ifft(freq, time2);
	}

	template<class Data>
	void normalise(Data &) {}
};

// ---------- main code

void testComplexFfts(int maxSize, double benchmarkSeconds) {
	if (benchmarkSeconds > 0) {
		signalsmith::plot::Plot2D fastSizePlot(200, 200);
		auto &fastSizeLine = fastSizePlot.line();
		for (int n = 1; n < 65536; ++n) {
			size_t fastN = signalsmith::linear::SplitFFT<double>::fastSizeAbove(n);
			fastSizeLine.add(std::log2(n), std::log2(fastN));
		}
		fastSizePlot.line().add(0, 0).add(16, 16);
		fastSizePlot.line().add(0, std::log2(1.25)).add(16 - std::log2(1.25), 16);
		fastSizePlot.write("fft-fast-sizes.svg");
	}
	
	RunPlot plot("ffts", benchmarkSeconds);

	auto simpleDouble = plot.runner<SimpleWrapper<double>>("Simple (double)");
	auto simpleFloat = plot.runner<SimpleWrapper<float>>("Simple (float)");
	auto simpleDoubleSplit = plot.runner<SimpleWrapperSplit<double>>("Simple (split double)");
	auto simpleFloatSplit = plot.runner<SimpleWrapperSplit<float>>("Simple (split float)");
	auto pow2Double = plot.runner<Pow2Wrapper<double>>("Pow2 (double)");
	auto pow2Float = plot.runner<Pow2Wrapper<float>>("Pow2 (float)");
	auto pow2DoubleSplit = plot.runner<Pow2WrapperSplit<double>>("Pow2 (split double)");
	auto pow2FloatSplit = plot.runner<Pow2WrapperSplit<float>>("Pow2 (split float)");
	auto splitDouble = plot.runner<SplitWrapper<double>>("Split (double)");
	auto splitFloat = plot.runner<SplitWrapper<float>>("Split (float)");
	auto splitDoubleSplit = plot.runner<SplitWrapperSplit<double>>("Split (split double)");
	auto splitFloatSplit = plot.runner<SplitWrapperSplit<float>>("Split (split float)");
	auto dspDouble = plot.runner<SignalsmithDSPWrapper<double>>("DSP library (double)");
	auto dspFloat = plot.runner<SignalsmithDSPWrapper<float>>("DSP library (float)");

	auto runSize = [&](int n, int pow3=0, int pow5=0, int pow7=0){
		std::cout << "                n = " << n << "\r" << std::flush;
		double refTime = 1e-8*n*(std::log2(n) + 1); // heuristic for expected computation time, just to compare different sizes

		RunData<double> dataDouble(n);
		RunData<float> dataFloat(n);
		for (size_t i = 0; i < 3; ++i) {
			dataDouble.complexToSplit(i);
			dataFloat.complexToSplit(i);
		}

		RunData<double> *refPtrDouble = nullptr;
		RunData<float> *refPtrFloat = nullptr;

		auto refDataDouble = dataDouble;
		auto refDataFloat = dataFloat;
		bool checkResults = true;
#if defined(__FAST_MATH__) && (__apple_build_version__ >= 16000000) && (__apple_build_version__ <= 16000099)
		// Apple Clang 16.0.0 is broken, and generates *completely* incorrect code for some SIMD operations
		checkResults = false;
#endif

		if (n < 256 && checkResults) {
			refPtrDouble = &refDataDouble;
			refPtrFloat = &refDataFloat;
			auto *timeDouble = refDataDouble.complex(0);
			auto *freqDouble = refDataDouble.complex(1);
			auto *time2Double = refDataDouble.complex(2);
			auto *timeFloat = refDataFloat.complex(0);
			auto *freqFloat = refDataFloat.complex(1);
			auto *time2Float = refDataFloat.complex(2);
			for (int f = 0; f < n; ++f) {
				std::complex<double> sumD = 0;
				std::complex<double> sumF = 0;
				for (int i = 0; i < n; ++i) {
					auto coeff = std::polar(1.0, -2*M_PI*i*f/n);
					sumD += coeff*timeDouble[i];
					std::complex<float> vf = timeFloat[i];
					sumF += coeff*std::complex<double>{vf.real(), vf.imag()};
				}
				freqDouble[f] = sumD;
				freqFloat[f] = {float(sumF.real()), float(sumF.imag())};
				// Round-trip should match, with scale factor N
				time2Double[f] = timeDouble[f]*double(n);
				time2Float[f] = timeFloat[f]*float(n);
			}
			for (size_t i = 0; i < 3; ++i) {
				refDataDouble.complexToSplit(i);
				refDataFloat.complexToSplit(i);
			}
		}

		if (pow3 + pow5 + pow7 == 0) {
			std::cout << "simple double           \r";
			simpleDouble.run(dataDouble, refTime, refPtrDouble);
			std::cout << "simple double (split)           \r";
			simpleDoubleSplit.run(dataDouble, refTime, refPtrDouble);
			std::cout << "simple float           \r";
			simpleFloat.run(dataFloat, refTime, refPtrFloat);
			std::cout << "simple float (split)           \r";
			simpleFloatSplit.run(dataFloat, refTime, refPtrFloat);
			std::cout << "Pow2 double           \r";
			pow2Double.run(dataDouble, refTime, refPtrDouble);
			std::cout << "Pow2 float           \r";
			pow2Float.run(dataFloat, refTime, refPtrFloat);
			std::cout << "Pow2 double (split)           \r";
			pow2DoubleSplit.run(dataDouble, refTime, refPtrDouble);
			std::cout << "Pow2 float (split)           \r";
			pow2FloatSplit.run(dataFloat, refTime, refPtrFloat);
		}
		if (pow5 + pow7 == 0) {
			std::cout << "DSP double           \r";
			dspDouble.run(dataDouble, refTime, refPtrDouble);
			std::cout << "DSP float           \r";
			dspFloat.run(dataFloat, refTime, refPtrFloat);
		}
		if (signalsmith::linear::SplitFFT<double>::fastSizeAbove(n) == size_t(n)) {
			std::cout << "split double           \r";
			splitDouble.run(dataDouble, refTime, refPtrDouble);
			std::cout << "split float           \r";
			splitFloat.run(dataFloat, refTime, refPtrFloat);
			std::cout << "split double (split)           \r";
			splitDoubleSplit.run(dataDouble, refTime, refPtrDouble);
			std::cout << "split float (split)           \r";
			splitFloatSplit.run(dataFloat, refTime, refPtrFloat);
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
	std::cout << "                                         \r" << std::flush;

	plot.plot.x.range(std::log, 1, maxSize).label("FFT size");
}

void testFfts(int maxSize, double benchmarkSeconds) {
	testComplexFfts(maxSize, benchmarkSeconds);
}
