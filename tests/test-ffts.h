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
	void normalise(Data &data) {
		data.complexToSplit(1);
		data.complexToSplit(2);
	}
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
	void normalise(Data &data) {
		data.complexToSplit(1);
		data.complexToSplit(2);
	}
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
	signalsmith::linear::SplitFFT<Sample, true> fft;

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
	void normalise(Data &data) {
		data.complexToSplit(1);
		data.complexToSplit(2);
	}
};
template<class Sample>
struct SplitWrapperSplit {
	signalsmith::linear::SplitFFT<Sample, true> fft;

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
struct UnSplitWrapper {
	signalsmith::linear::SplitFFT<Sample, false> fft;

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
	void normalise(Data &data) {
		data.complexToSplit(1);
		data.complexToSplit(2);
	}
};
template<class Sample>
struct UnSplitWrapperSplit {
	signalsmith::linear::SplitFFT<Sample, false> fft;

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
	void normalise(Data &data) {
		data.complexToSplit(1);
		data.complexToSplit(2);
	}
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
	void normalise(Data &data) {
		data.complexToSplit(1);
		data.complexToSplit(2);
	}
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

	auto splitFloat = plot.runner<SplitWrapper<float>>("Split (float)");
	auto splitDouble = plot.runner<SplitWrapper<double>>("Split (double)");
	auto splitFloatSplit = plot.runner<SplitWrapperSplit<float>>("Split (split float)");
	auto splitDoubleSplit = plot.runner<SplitWrapperSplit<double>>("Split (split double)");
#ifdef INCLUDE_UNSPLIT_SPLITFFT
	auto unsplitFloat = plot.runner<UnSplitWrapper<float>>("Un-Split (float)");
	auto unsplitDouble = plot.runner<UnSplitWrapper<double>>("Un-Split (double)");
	auto unsplitFloatSplit = plot.runner<UnSplitWrapperSplit<float>>("Un-Split (split float)");
	auto unsplitDoubleSplit = plot.runner<UnSplitWrapperSplit<double>>("Un-Split (split double)");
#endif
	auto pow2Float = plot.runner<Pow2Wrapper<float>>("Pow2 (float)");
	auto pow2Double = plot.runner<Pow2Wrapper<double>>("Pow2 (double)");
	auto pow2FloatSplit = plot.runner<Pow2WrapperSplit<float>>("Pow2 (split float)");
	auto pow2DoubleSplit = plot.runner<Pow2WrapperSplit<double>>("Pow2 (split double)");

	auto simpleFloat = plot.runner<SimpleWrapper<float>>("Simple (float)");
	auto simpleDouble = plot.runner<SimpleWrapper<double>>("Simple (double)");
	auto simpleFloatSplit = plot.runner<SimpleWrapperSplit<float>>("Simple (split float)");
	auto simpleDoubleSplit = plot.runner<SimpleWrapperSplit<double>>("Simple (split double)");
	auto dspFloat = plot.runner<SignalsmithDSPWrapper<float>>("DSP library (float)");
	auto dspDouble = plot.runner<SignalsmithDSPWrapper<double>>("DSP library (double)");

	auto runSize = [&](int n, int pow3=0, int pow5=0, int pow7=0, bool alwaysRunSplit=false){
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
		checkResults = (n < 48); // below this, it doesn't use the broken SIMD operations
#endif

		if (checkResults && n < 256) {
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
		if (alwaysRunSplit || signalsmith::linear::SplitFFT<double>::fastSizeAbove(n) == size_t(n)) {
			std::cout << "split double           \r";
			splitDouble.run(dataDouble, refTime, refPtrDouble);
			std::cout << "split float           \r";
			splitFloat.run(dataFloat, refTime, refPtrFloat);
			std::cout << "split double (split)           \r";
			splitDoubleSplit.run(dataDouble, refTime, refPtrDouble);
			std::cout << "split float (split)           \r";
			splitFloatSplit.run(dataFloat, refTime, refPtrFloat);
#ifdef INCLUDE_UNSPLIT_SPLITFFT
			std::cout << "unsplit double           \r";
			unsplitDouble.run(dataDouble, refTime, refPtrDouble);
			std::cout << "unsplit float           \r";
			unsplitFloat.run(dataFloat, refTime, refPtrFloat);
			std::cout << "unsplit double (split)           \r";
			unsplitDoubleSplit.run(dataDouble, refTime, refPtrDouble);
			std::cout << "unsplit float (split)           \r";
			unsplitFloatSplit.run(dataFloat, refTime, refPtrFloat);
#endif
		}
	};

	if (benchmarkSeconds == 0) {
		// test the split FFT for 1-64, but it's not expected to be fast so don't benchmark it
		for (int n = 1; n < 64; ++n) {
			runSize(n, -1, -1, -1, true);
		}
	}
	
	for (int n = 1; n <= maxSize; n *= 2) {
		if (n/16) runSize(n*9/16, 2);
		if (n/8) runSize(n*5/8, 0, 1);
		if (n/4) runSize(n*3/4, 1);
		if (n/8) runSize(n*7/8, 0, 0, 1);
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

template<typename Sample>
struct SplitRunner {
	using Complex = std::complex<Sample>;

	signalsmith::linear::SplitFFT<Sample, true> fft;
	
	std::vector<Sample> timeR, timeI, freqR, freqI, time2R, time2I;
	std::vector<Complex> time, freq, time2;
	
	void run(size_t n, signalsmith::plot::Line2D &line, double benchmarkSeconds) {
		time.resize(n);
		freq.resize(n);
		time2.resize(n);

		Stopwatch stopwatch;
		fft.resize(n);

		double totalTime = 0;
		std::vector<double> stepTimes(fft.steps()*2);
		while (totalTime < benchmarkSeconds) {
			for (size_t s = 0; s < fft.steps(); ++s) {
				stopwatch.start();
				for (size_t r = 0; r < 10; ++r) {
					fft.fft(s, time.data(), freq.data());
				}
				double t = stopwatch.lap();
				totalTime += t;
				stepTimes[s] += t;
			}
			for (size_t s = 0; s < fft.steps(); ++s) {
				stopwatch.start();
				for (size_t r = 0; r < 10; ++r) {
					fft.ifft(s, freq.data(), time2.data());
				}
				double t = stopwatch.lap();
				totalTime += t;
				stepTimes[s + fft.steps()] += t;
			}
		}
		Sample averageTime = totalTime/stepTimes.size();
		for (size_t s = 0; s < stepTimes.size(); ++s) {
			line.add(double(s)/stepTimes.size(), stepTimes[s]/averageTime);
			line.add(double(s + 1)/stepTimes.size(), stepTimes[s]/averageTime);
		}
	}

	void runSplit(size_t n, signalsmith::plot::Line2D &line, double benchmarkSeconds) {
		timeR.resize(n);
		timeI.resize(n);
		freqR.resize(n);
		freqI.resize(n);
		time2R.resize(n);
		time2I.resize(n);

		Stopwatch stopwatch;
		fft.resize(n);

		double totalTime = 0;
		std::vector<double> stepTimes(fft.steps()*2);
		while (totalTime < benchmarkSeconds) {
			for (size_t s = 0; s < fft.steps(); ++s) {
				stopwatch.start();
				for (size_t r = 0; r < 10; ++r) {
					fft.fft(s, timeR.data(), timeI.data(), freqR.data(), freqI.data());
				}
				double t = stopwatch.lap();
				totalTime += t;
				stepTimes[s] += t;
			}
			for (size_t s = 0; s < fft.steps(); ++s) {
				stopwatch.start();
				for (size_t r = 0; r < 10; ++r) {
					fft.ifft(s, freqR.data(), freqI.data(), time2R.data(), time2I.data());
				}
				double t = stopwatch.lap();
				totalTime += t;
				stepTimes[s + fft.steps()] += t;
			}
		}
		Sample averageTime = totalTime/stepTimes.size();
		for (size_t s = 0; s < stepTimes.size(); ++s) {
			line.add(double(s)/stepTimes.size(), stepTimes[s]/averageTime);
			line.add(double(s + 1)/stepTimes.size(), stepTimes[s]/averageTime);
		}
	}
};

void testComplexFftSplits(size_t maxSize, double benchmarkSeconds) {
	std::cout << "\nfft splits\n----------\n";
	benchmarkSeconds *= 1000;

	signalsmith::plot::Figure figure;
	auto &floatPlot = figure(0, 0).plot(250, 200).title("float");
	floatPlot.y.major(0).label("time").minor(1, "100%").minor(2, "200%");
	floatPlot.x.blank();
	auto &floatSplitPlot = figure(1, 0).plot(250, 200).title("float (split-complex)");
	floatSplitPlot.x.copyFrom(floatPlot.x);
	floatSplitPlot.y.copyFrom(floatPlot.y).flip().label("");
	auto &doublePlot = figure(0, 1).plot(250, 200).title("double");
	doublePlot.x.copyFrom(floatPlot.x);
	doublePlot.y.copyFrom(floatPlot.y);
	auto &doubleSplitPlot = figure(1, 1).plot(250, 200).title("double (split-complex)");
	doubleSplitPlot.x.copyFrom(floatPlot.x);
	doubleSplitPlot.y.copyFrom(floatPlot.y).flip().label("");

	SplitRunner<float> floatRunner;
	SplitRunner<double> doubleRunner;
	
	auto runSize = [&](size_t n){
		std::cout << "\tn = " << n << "      \r" << std::flush;
		{
			auto &line = floatPlot.fill().fillToY(0);
			floatRunner.run(n, line, benchmarkSeconds);
		}
		{
			auto &line = floatSplitPlot.fill().fillToY(0);
			floatRunner.runSplit(n, line, benchmarkSeconds);
		}
		{
			auto &line = doublePlot.fill().fillToY(0);
			doubleRunner.run(n, line, benchmarkSeconds);
		}
		{
			auto &line = doubleSplitPlot.fill().fillToY(0);
			doubleRunner.runSplit(n, line, benchmarkSeconds);
		}
	};

	runSize(256);
	for (size_t m = 512; m <= maxSize; m *= 2) {
		runSize(m/4*3);
		runSize(m);
		runSize(m/4*5);
	}
	std::cout << "\n";
	
	figure.write("fft-splits.svg");
}

template<class Sample, bool modified>
void testRealFft(int size) {
	signalsmith::linear::FFT<Sample> complexFft(size);
	signalsmith::linear::RealFFT<Sample, true, modified> realFft(size);
	
	RunData<Sample> data(size);
	auto *realTime = data.real(0);
	auto *realFreqC = data.complex(0);
	auto realFreqS = data.split(0);
	auto *realTime2C = data.real(1);
	auto realTime2S = data.real(2);
	
	auto *complexTime = data.complex(1);
	auto *complexFreq = data.complex(2);
	auto *complexTime2 = data.complex(3);
	
	for (int i = 0; i < size; ++i) {
		complexTime[i] = realTime[i];
		// So we can see in the debug output whether they've changde
		realFreqC[i] += 200;
		realFreqS.real[i] += 300;
		realFreqS.imag[i] += 400;
	}
	
	if (modified) {
		// Half-bin shift downwards
		for (int i = 0; i < size; ++i) {
			Sample shiftPhase = -M_PI*i/size;
			complexTime[i] *= std::polar(Sample(1), shiftPhase);
		}
	}

	// Interleaved complex input
	realFft.fft(realTime, realFreqC);
	realFft.ifft(realFreqC, realTime2C);
	// Split-complex input
	realFft.fft(realTime, realFreqS.real, realFreqS.imag);
	realFft.ifft(realFreqS.real, realFreqS.imag, realTime2S);
	
	complexFft.fft(complexTime, complexFreq);
	complexFft.ifft(complexFreq, complexTime2);

	if (modified) {
		// Half-bin shift upwards
		for (int i = 0; i < size; ++i) {
			Sample shiftPhase = M_PI*i/size;
			complexTime2[i] *= std::polar(Sample(1), shiftPhase);
		}
	} else {
		complexFreq[0].imag(complexFreq[size/2].real()); // Nyquist packing
	}

	double error = 0;
	for (int i = 0; i < size/2; ++i) {
		error += std::norm(complexFreq[i] - realFreqC[i]);
		error += std::norm(complexFreq[i] - realFreqS[i]);
	}
	for (int i = 0; i < size; ++i) {
		error += std::norm(complexTime2[i] - realTime2C[i]);
		error += std::norm(complexTime2[i] - realTime2S[i]);
	}
	error = std::sqrt(error/size);
	
	if (error >= size*0.001) {
		LOG_EXPR(size);
		LOG_EXPR(modified);
		LOG_EXPR(error);
		data.log();
		abort();
	}
}

void testRealFfts(int, double) {
	std::cout << "Real FFTs\n---------\n";

	int sizeEnd = 128;
#if defined(__FAST_MATH__) && (__apple_build_version__ >= 16000000) && (__apple_build_version__ <= 16000099)
		// Apple Clang 16.0.0 is broken, and generates *completely* incorrect code for some SIMD operations
		sizeEnd = 16; // below this, it doesn't use the broken SIMD operations
#endif

	for (int size = 2; size < sizeEnd; size += 2) {
		testRealFft<float, false>(size);
		testRealFft<double, false>(size);
		testRealFft<float, true>(size);
		testRealFft<double, true>(size);
	}
}

void testFfts(int maxSize, double benchmarkSeconds) {
	testComplexFfts(maxSize, benchmarkSeconds);
	if (benchmarkSeconds > 0) {
		testComplexFftSplits(maxSize, benchmarkSeconds);
	}
	
	testRealFfts(maxSize, benchmarkSeconds);
}
