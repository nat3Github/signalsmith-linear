#include "../stft.h"

#include "./test-runner.h"

#include "plot/plot.h"

template<class Sample>
void testStft(size_t channels, size_t blockSize, size_t minInterval, size_t maxInterval) {
	signalsmith::plot::Figure figure;

	double debugTime = 0;
	auto debugTick = [&](){
		double t = debugTime;
		debugTime += 0.1;
		return t;
	};

	std::default_random_engine randomEngine(24680);
	std::uniform_int_distribution<size_t> intervalDist{minInterval, maxInterval};
	
	size_t length = blockSize*10;
	RunData<Sample> data(length);
	
	std::vector<Sample *> input, output;
	for (size_t c = 0; c < channels; ++c) {
		input.push_back(data.real(c));
		for (size_t i = 0; i < length; ++i) {
			input.back()[i] += 2*std::sin(i*2.0/blockSize);
		}
		output.push_back(data.real(c + channels));
	}
	
	signalsmith::linear::DynamicSTFT<Sample> stft;
	stft.configure(channels, channels, blockSize);
	
	const size_t plotChannel = 0;
	
	auto &resultPlot = figure(0, -1).plot(800, 200);
	resultPlot.y.major(0);
	auto &outputLine = resultPlot.line();
	auto &inputLine = resultPlot.line();
	resultPlot.legend(1, 1).add(outputLine, "output").add(inputLine, "input");
	for (size_t i = 0; i < length; ++i) {
		inputLine.add(i, input[plotChannel][i]);
	}

//	auto &inputPlot = figure(0, 0)(0, 0).plot(200, 200).title("input");
//	inputPlot.y.major(0);
//	auto &historyLine = inputPlot.line();
//	auto &inputTimeLine = inputPlot.line();

//	auto &outputPlot = figure(0, 0)(1, 0).plot(200, 200).title("output");
//	outputPlot.y.major(0);
//	auto &sumLine = outputPlot.line();
//	auto &outputTimeLine = outputPlot.line();
//	auto &windowProductLine = outputPlot.line();
//	outputPlot.legend(2, 1).add(sumLine, "rolling buffer").add(outputTimeLine, "current block (centred)").add(windowProductLine, "window product");

	auto &windowPlot = figure(0, 0)(0, 1).plot(200, 200).title("window");
	windowPlot.y.major(0);
	auto &aWindowLine = windowPlot.line();
	auto &sWindowLine = windowPlot.line();
	windowPlot.legend(-1, 1).add(aWindowLine, "analysis").add(sWindowLine, "synthesis");

	auto &spectrumPlot = figure(0, 0)(1, 1).plot(200, 200).title("spectrum");
	spectrumPlot.y.major(0);
	auto &realLine = spectrumPlot.line();
	auto &imagLine = spectrumPlot.line();
	
	size_t start = 0;
	while (start + blockSize < length) {
		size_t interval = intervalDist(randomEngine);
		
		// Randomly mess the windows up to make sure we still have perfect reconstruction
		{
			std::uniform_real_distribution<double> dist{0, 1};
			Sample factor = 0.5 + dist(randomEngine);
			size_t windowModWidth = dist(randomEngine)*blockSize/2;
			size_t windowModPos = dist(randomEngine)*(blockSize - windowModWidth);
			for (size_t i = 0; i < windowModWidth; i++) {
				stft.analysisWindow()[i + windowModPos] *= factor;
			}
		}
		{
			std::uniform_real_distribution<double> dist{0, 1};
			Sample factor = 0.5 + dist(randomEngine);
			size_t windowModWidth = dist(randomEngine)*blockSize/2;
			size_t windowModPos = dist(randomEngine)*(blockSize - windowModWidth);
			for (size_t i = 0; i < windowModWidth; i++) {
				stft.synthesisWindow()[i + windowModPos] *= factor;
			}
		}
		for (size_t i = 0; i < blockSize; ++i) {
			aWindowLine.add(i, stft.analysisWindow()[i]);
			sWindowLine.add(i, stft.synthesisWindow()[i]);
		}
		windowPlot.toFrame(debugTick());
		
		for (size_t c = 0; c < channels; ++c) {
			stft.readOutput(c, 0, interval, output[c] + start);
		}
//		for (size_t i = 0; i < blockSize; ++i) {
//			sumLine.add(i, stft.sumBuffer[i]);
//			windowProductLine.add(i, stft.sumWindowProducts[i]);
//		}
//		for (size_t i = 0; i < stft.timeBuffer.size(); ++i) {
//			outputTimeLine.add(i, stft.timeBuffer[i]);
//		}
//		sumLine.marker(stft.outputPos, 0);
//		outputPlot.toFrame(debugTick());

		for (size_t c = 0; c < channels; ++c) {
			stft.writeInput(c, 0, interval, input[c] + start);
		}
		stft.moveInput(interval);

		stft.analyse();
		
//		for (size_t i = 0; i < blockSize; ++i) {
//			historyLine.add(i, stft.inputBuffer[i]);
//		}
//		historyLine.marker(stft.inputPos, 0);
//		for (size_t i = 0; i < stft.timeBuffer.size(); ++i) {
//			inputTimeLine.add(i, stft.timeBuffer[i]);
//		}
//		inputPlot.toFrame(debugTick());
		
		for (size_t f = 0; f < stft.bands(); ++f) {
			realLine.add(f, stft.spectrum(0)[f].real());
			imagLine.add(f, stft.spectrum(0)[f].imag());
		}
		spectrumPlot.toFrame(debugTick());
		
		stft.synthesise(interval);

//		for (size_t i = 0; i < blockSize; ++i) {
//			sumLine.add(i, stft.sumBuffer[i]);
//			windowProductLine.add(i, stft.sumWindowProducts[i]);
//		}
//		for (size_t i = 0; i < stft.timeBuffer.size(); ++i) {
//			outputTimeLine.add(i, stft.timeBuffer[i]);
//		}
//		sumLine.marker(stft.outputPos, 0);
//		outputPlot.toFrame(debugTick());

		start += interval;

		for (size_t i = 0; i < start; ++i) {
			outputLine.add(i, output[plotChannel][i]);
		}
		outputLine.toFrame(debugTick());
	}
	
	double error = 0;
	for (size_t i = 0; i < start - blockSize; ++i) {
		for (size_t c = 0; c < channels; ++c) {
			error += std::abs(output[c][i + blockSize] - input[c][i]);
		}
	}
	figure.loopFrame(debugTick());
	figure.write("stft-debug.svg");

	if (error > length*0.001) {
		LOG_EXPR(error);
		abort();
	}
}

void testStfts(int, double) {
	std::cout << "STFT\n----\n";
	testStft<double>(2, 256, 32, 192);
}
