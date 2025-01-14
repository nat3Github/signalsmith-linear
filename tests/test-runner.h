#pragma once

#include "./stopwatch.h"
#include "plot/plot.h"

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
	
	Sample rA, rB, rC;
	Complex cA, cB, cC;
	std::vector<Sample> rvA, rvB, rvC, svAr, svAi, svBr, svBi, svCr, svCi;
	std::vector<Complex> cvA, cvB, cvC;
	Sample *rpA, *rpB, *rpC;
	Complex *cpA, *cpB, *cpC;
	signalsmith::linear::SplitComplex<Sample> sA, sB, sC;
	
	RunData(int size, int maxSize=-1, int seed=0) : size(size), maxSize(maxSize >= 0 ? maxSize : size), rvA(size), rvB(size), rvC(size), svAr(size), svAi(size), svBr(size), svBi(size), svCr(size), svCi(size), cvA(size), cvB(size), cvC(size), rpA(rvA.data()), rpB(rvB.data()), rpC(rvC.data()), cpA(cvA.data()), cpB(cvB.data()), cpC(cvC.data()), sA({svAr.data(), svAi.data()}), sB({svBr.data(), svBi.data()}), sC({svCr.data(), svCi.data()}) {
		randomise(seed);
	}
	RunData(const RunData &other) : size(other.size), maxSize(other.maxSize), rA(other.rA), rB(other.rB), rC(other.rC), cA(other.cA), cB(other.cB), cC(other.cC), rvA(other.rvA), rvB(other.rvB), rvC(other.rvC), svAr(other.svAr), svAi(other.svAi), svBr(other.svBr), svBi(other.svBi), svCr(other.svCr), svCi(other.svCi), cvA(other.cvA), cvB(other.cvB), cvC(other.cvC), rpA(rvA.data()), rpB(rvB.data()), rpC(rvC.data()), cpA(cvA.data()), cpB(cvB.data()), cpC(cvC.data()), sA({svAr.data(), svAi.data()}), sB({svBr.data(), svBi.data()}), sC({svCr.data(), svCi.data()}) {}
	
	double distance(const RunData<Sample> &other) const {
		double error2 = 0;
		
		for (int i = 0; i < size; ++i) {
			error2 += std::norm(rvA[i] - other.rvA[i]);
			error2 += std::norm(rvB[i] - other.rvB[i]);
			error2 += std::norm(rvC[i] - other.rvC[i]);
			error2 += std::norm(cvA[i] - other.cvA[i]);
			error2 += std::norm(cvB[i] - other.cvB[i]);
			error2 += std::norm(cvC[i] - other.cvC[i]);
			error2 += std::norm(sA.get(i) - other.sA.get(i));
			error2 += std::norm(sB.get(i) - other.sB.get(i));
			error2 += std::norm(sC.get(i) - other.sC.get(i));
		}
		error2 /= size;

		error2 += std::norm(rA - other.rA);
		error2 += std::norm(rB - other.rB);
		error2 += std::norm(rC - other.rC);
		error2 += std::norm(cA - other.cA);
		error2 += std::norm(cB - other.cB);
		error2 += std::norm(cC - other.cC);
		
		return std::sqrt(error2);
	}
	
private:
	void randomise(int seed) {
		std::default_random_engine randomEngine(seed);

		std::uniform_real_distribution<Sample> dist{-1, 1};
		rA = dist(randomEngine);
		rB = dist(randomEngine);
		rC = dist(randomEngine);
		for (auto &v : rvA) v = dist(randomEngine);
		for (auto &v : rvB) v = dist(randomEngine);
		for (auto &v : rvC) v = dist(randomEngine);
		cA = {dist(randomEngine), dist(randomEngine)};
		cB = {dist(randomEngine), dist(randomEngine)};
		cC = {dist(randomEngine), dist(randomEngine)};
		for (auto &v : cvA) v = {dist(randomEngine), dist(randomEngine)};
		for (auto &v : cvB) v = {dist(randomEngine), dist(randomEngine)};
		for (auto &v : cvC) v = {dist(randomEngine), dist(randomEngine)};
	}
};

double nToX(int n) {
	return (n >= 1) ? std::log(n) + 1 : n;
}

template<class Wrapper>
struct Runner {
	double benchmarkSeconds = 0.05;
	double benchmarkChunk = 0.01;

	std::string name;
	signalsmith::plot::Line2D &line;
	Stopwatch stopwatch{false};
	Wrapper wrapper;
	
	Runner(const std::string &name, signalsmith::plot::Line2D &line, signalsmith::plot::Legend &legend) : name(name), line(line) {
		legend.add(line, name);
	}
	Runner(Runner &&other) : name(other.name), line(other.line), wrapper(other.wrapper) {}

	template<class Data>
	double run(const Data &data, double refTime, const Data *refData, int strideIn, int strideOut) {
		Data copy = data;
		run(copy, refTime, refData, strideIn, strideOut);
	}

	template<class Data>
	double run(Data &data, double refTime, const Data *refData, int strideIn, int strideOut) {
		wrapper.prepare(data.size, data.maxSize);
		size_t rounds = 0, roundStep = 1;

		double error = 0;
		if (refData != nullptr) {
			wrapper.run(data, strideIn, strideOut);
			error = data.distance(*refData);
		}

		double dummySum = 1;
		double seconds = 0;
		while (seconds < benchmarkSeconds) {
			stopwatch.start();
			
			for (size_t r = 0; r < roundStep; ++r) {
				wrapper.run(data, strideIn, strideOut);
				dummySum += data.rvA[0];
			}

			double lap = stopwatch.seconds(stopwatch.lap());
			if (lap < benchmarkChunk) {
				roundStep *= 2;
			} else {
				seconds += lap;
				rounds += roundStep;
			}
		}
		double rps = rounds/seconds;
		double scaledRps = rps*refTime;
		line.add(nToX(data.size), scaledRps);

		std::cout << "\t" << data.size << "\t" << name;
		for (size_t c = name.size(); c < 23; ++c) std::cout << " ";
		std::cout << "\tspeed: " << scaledRps;
		
		if (refData != nullptr) {
			std::cout << "\terror: " << error << "\n";
			if (error > 0.0001*data.size) { // sanity check
				std::cout << "\nrA=" << data.rA << ", rB=" << data.rB << ", rC=" << data.rC << ", cA=" << data.cA << ", cB=" << data.cB << ", cC=" << data.cC << "\n";
				for (int i = 0; i < data.size && i < 10; ++i) {
					std::cout << data.rvA[i] << "\t" << data.rvB[i] << "\t" << data.rvC[i] << "\t" << data.cvA[i] << "\t" << data.cvB[i] << "\t" << data.cvC[i] << "\n";
				}

				auto &rData = *refData;
				std::cout << "\nrA=" << rData.rA << ", rB=" << rData.rB << ", rC=" << rData.rC << ", cA=" << rData.cA << ", cB=" << rData.cB << ", cC=" << rData.cC << "\n";
				for (int i = 0; i < data.size && i < 10; ++i) {
					std::cout << rData.rvA[i] << "\t" << rData.rvB[i] << "\t" << rData.rvC[i] << "\t" << rData.cvA[i] << "\t" << rData.cvB[i] << "\t" << rData.cvC[i] << "\n";
				}
				abort();
			}
		} else {
			volatile double d = dummySum;
			if (d != dummySum) {
				std::cout << "uses dummySum so the compiler can't ever remove wrapper.run()";
			}
			std::cout << "\n";
		}
		
		return scaledRps;
	}
};

struct RunPlot {
	using Plot = signalsmith::plot::Plot2D;

	std::string name;
	double benchmarkSeconds;
	std::unique_ptr<Plot> plotPtr = nullptr;
	Plot &plot;
	signalsmith::plot::Legend &legend;
	
	RunPlot(const std::string &name, double benchmarkSeconds=0.05) : name(name), benchmarkSeconds(benchmarkSeconds), plotPtr(new Plot(600, 200)), plot(*plotPtr), legend(plot.legend(0, 1)) {
		std::cout << "\n" << name << "\n";
		for (size_t i = 0; i < name.size(); ++i) std::cout << "-";
		std::cout << "\n";
		plot.x.label(name);
	}
	RunPlot(const std::string &name, Plot &plot, double benchmarkSeconds=0.05) : name(name), benchmarkSeconds(benchmarkSeconds), plot(plot), legend(plot.legend(0, 1)) {
		plot.x.label(name);
	}
	~RunPlot() {
		if (benchmarkSeconds) {
			plot.y.major(0); // auto-scaled range includes 0
			plot.y.blankLabels().label("speed"); // values don't matter, only the comparison
			plot.write(name + ".svg");
		}
	}

	bool firstTick = true;
	void tick(int n, int pow=2) {
		int log = std::round(std::log(n)/std::log(pow));
		bool isPower = n == int(std::pow(pow, log));
		std::string direct = std::to_string(n);
		std::string sci = std::to_string(pow) + "^" + std::to_string(log);
		std::string &label = (isPower && direct.size() <= sci.size()) ? direct : sci;
		if (firstTick) {
			firstTick = false;
			plot.x.major(nToX(n), label);
		} else {
			plot.x.tick(nToX(n), label);
		}
	}

	template<class Wrapper>
	Runner<Wrapper> runner(std::string name) {
		Runner<Wrapper> result{name, plot.line(), legend};
		result.benchmarkSeconds = benchmarkSeconds;
		result.benchmarkChunk = benchmarkSeconds/5;
		return result;
	}
};
