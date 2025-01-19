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

	const size_t size;
	std::vector<std::vector<Sample>> realVectors;
	std::vector<std::vector<Complex>> complexVectors;
	
	RunData(size_t size, int seed=0) : size(size), seed(seed) {}
	
	Sample * real(size_t index) {
		while (index >= realVectors.size()) {
			std::default_random_engine engine(seed + realVectors.size());
			std::uniform_real_distribution<Sample> dist{-1, 1};
			
			realVectors.emplace_back(size);
			for (auto &v : realVectors.back()) v = dist(engine);
		}
		return realVectors[index].data();
	}
	Complex * complex(size_t index) {
		while (index >= complexVectors.size()) {
			std::default_random_engine engine(seed + realVectors.size());
			std::uniform_real_distribution<Sample> dist{-1, 1};
			
			complexVectors.emplace_back(size);
			for (auto &v : complexVectors.back()) {
				v = {dist(engine), dist(engine)};
			}
		}
		return realVectors[index].data();
	}
	
	double distance(const RunData<Sample> &other) const {
		double error2 = 0;
		
		for (size_t vi = 0; vi < realVectors.size(); ++vi) {
			auto &thisVector = realVectors[vi];
			auto &otherVector = other.realVectors[vi];
			for (int i = 0; i < size; ++i) {
				auto diff = thisVector[i] - otherVector[i];
				error2 += diff*diff;
			}
		}
		for (size_t vi = 0; vi < complexVectors.size(); ++vi) {
			auto &thisVector = complexVectors[vi];
			auto &otherVector = other.complexVectors[vi];
			for (int i = 0; i < size; ++i) {
				error2 += std::norm(thisVector[i] - otherVector[i]);
			}
		}
		
		return std::sqrt(error2/size);
	}
	
private:
	int seed;
};

template<class Fn>
double runBenchmark(double benchmarkSeconds, Fn &&fn) {
	double benchmarkChunk = benchmarkSeconds/5;

	double seconds = 0;
	size_t rounds = 0, roundStep = 1;
	Stopwatch stopwatch{false};
	while (seconds < benchmarkSeconds) {
		stopwatch.start();
		
		for (size_t r = 0; r < roundStep; ++r) {
			fn();
		}

		double lap = stopwatch.seconds(stopwatch.lap());
		if (lap < benchmarkChunk) {
			roundStep *= 2;
		} else {
			seconds += lap;
			rounds += roundStep;
		}
	}
	
	return rounds/seconds;
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
};
