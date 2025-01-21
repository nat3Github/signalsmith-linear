#pragma once

#include "./stopwatch.h"
#include "plot/plot.h"

#include <complex>
#include <vector>
#include <random>

template<typename Sample>
struct RunData {
	using Complex = std::complex<Sample>;

	const size_t size;
	std::vector<std::vector<Sample>> realVectors;
	std::vector<std::vector<Sample>> positiveVectors;
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
	Sample * positive(size_t index) {
		while (index >= positiveVectors.size()) {
			std::default_random_engine engine(seed + positiveVectors.size());
			std::uniform_real_distribution<Sample> dist{0, 1};
			
			positiveVectors.emplace_back(size);
			for (auto &v : positiveVectors.back()) {
				v = dist(engine);
				while (v <= 0) v = dist(engine);
			}
		}
		return positiveVectors[index].data();
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
		return complexVectors[index].data();
	}
	
	double distance(const RunData<Sample> &other) const {
		double error2 = 0;
		
		for (size_t vi = 0; vi < realVectors.size(); ++vi) {
			auto &thisVector = realVectors[vi];
			auto &otherVector = other.realVectors[vi];
			for (size_t i = 0; i < size; ++i) {
				auto diff = thisVector[i] - otherVector[i];
				error2 += diff*diff;
			}
		}
		for (size_t vi = 0; vi < positiveVectors.size(); ++vi) {
			auto &thisVector = positiveVectors[vi];
			auto &otherVector = other.positiveVectors[vi];
			for (size_t i = 0; i < size; ++i) {
				auto diff = thisVector[i] - otherVector[i];
				error2 += diff*diff;
			}
		}
		for (size_t vi = 0; vi < complexVectors.size(); ++vi) {
			auto &thisVector = complexVectors[vi];
			auto &otherVector = other.complexVectors[vi];
			for (size_t i = 0; i < size; ++i) {
				error2 += std::norm(thisVector[i] - otherVector[i]);
			}
		}
		
		return std::sqrt(error2/size);
	}
	
	void log() const {
		for (size_t i = 0; i < realVectors.size(); ++i) {
			std::cout << "\tr" << i;
		}
		for (size_t i = 0; i < complexVectors.size(); ++i) {
			std::cout << "\tc" << i;
		}
		std::cout << "\n";
		for (size_t i = 0; i < size; ++i) {
			for (auto &vector : realVectors) {
				std::cout << "\t" << vector[i];
			}
			for (auto &vector : complexVectors) {
				std::cout << "\t" << vector[i];
			}
			std::cout << "\n";
		}
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
	
	RunPlot(const std::string &name, double benchmarkSeconds=0.05) : name(name), benchmarkSeconds(benchmarkSeconds), plotPtr(new Plot(450, 200)), plot(*plotPtr), legend(plot.legend(1.6, 1)) {
		std::cout << "\n" << name << "\n";
		for (size_t i = 0; i < name.size(); ++i) std::cout << "-";
		std::cout << "\n";
		plot.x.label(name);
	}
	RunPlot(const std::string &name, Plot &plot, double benchmarkSeconds=0.05) : name(name), benchmarkSeconds(benchmarkSeconds), plot(plot), legend(plot.legend(1.6, 1)) {
		plot.x.label(name);
	}
	~RunPlot() {
		if (benchmarkSeconds) {
			plot.y.major(0); // auto-scaled range includes 0
			plot.y.blankLabels().label("speed"); // values don't matter, only the comparison
			plot.write(name + ".svg");
		}
	}
	
	template<class PrepareAndRun>
	struct Runner {
		Runner(RunPlot &runPlot, std::string name) : runPlot(runPlot), line(runPlot.plot.line()) {
			runPlot.legend.add(line, name);
		}
		
		template<class Data>
		void run(Data data, double refTime, Data *maybeRefData, double errorLimit=1e-4) {
			PrepareAndRun obj;
			obj.prepare(data.size, data.size);
			obj.run(data);
			
			if (maybeRefData) {
				double error = data.distance(*maybeRefData);
				if (error > errorLimit) {
					std::cout << "\nsize = " << data.size << ", error = " << error << "\n";
					std::cout << "\nReference:\n";
					maybeRefData->log();
					std::cout << "\nLinear:\n";
					data.log();
					abort();
				}
			}
			
			if (runPlot.benchmarkSeconds) {
				double speed = runBenchmark(runPlot.benchmarkSeconds, [&](){
					obj.run(data);
				});
				line.add(data.size, speed*refTime);
			}
		}
	private:
		RunPlot &runPlot;
		signalsmith::plot::Line2D &line;
	};
	
	template<class PrepareAndRun>
	Runner<PrepareAndRun> runner(std::string name) {
		return {*this, name};
	}

	bool firstTick = true;
	void tick(int n, int base=0) {
		std::string nString = std::to_string(n);
		if (base > 0) {
			int log = std::round(std::log(n)/std::log(base));
			int pow = std::round(std::pow(base, log));
			if (pow == n) {
				std::string powString = std::to_string(base) + "^" + std::to_string(log);
				if (powString.size() + 1 < nString.size()) nString = powString;
			}
		}
		tick(n, nString);
	}
	void tick(int n, std::string nString) {
		if (firstTick) {
			firstTick = false;
			plot.x.major(n, nString);
		} else {
			plot.x.tick(n, nString);
		}
	}
};
