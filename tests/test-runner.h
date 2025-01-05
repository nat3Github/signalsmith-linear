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
	std::vector<Sample> rvA, rvB, rvC;
	std::vector<Complex> cvA, cvB, cvC;
	Sample *rpA, *rpB, *rpC;
	Complex *cpA, *cpB, *cpC;
	
	RunData(int size, int maxSize, int seed=0) : size(size), maxSize(maxSize), rvA(size), rvB(size), rvC(size), cvA(size), cvB(size), cvC(size), rpA(rvA.data()), rpB(rvB.data()), rpC(rvC.data()), cpA(cvA.data()), cpB(cvB.data()), cpC(cvC.data()), randomEngine(seed) {
		randomise();
	}
	
	void randomise() {
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
	
	double distance(const RunData<Sample> &other) {
		double error2 = 0;
		
		for (int i = 0; i < size; ++i) {
			error2 += std::norm(rvA[i] - other.rvA[i]);
			error2 += std::norm(rvB[i] - other.rvB[i]);
			error2 += std::norm(rvC[i] - other.rvC[i]);
			error2 += std::norm(cvA[i] - other.cvA[i]);
			error2 += std::norm(cvB[i] - other.cvB[i]);
			error2 += std::norm(cvC[i] - other.cvC[i]);
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
	std::default_random_engine randomEngine;
};

double nToX(int n) {
	return (n >= 1) ? std::log(n) + 1 : n;
}

template<class Wrapper>
struct Runner {
	static constexpr double testSeconds = 0.05;//0.5;
	static constexpr double testChunk = 0.01;//0.1;

	std::string name;
	signalsmith::plot::Line2D &line;
	Stopwatch stopwatch{false};
	Wrapper wrapper;
	
	Runner(const std::string &name, signalsmith::plot::Line2D &line, signalsmith::plot::Legend &legend) : name(name), line(line) {
		legend.add(line, name);
	}
	Runner(Runner &&other) : name(other.name), line(other.line), wrapper(other.wrapper) {}

	template<class Data>
	void run(const Data &data, double refTime=1, const Data *refData=nullptr) {
		Data copy = data;
		run(copy, refTime, refData);
	}

	template<class Data>
	void run(Data &data, double refTime=1, const Data *refData=nullptr) {
		wrapper.prepare(data.size, data.maxSize);
		size_t rounds = 0, roundStep = 1;

		double error = 0;
		if (refData != nullptr) {
			wrapper.run(data);
			error = data.distance(*refData);
		}

		double dummySum = 1;
		double seconds = 0;
		while (seconds < testSeconds) {
			stopwatch.start();
			
			for (size_t r = 0; r < roundStep; ++r) {
				wrapper.run(data);
				dummySum += data.rvA[0];
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
		double scaledRps = rps*refTime;
		line.add(nToX(data.size), scaledRps);

		std::cout << data.size << "\t" << name;
		for (size_t c = name.size(); c < 23; ++c) std::cout << " ";
		std::cout << "\tspeed: " << scaledRps;
		
		if (refData != nullptr) {
			std::cout << "\terror: " << error << "\n";
			if (error > 0.000001*data.size) abort(); // sanity check
		} else {
			volatile double d = dummySum;
			if (d != dummySum) {
				std::cout << "uses dummySum so the compiler can't ever remove wrapper.run()";
			}
			std::cout << "\n";
		}
	}
};

struct RunPlot {
	std::string name;
	signalsmith::plot::Figure figure;
	signalsmith::plot::Plot2D &plot;
	signalsmith::plot::Legend &legend;
	
	RunPlot(const std::string &name) : name(name), plot(figure.plot(800, 250)), legend(plot.legend(0, 1)) {
		plot.x.label("size");
	}
	~RunPlot() {
		plot.y.major(0); // auto-scaled range includes 0
		plot.y.blankLabels().label("speed"); // values don't matter, only the comparison
		figure.write(name + ".svg");
	}

	bool firstTick = true;
	void tick(int n) {
		if (firstTick) {
			firstTick = false;
			plot.x.major(nToX(n), std::to_string(n));
		} else {
			plot.x.tick(nToX(n), std::to_string(n));
		}
	}

	template<class Wrapper>
	Runner<Wrapper> runner(std::string name) {
		return {name, plot.line(), legend};
	}
};
