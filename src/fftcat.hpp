#ifndef SRC_FFTCAT_HPP_
#define SRC_FFTCAT_HPP_

#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include "aquila/global.h"
#include "aquila/functions.h"
#include "aquila/transform/FftFactory.h"
#include "aquila/source/FramesCollection.h"
#include "aquila/tools/TextPlot.h"

using namespace Aquila;
using std::vector;

static TextPlot plot;

void op_fft(std::shared_ptr<Fft> signalFFT, const vector<double>& source, vector<double>& target, size_t bufferSize, size_t sampleRate) {
	SignalSource in(source, sampleRate);
	FramesCollection frames(in, bufferSize);

	for (auto frame : frames) {
		SpectrumType signalSpectrum = signalFFT->fft(frame.toArray());

		for (std::size_t i = 0; i < bufferSize; ++i) {
			const double& re = signalSpectrum[i].real();
			const double& im = signalSpectrum[i].imag();
			std::cout.write(reinterpret_cast<const char*>(&re), sizeof(double));
			std::cout.write(reinterpret_cast<const char*>(&im), sizeof(double));
		}
	}
}

void op_iplot(const SpectrumType& source) {
	plot.setTitle("Input FFT");
	plot.plotSpectrum(source);
}

void op_plot(std::shared_ptr<Fft> signalFFT, const vector<double>& source, vector<double>& target, size_t bufferSize, size_t sampleRate) {
	plot.setTitle("Output FFT");
	SignalSource in(source, sampleRate);
	FramesCollection frames(in, bufferSize);

	for (auto frame : frames) {
		const SpectrumType& signalSpectrum = signalFFT->fft(frame.toArray());
		plot.plotSpectrum(signalSpectrum);
	}
}

template<typename T>
void op_ifft(std::shared_ptr<Fft> signalFFT, const SpectrumType& source, vector<double>& target, size_t bufferSize) {
	signalFFT->ifft(source, target.data());

	for (std::size_t i = 0; i < bufferSize; ++i) {
		const T& sample = target[i];
		std::cout.write(reinterpret_cast<const char*>(&sample), sizeof(T));
	}
}

template<typename T>
void ifftcat(std::vector<std::istream*> streams, size_t bufferSize, uint32_t sampleRate, bool plot) {
	SpectrumType source(bufferSize);
	vector<double> target(bufferSize);
	auto signalFFT = FftFactory::getFft(bufferSize);

	for (std::istream* is : streams) {
		double real;
		double imag;
		while (is->good()) {
			for (size_t i = 0; i < bufferSize; ++i) {
				if (is->good()) {
					is->read(reinterpret_cast<char*>(&real), sizeof(double));
					is->read(reinterpret_cast<char*>(&imag), sizeof(double));
					source[i] = { real, imag };
				} else {
					source.resize(i);
					target.resize(i);
					bufferSize = i;
				}
			}

			if(plot) {
				op_iplot(source);
			} else {
				op_ifft<T>(signalFFT, source, target, bufferSize);
			}
		}
	}
}

template<typename T>
void fftcat(std::vector<std::istream*> streams, size_t bufferSize, uint32_t sampleRate, bool plot) {
	vector<double> source(bufferSize);
	vector<double> target(bufferSize);
	auto signalFFT = FftFactory::getFft(bufferSize);

	for (std::istream* is : streams) {
		T sample;

		while (is->good()) {
			for (size_t i = 0; i < bufferSize; ++i) {
				if (is->good()) {
					is->read(reinterpret_cast<char*>(&sample), sizeof(T));
					source[i] = sample;
				} else {
					source.resize(i);
					target.resize(i);
					bufferSize = i;
				}
			}

			if(plot) {
				op_plot(signalFFT, source, target, bufferSize, sampleRate);
			} else {
				op_fft(signalFFT, source, target, bufferSize, sampleRate);
			}
		}
	}
}

#endif /* SRC_FFTCAT_HPP_ */
