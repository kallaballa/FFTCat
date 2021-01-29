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


class FFTFilter {
public:
	FrequencyType lowFreq_;
	FrequencyType highFreq_;
};

SpectrumType make_high_pass(const size_t& bufferSize, const size_t& sampleRate, const FFTFilter& filter) {
  SpectrumType filterSpectrum(bufferSize);
	for (std::size_t i = 0; i < bufferSize; ++i) {
		if (i < (bufferSize * filter.highFreq_ / sampleRate)) {
			filterSpectrum[i] = 0.0;
		} else {
			filterSpectrum[i] = 1.0;
		}
	}
	return filterSpectrum;
}

SpectrumType make_low_pass(const size_t& bufferSize, const size_t& sampleRate, const FFTFilter& filter) {
  SpectrumType filterSpectrum(bufferSize);
	for (std::size_t i = 0; i < bufferSize; ++i) {
		if (i < (bufferSize * filter.lowFreq_ / sampleRate)) {
			filterSpectrum[i] = 1.0;
		} else {
			filterSpectrum[i] = 0.0;
		}
	}
	return filterSpectrum;
}

void op_apply_filters(size_t bufferSize, uint32_t sampleRate, SpectrumType& source, const std::vector<FFTFilter>& filters) {
	for(const FFTFilter& f : filters) {
		SpectrumType highpass = make_high_pass(bufferSize, sampleRate, f);
		SpectrumType lowpass = make_low_pass(bufferSize, sampleRate, f);
		std::transform(
				std::begin(source),
				std::end(source),
				std::begin(highpass),
				std::begin(source),
				[] (Aquila::ComplexType x, Aquila::ComplexType y) { return x * y; }
		);
		std::transform(
				std::begin(source),
				std::end(source),
				std::begin(lowpass),
				std::begin(source),
				[] (Aquila::ComplexType x, Aquila::ComplexType y) { return x * y; }
		);
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

template<typename Tfftdata>
void op_fft(std::shared_ptr<Fft> signalFFT, const vector<double>& source, vector<double>& target, size_t bufferSize, size_t sampleRate, const std::vector<FFTFilter>& filters) {
	SignalSource in(source, sampleRate);
	FramesCollection frames(in, bufferSize);

	for (auto frame : frames) {
		SpectrumType signalSpectrum = signalFFT->fft(frame.toArray());
		op_apply_filters(bufferSize, sampleRate, signalSpectrum, filters);

		for (std::size_t i = 0; i < bufferSize; ++i) {
			Tfftdata re = signalSpectrum[i].real();
			Tfftdata im = signalSpectrum[i].imag();
			std::cout.write(reinterpret_cast<const char*>(&re), sizeof(Tfftdata));
			std::cout.write(reinterpret_cast<const char*>(&im), sizeof(Tfftdata));
		}
	}
}
template<typename Tsample>
void op_ifft(std::shared_ptr<Fft> signalFFT, SpectrumType& source, vector<double>& target, size_t bufferSize, size_t sampleRate, const std::vector<FFTFilter>& filters) {
	op_apply_filters(bufferSize, sampleRate, source, filters);
	signalFFT->ifft(source, target.data());

	for (std::size_t i = 0; i < bufferSize; ++i) {
		const Tsample& sample = target[i];
		std::cout.write(reinterpret_cast<const char*>(&sample), sizeof(Tsample));
	}
}

template<typename Tsample, typename Tfftdata>
void ifftcat(std::vector<std::istream*> streams, size_t bufferSize, uint32_t sampleRate, bool plot, std::vector<FFTFilter> filters) {
	SpectrumType source(bufferSize);
	vector<double> target(bufferSize);
	auto signalFFT = FftFactory::getFft(bufferSize);

	for (std::istream* is : streams) {
		Tfftdata real;
		Tfftdata imag;
		while (is->good()) {
			for (size_t i = 0; i < bufferSize; ++i) {
				if (is->good()) {
					is->read(reinterpret_cast<char*>(&real), sizeof(Tfftdata));
					is->read(reinterpret_cast<char*>(&imag), sizeof(Tfftdata));
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
				op_ifft<Tsample>(signalFFT, source, target, bufferSize, sampleRate, filters);
			}
		}
	}
}

template<typename Tsample, typename Tfftdata>
void fftcat(std::vector<std::istream*> streams, size_t bufferSize, uint32_t sampleRate, bool plot, std::vector<FFTFilter> filters) {
	vector<double> source(bufferSize);
	vector<double> target(bufferSize);
	auto signalFFT = FftFactory::getFft(bufferSize);

	for (std::istream* is : streams) {
		Tsample sample;

		while (is->good()) {
			for (size_t i = 0; i < bufferSize; ++i) {
				if (is->good()) {
					is->read(reinterpret_cast<char*>(&sample), sizeof(Tsample));
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
				op_fft<Tfftdata>(signalFFT, source, target, bufferSize, sampleRate, filters);
			}
		}
	}
}

#endif /* SRC_FFTCAT_HPP_ */
