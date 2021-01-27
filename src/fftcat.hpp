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

using std::vector;

template<typename T>
void cat_fft(std::vector<std::istream*> streams, size_t bufferSize, uint32_t sampleRate) {
	using namespace Aquila;
	vector<double> source(bufferSize);
	vector<double> target(bufferSize);
	auto signalFFT = FftFactory::getFft(bufferSize);

	for(std::istream* is : streams) {
		T sample;

		for(size_t i = 0; i < bufferSize; ++i) {
			if(is->good()) {
				is->read(reinterpret_cast<char*>(&sample), sizeof(T));
				source[i] = sample;
			}	else {
				return;
			}
		}

		SignalSource in(source, sampleRate);
		FramesCollection frames(in,bufferSize);

		for (auto frame : frames) {
			SpectrumType signalSpectrum = signalFFT->fft(frame.toArray());

			for (std::size_t i = 0; i < bufferSize; ++i) {
				for (size_t j = 0; j < bufferSize; ++j) {
					const T& mag = abs(signalSpectrum[j]) / bufferSize;
					std::cout.write(reinterpret_cast<const char*>(&mag), sizeof(T));
				}
			}
		}
	}
}

template<typename T>
void plot_fft(std::vector<std::istream*> streams, size_t bufferSize, uint32_t sampleRate) {
	using namespace Aquila;
	vector<double> source(bufferSize);
	vector<double> target(bufferSize);
	auto signalFFT = FftFactory::getFft(bufferSize);
	TextPlot plot;
	plot.setTitle("Signal");

	for(std::istream* is : streams) {

		T sample;
		while (is->good()) {
			for (size_t i = 0; i < bufferSize; ++i) {
				if (is->good()) {
					is->read(reinterpret_cast<char*>(&sample), sizeof(T));
					source[i] = sample;
				} else {
					source.resize(i);
				}
			}

			SignalSource in(source, sampleRate);
			FramesCollection frames(in, bufferSize);
			for (auto frame : frames) {
				const SpectrumType& signalSpectrum = signalFFT->fft(frame.toArray());
				plot.plotSpectrum(signalSpectrum);
			}
		}
	}
}



#endif /* SRC_FFTCAT_HPP_ */
