#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <iostream>
#include <vector>
#include <boost/program_options.hpp>
#include <fstream>

#include "fftcat.hpp"

namespace po = boost::program_options;

using std::string;
using std::cerr;
using std::endl;
using std::vector;

//does report true for 0
bool is_power_of_two(const size_t& x) {
	return (x & (x - 1)) == 0;
}

int main(int argc, char** argv) {
  std::vector<string> files;
  size_t fftSize = 256;
  uint32_t sampleRate = 44100;
  uint32_t sampleLen = 4;
	FrequencyType lowPassFreq = std::numeric_limits<FrequencyType>::max();
	FrequencyType highPassFreq = 0;

  po::options_description genericDesc("Options");
  genericDesc.add_options()("help,h", "Produce help message")
		("inverse,i", "Calculate the ifft instead of the fft")
		("plot,p", "Write a signal plot to the terminal instead of the raw fft stream")
		("fftsize,z", po::value<size_t>(&fftSize)->default_value(fftSize),"The number of bins for the fft. must be a power of 2")
		("samplerate,r", po::value<uint32_t>(&sampleRate)->default_value(sampleRate),"The i/o sample rate in hertz")
  	("samplelen,s", po::value<uint32_t>(&sampleLen)->default_value(sampleLen),"The i/o sample length in bytes")
		("lowpass,l", po::value<FrequencyType>(&lowPassFreq)->default_value(lowPassFreq),"The frequency of the low pass filter in hertz")
		("highpass,k", po::value<FrequencyType>(&highPassFreq)->default_value(highPassFreq),"The frequency of the high pass filter in hertz");

	po::options_description hidden("Hidden options");
  hidden.add_options()("files", po::value<std::vector<string>>(&files), "files");

  po::options_description cmdline_options;
  cmdline_options.add(genericDesc).add(hidden);

  po::positional_options_description p;
  p.add("files", -1);

  po::options_description visible;
  visible.add(genericDesc);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cerr << "Usage: fftcat [options] [FILES...]" << std::endl;
    std::cerr << visible;
    return 1;
  }

  std::vector<std::istream*> streams;
  if(files.empty() || (files.size() == 1 && files[0] == "-")) {
  	streams.push_back(&std::cin);
  } else {
  	for(const string& f : files) {
  		streams.push_back(new std::ifstream(f));
  		if(!streams.back()->good()) {
  			std::cerr << "Can't open file: " << f << std::endl;
  			return 3;
  		}
  	}
  }

  if(sampleLen == 0 || !is_power_of_two(sampleLen)) {
  	std::cerr << "sampleLen has to be a power of 2" << std::endl;
  	return 2;
  }

  if(sampleLen > 4) {
  	std::cerr << "sampleLen has to be <= 4" << std::endl;
  	return 2;
  }

  if(fftSize == 0 || !is_power_of_two(fftSize)) {
  	std::cerr << "fftsize has to be a power of 2" << std::endl;
  	return 2;
  }

  bool doPlot = vm.count("plot");

	std::vector<FFTFilter> filters;
	if(lowPassFreq < std::numeric_limits<FrequencyType>::max()) {
		filters.push_back({LOW_PASS, lowPassFreq});
	}

	if(highPassFreq > 0) {
		filters.push_back({HIGH_PASS, highPassFreq});
	}

	//TODO type-selection based on the platform.
  if(vm.count("inverse")) {
		if(sampleLen == 1) {
			ifftcat<uint8_t, float_t>(streams, fftSize, sampleRate, doPlot, filters);
		} else if(sampleLen == 2) {
			ifftcat<uint16_t, float_t>(streams, fftSize, sampleRate, doPlot, filters);
		} else if(sampleLen == 4) {
			ifftcat<uint32_t, float_t>(streams, fftSize, sampleRate, doPlot, filters);
		} else if(sampleLen == 8) {
			ifftcat<uint64_t, double_t>(streams, fftSize, sampleRate, doPlot, filters);
		} else {
			std::cerr << "sample length must be a power of 2 and less than 16" << std::endl;
			return 2;
		}
  } else {
		if(sampleLen == 1) {
			fftcat<uint8_t, float_t>(streams, fftSize, sampleRate, doPlot, filters);
		} else if(sampleLen == 2) {
			fftcat<uint16_t, float_t>(streams, fftSize, sampleRate, doPlot, filters);
		} else if(sampleLen == 4) {
			fftcat<uint32_t, float_t>(streams, fftSize, sampleRate, doPlot, filters);
		} else if(sampleLen == 8) {
			fftcat<uint64_t, double_t>(streams, fftSize, sampleRate, doPlot, filters);
		} else {
			std::cerr << "sample length must be a power of 2 and less than 16" << std::endl;
			return 2;
		}
  }

	return 0;
}
