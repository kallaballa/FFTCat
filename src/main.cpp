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

int main(int argc, char** argv) {
  std::vector<string> files;
  size_t bufferSize = 256;
  uint32_t sampleRate = 44100;
  uint32_t sampleLen = 4;
  string frequencies;

  po::options_description genericDesc("Options");
  genericDesc.add_options()("help,h", "Produce help message")
		("inverse,i", "Calculate the ifft instead of the fft")
		("plot,p", "Write a signal plot to the terminal instead of the raw fft stream")
		("buffersize,b", po::value<size_t>(&bufferSize)->default_value(bufferSize),"The i/o buffer size in bytes")
		("samplerate,r", po::value<uint32_t>(&sampleRate)->default_value(sampleRate),"The i/o sample rate in hertz")
  	("samplelen,s", po::value<uint32_t>(&sampleLen)->default_value(sampleLen),"The i/o sample length in bytes")
		("filter,f", po::value<string>(&frequencies)->default_value(frequencies),"The low and high pass frequency of the filter in hertz delimited by a colon character");

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

  bool doPlot = vm.count("plot");

	std::vector<FFTFilter> filters;
	FrequencyType lowPassFreq = std::numeric_limits<FrequencyType>::max();
	FrequencyType highPassFreq = 0;
  if(!frequencies.empty()) {
  	auto pos = frequencies.find(',');
  	if(pos != string::npos) {
  		lowPassFreq = stoi(frequencies.substr(0, pos));
  		highPassFreq = stoi(frequencies.substr(pos + 1));
  	}

		if(lowPassFreq < std::numeric_limits<FrequencyType>::max() && highPassFreq > 0) {
			filters.push_back({lowPassFreq, highPassFreq});
		}
  }

  //TODO type-selection based on the platform.
  if(vm.count("inverse")) {
		if(sampleLen == 1) {
			ifftcat<uint8_t, float_t>(streams, bufferSize, sampleRate, doPlot, filters);
		} else if(sampleLen == 2) {
			ifftcat<uint16_t, float_t>(streams, bufferSize, sampleRate, doPlot, filters);
		} else if(sampleLen == 4) {
			ifftcat<uint32_t, float_t>(streams, bufferSize, sampleRate, doPlot, filters);
		} else if(sampleLen == 8) {
			ifftcat<uint64_t, double_t>(streams, bufferSize, sampleRate, doPlot, filters);
		} else {
			std::cerr << "sample length must be a power of 2 and less than 16" << std::endl;
			return 2;
		}
  } else {
		if(sampleLen == 1) {
			fftcat<uint8_t, float_t>(streams, bufferSize, sampleRate, doPlot, filters);
		} else if(sampleLen == 2) {
			fftcat<uint16_t, float_t>(streams, bufferSize, sampleRate, doPlot, filters);
		} else if(sampleLen == 4) {
			fftcat<uint32_t, float_t>(streams, bufferSize, sampleRate, doPlot, filters);
		} else if(sampleLen == 8) {
			fftcat<uint64_t, double_t>(streams, bufferSize, sampleRate, doPlot, filters);
		} else {
			std::cerr << "sample length must be a power of 2 and less than 16" << std::endl;
			return 2;
		}
  }

	return 0;
}
