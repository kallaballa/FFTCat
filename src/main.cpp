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
  size_t bufferSize = 1024;
  uint32_t sampleRate = 44100;
  uint32_t sampleLen = 1;

  po::options_description genericDesc("Options");
  genericDesc.add_options()("help,h", "Produce help message")
		("plot,p", "Write a signal plot to the terminal instead of the raw fft stream")
		("buffersize,b", po::value<size_t>(&bufferSize)->default_value(bufferSize),"The i/o buffer size in bytes")
		("samplerate,s", po::value<uint32_t>(&sampleRate)->default_value(sampleRate),"The i/o sample rate in hertz")
  	("samplelen,l", po::value<uint32_t>(&sampleLen)->default_value(sampleLen),"The i/o sample length in bytes");

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
    std::cerr << "Usage: fftcat [options] FILES..." << std::endl;
    std::cerr << visible;
    return 1;
  }


  std::vector<std::istream*> streams;

  if(files.empty() || (files.size() == 1 && files[0] == "-")) {
  	streams.push_back(&std::cin);
  } else {
  	for(const string& f : files)
  		streams.push_back(new std::ifstream(f));
  }

  if(vm.count("plot")) {
    if(sampleLen == 1) {
      plot_fft<char>(streams, bufferSize, sampleRate);
    } else if(sampleLen == 2) {
    	plot_fft<uint16_t>(streams, bufferSize, sampleRate);
    } else if(sampleLen == 4) {
    	plot_fft<uint32_t>(streams, bufferSize, sampleRate);
    } else if(sampleLen == 8) {
    	plot_fft<uint64_t>(streams, bufferSize, sampleRate);
    } else {
    	std::cerr << "sample length must be a power of 2" << std::endl;
    	return 2;
    }
  } else {
		if(sampleLen == 1) {
			cat_fft<char>(streams, bufferSize, sampleRate);
		} else if(sampleLen == 2) {
			cat_fft<uint16_t>(streams, bufferSize, sampleRate);
		} else if(sampleLen == 4) {
			cat_fft<uint32_t>(streams, bufferSize, sampleRate);
		} else if(sampleLen == 8) {
			cat_fft<uint64_t>(streams, bufferSize, sampleRate);
		} else {
			std::cerr << "sample length must be a power of 2" << std::endl;
			return 2;
		}
  }

  for(std::istream* is : streams) {
  	delete is;
  }
  return 0;
}
