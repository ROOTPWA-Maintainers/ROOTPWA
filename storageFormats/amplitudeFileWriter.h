#ifndef AMPLITUDEFILEWRITER_H
#define AMPLITUDEFILEWRITER_H

#include <complex>

#include "amplitudeMetadata.h"
#include "hashCalculator.h"


class TFile;

namespace rpwa {

	class amplitudeTreeLeaf;

	class amplitudeFileWriter {

	  public:

		amplitudeFileWriter();
		~amplitudeFileWriter();

		bool initialize(TFile&                                  outputFile,
		                const std::vector<rpwa::eventMetadata*> eventMetadata,
		                const std::string&                      keyfileContent,
		                const std::string&                      objectBasename,
		                const int&                              splitlevel = 99,
		                const int&                              buffsize = 256000);

		void addAmplitude(const std::complex<double>& amplitude);
		void addAmplitudes(const std::vector<std::complex<double> >& amplitudes);

		void reset();

		bool finalize();

		const bool& initialized() { return _initialized; }

	  private:

		bool _initialized;
		TFile* _outputFile;
		rpwa::amplitudeMetadata _metadata;
		rpwa::amplitudeTreeLeaf* _ampTreeLeaf;
		hashCalculator _hashCalculator;

	};

}


#endif
