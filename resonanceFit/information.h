///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2016 Sebastian Uhl (TUM)
//
//    This file is part of ROOTPWA
//
//    ROOTPWA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ROOTPWA is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ROOTPWA.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      storage for the general information required by the resonance fit
//      - fit results to be used
//      - waves and used mass ranges
//
//-------------------------------------------------------------------------


#ifndef RESONANCEFIT_INFORMATION_HH
#define RESONANCEFIT_INFORMATION_HH

#include <iostream>
#include <string>
#include <vector>

namespace rpwa {

	namespace resonanceFit {

		class information {

		public:

			class bin {

			public:

				bin(const std::string& fileName,
				    const double tPrimeMean,
				    const double rescaleErrors = 1.0,
				    const std::vector<std::string>& sysFileNames = std::vector<std::string>());
				~bin() {}

				const std::string& fileName() const { return _fileName; }
				double tPrimeMean() const { return _tPrimeMean; }
				double rescaleErrors() const { return _rescaleErrors; }
				const std::vector<std::string>& sysFileNames() const { return _sysFileNames; }

				std::ostream& print(std::ostream& out = std::cout, const bool newLine = true) const;

			private:

				std::string _fileName;
				double _tPrimeMean;
				double _rescaleErrors;
				std::vector<std::string> _sysFileNames;

			};

			class wave {

			public:

				wave(const std::string& waveName,
				     const std::pair<double, double>& massLimits = std::pair<double, double>(-1.0, -1.0),
				     const std::vector<std::string>& waveNameAlternatives = std::vector<std::string>());
				~wave() {}

				const std::string& waveName() const { return _waveName; }
				const std::pair<double, double>& massLimits() const { return _massLimits; }
				const std::vector<std::string>& waveNameAlternatives() const { return _waveNameAlternatives; }

				std::ostream& print(std::ostream& out = std::cout, const bool newLine = true) const;

			private:

				std::string _waveName;
				std::pair<double, double> _massLimits;
				std::vector<std::string> _waveNameAlternatives;

			};

			information(const std::vector<rpwa::resonanceFit::information::bin>& bins,
			            const std::vector<rpwa::resonanceFit::information::wave>& waves);
			~information() {}

			size_t nrBins() const { return _bins.size(); }
			const std::vector<rpwa::resonanceFit::information::bin>& bins() const { return _bins; }
			const rpwa::resonanceFit::information::bin& getBin(const size_t idxBin) const { return _bins[idxBin]; }

			size_t nrWaves() const { return _waves.size(); }
			const std::vector<rpwa::resonanceFit::information::wave>& waves() const { return _waves; }
			const rpwa::resonanceFit::information::wave& getWave(const size_t idxWave) const { return _waves[idxWave]; }

			std::ostream& print(std::ostream& out = std::cout, const bool newLine = true) const;

		private:

			std::vector<rpwa::resonanceFit::information::bin> _bins;
			std::vector<rpwa::resonanceFit::information::wave> _waves;

		};

		std::ostream& operator<< (std::ostream& out, const rpwa::resonanceFit::information& information);
		std::ostream& operator<< (std::ostream& out, const rpwa::resonanceFit::information::bin& bin);
		std::ostream& operator<< (std::ostream& out, const rpwa::resonanceFit::information::wave& wave);

	} // end namespace resonanceFit

} // end namespace rpwa


inline
std::ostream&
rpwa::resonanceFit::operator<< (std::ostream& out, const rpwa::resonanceFit::information& information)
{
	return information.print(out, false);
}


inline
std::ostream&
rpwa::resonanceFit::operator<< (std::ostream& out, const rpwa::resonanceFit::information::bin& bin)
{
	return bin.print(out, false);
}


inline
std::ostream&
rpwa::resonanceFit::operator<< (std::ostream& out, const rpwa::resonanceFit::information::wave& wave)
{
	return wave.print(out, false);
}


#endif // RESONANCEFIT_INFORMATION_HH
