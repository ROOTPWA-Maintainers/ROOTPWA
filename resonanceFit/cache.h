///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2015 Sebastian Uhl (TUM)
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
//      cache results during the resonance fit
//
//-------------------------------------------------------------------------


#ifndef RESONANCEFIT_CACHE_HH
#define RESONANCEFIT_CACHE_HH

#include <complex>
#include <iostream>

#include <boost/multi_array.hpp>

namespace rpwa {

	namespace resonanceFit {

		class cache {

		public:

			cache(const size_t maxWaves,
			      const size_t maxComponents,
			      const size_t maxChannels,
			      const size_t maxBins,
			      const size_t maxMassBins);
			~cache() {}

			std::complex<double> getCoupling(const size_t idxComponent, const size_t idxChannel, const size_t idxBin, const size_t idxMassBin) const { return _couplings[idxComponent][idxChannel][idxBin][idxMassBin]; }
			std::complex<double> getComponent(const size_t idxComponent, const size_t idxBin, const size_t idxMassBin) const { return _components[idxComponent][idxBin][idxMassBin]; }
			std::complex<double> getProdAmp(const size_t idxWave, const size_t idxBin, const size_t idxMassBin) const { return _prodAmps[idxWave][idxBin][idxMassBin]; }

			void setCoupling(const size_t idxComponent, const size_t idxChannel, const size_t idxBin, const size_t idxMassBin, const std::complex<double> coupling);
			void setComponent(const size_t idxComponent, const size_t idxBin, const size_t idxMassBin, const std::complex<double> component);
			void setProdAmp(const size_t idxWave, const size_t idxBin, const size_t idxMassBin, const std::complex<double> prodAmp);

			std::ostream& print(std::ostream& out = std::cout) const;

		private:

			boost::multi_array<std::complex<double>, 4> _couplings;
			boost::multi_array<std::complex<double>, 3> _components;
			boost::multi_array<std::complex<double>, 3> _prodAmps;

		};

		std::ostream& operator<< (std::ostream& out, const rpwa::resonanceFit::cache& cache);

	} // end namespace resonanceFit

} // end namespace rpwa


inline
std::ostream&
rpwa::resonanceFit::operator<< (std::ostream& out, const rpwa::resonanceFit::cache& cache)
{
	return cache.print(out);
}


#endif // RESONANCEFIT_CACHE_HH
