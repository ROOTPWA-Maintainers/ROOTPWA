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
//      storage for the data required by the resonance fit
//      - vectors of production amplitudes
//      - spin-density matrices
//      - the corresponding covariances
//      - phase-space integrals of waves
//
//-------------------------------------------------------------------------


#ifndef RESONANCEFIT_DATA_HH
#define RESONANCEFIT_DATA_HH

#include <boost/multi_array.hpp>

namespace rpwa {

	namespace resonanceFit {

		class data {

		public:

			data(const std::vector<size_t>& nrMassBins,
			     const boost::multi_array<double, 2>& massBinCenters);
			~data() {}

			size_t nrBins() const { return _nrMassBins.size(); }
			size_t maxMassBins() const { return *(std::max_element(_nrMassBins.begin(), _nrMassBins.end())); }

			const std::vector<size_t>& nrMassBins() const { return _nrMassBins; }
			const boost::multi_array<double, 2>& massBinCenters() const { return _massBinCenters; }

		private:

			std::vector<size_t> _nrMassBins;
			boost::multi_array<double, 2> _massBinCenters;

		};

	} // end namespace resonanceFit

} // end namespace rpwa

#endif // RESONANCEFIT_DATA_HH
