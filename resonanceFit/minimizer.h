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
//      abstract base class for wrappers around minimizers
//
//-------------------------------------------------------------------------


#ifndef RESONANCEFIT_MINIMIZER_HH
#define RESONANCEFIT_MINIMIZER_HH

#include <map>

namespace rpwa {

	namespace resonanceFit {

		class cache;
		class parameters;

		class minimizer {

		public:

			minimizer() {}
			virtual ~minimizer() {}

			virtual std::map<std::string, double> minimize(rpwa::resonanceFit::parameters& fitParameters,
			                                               rpwa::resonanceFit::parameters& fitParametersError,
			                                               rpwa::resonanceFit::cache& cache) = 0;

		};

	} // end namespace resonanceFit

} // end namespace rpwa


#endif // RESONANCEFIT_MINIMIZER_HH
