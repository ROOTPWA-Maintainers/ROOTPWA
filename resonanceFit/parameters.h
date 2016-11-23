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
//      parameters for resonance fit
//
//-------------------------------------------------------------------------


#ifndef RESONANCEFIT_PARAMETERS_HH
#define RESONANCEFIT_PARAMETERS_HH

#include <complex>
#include <iostream>

#include <boost/multi_array.hpp>

namespace rpwa {

	namespace resonanceFit {

		class parameters {

		public:

			parameters();
			parameters(const size_t maxComponents,
			           const size_t maxChannels,
			           const size_t maxParameters,
			           const size_t maxBins);
			~parameters() {}

			std::complex<double> getBranching(const size_t idxComponent, const size_t idxChannel) const { return _branchings[idxComponent][idxChannel]; };
			std::complex<double> getCoupling(const size_t idxComponent, const size_t idxChannel, const size_t idxBin) const { return _couplings[idxComponent][idxChannel][idxBin]; };
			double getParameter(const size_t idxComponent, const size_t idxParameter) const { return _parameters[idxComponent][idxParameter]; };
			const double* getParameters(const size_t idxComponent) const { return _parameters[idxComponent].origin(); };

			void setBranching(const size_t idxComponent, const size_t idxChannel, const std::complex<double> branching) { _branchings[idxComponent][idxChannel] = branching; }
			void setCoupling(const size_t idxComponent, const size_t idxChannel, const size_t idxBin, const std::complex<double> coupling) { _couplings[idxComponent][idxChannel][idxBin] = coupling; }
			void setParameter(const size_t idxComponent, const size_t idxParameter, const double parameter) { _parameters[idxComponent][idxParameter] = parameter; }

			void resize(const size_t maxComponents,
			            const size_t maxChannels,
			            const size_t maxParameters,
			            const size_t maxBins);

			std::ostream& print(std::ostream& out = std::cout, const bool newLine = true) const;

		private:

			const bool _fixedSize;

			boost::multi_array<std::complex<double>, 2> _branchings;
			boost::multi_array<std::complex<double>, 3> _couplings;
			boost::multi_array<double, 2> _parameters;

		};

		std::ostream& operator<< (std::ostream& out, const rpwa::resonanceFit::parameters& parameters);

	} // end namespace resonanceFit

} // end namespace rpwa


inline
std::ostream&
rpwa::resonanceFit::operator<< (std::ostream& out, const rpwa::resonanceFit::parameters& parameters)
{
	return parameters.print(out, false);
}


#endif // RESONANCEFIT_PARAMETERS_HH
