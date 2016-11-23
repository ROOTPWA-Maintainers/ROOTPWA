///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010-2012 Sebastian Neubert (TUM)
//    Copyright 2015-2016 Sebastian Uhl (TUM)
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
//      final-state mass-dependence of resonance fit
//      - use ROOT's TFormula for user-definable functions
//
//-------------------------------------------------------------------------


#ifndef RESONANCEFIT_FSMD_HH
#define RESONANCEFIT_FSMD_HH

#include <complex>
#include <iostream>
#include <limits>
#include <vector>

#include <boost/multi_array.hpp>

#include "parameter.h"

namespace YAML {
	class Emitter;
	class Node;
}
class TFormula;

namespace rpwa {

	namespace resonanceFit {

		class cache;
		class parameters;

		class fsmd {

		public:

			fsmd(const size_t id,
			     const std::vector<size_t>& nrMassBins,
			     const boost::multi_array<double, 2>& massBinCenters,
			     const std::shared_ptr<TFormula>& function,
			     const boost::multi_array<rpwa::resonanceFit::parameter, 1>& parameters);
			fsmd(const size_t id,
			     const std::vector<size_t>& nrMassBins,
			     const boost::multi_array<double, 2>& massBinCenters,
			     const std::vector<std::shared_ptr<TFormula> >& functions,
			     const boost::multi_array<rpwa::resonanceFit::parameter, 2>& parameters);
			~fsmd();

			size_t getId() const { return _id; }

			bool isSameFunctionForAllBins() const { return _sameFunctionForAllBins; }

			const std::shared_ptr<TFormula>& getFunction(const size_t idxBin) const { return _functions[idxBin]; }

			size_t getNrBins() const { return _nrBins; }

			size_t getNrParameters(const size_t idxBin) const { return _nrParameters[idxBin]; }
			size_t getParameterIndex(const size_t idxBin) const { return _parametersIndex[idxBin]; }
			size_t importParameters(const double* par,
			                        rpwa::resonanceFit::parameters& fitParameters,
			                        rpwa::resonanceFit::cache& cache) const;

			const rpwa::resonanceFit::parameter& getParameter(const size_t idxBin, const size_t idxParameter) const { return _parameters[idxBin][idxParameter]; }

			std::complex<double> val(const rpwa::resonanceFit::parameters& fitParameters,
			                         rpwa::resonanceFit::cache& cache,
			                         const size_t idxBin,
			                         const double mass,
			                         const size_t idxMass = std::numeric_limits<size_t>::max()) const;

			std::ostream& print(std::ostream& out = std::cout, const bool newLine = true) const;

		private:

			const size_t _id;

			bool _sameFunctionForAllBins;

			size_t _nrBins;
			size_t _maxParameters;

			std::vector<std::vector<size_t> > _binsEqualValues;

			std::vector<std::shared_ptr<TFormula> > _functions;

			std::vector<size_t> _nrParameters;
			std::vector<size_t> _parametersIndex;

			boost::multi_array<rpwa::resonanceFit::parameter, 2> _parameters;

		};

		std::ostream& operator<< (std::ostream& out, const rpwa::resonanceFit::fsmd& fsmd);

	} // end namespace resonanceFit

} // end namespace rpwa


inline
std::ostream&
rpwa::resonanceFit::operator<< (std::ostream& out, const rpwa::resonanceFit::fsmd& fsmd)
{
	return fsmd.print(out, false);
}


#endif // RESONANCEFIT_FSMD_HH
