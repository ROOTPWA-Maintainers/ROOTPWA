///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010-2012 Sebastian Neubert (TUM)
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
//      final-state mass-dependence of resonance fit
//      - use ROOT's TFormula for user-definable functions
//
//-------------------------------------------------------------------------


#ifndef MASSDEPFITFSMD_HH
#define MASSDEPFITFSMD_HH

#include <iostream>
#include <limits>
#include <vector>

namespace libconfig {
	class Setting;
}
namespace ROOT {
	namespace Math {
		class Minimizer;
	}
}
class TFormula;

namespace rpwa {

	namespace massDepFit {

		class parameters;

		class fsmd {

		public:

			fsmd(const size_t id);
			~fsmd();

			bool init(const libconfig::Setting* configComponent,
			          rpwa::massDepFit::parameters& fitParameters,
			          const std::vector<double>& massBinCenters,
			          const bool debug);

			bool update(const libconfig::Setting* configComponent,
			            const rpwa::massDepFit::parameters& fitParameters,
			            const ROOT::Math::Minimizer* minimizer,
			            const bool debug) const;

			size_t getNrParameters() const { return _nrParameters; }
			size_t importParameters(const double* par,
			                        rpwa::massDepFit::parameters& fitParameters);

			bool getParameterFixed(const size_t idxParameter) const { return _parametersFixed[idxParameter]; }
			double getParameterLimitLower(const size_t idxParameter) const { return _parametersLimitLower[idxParameter]; }
			bool getParameterLimitedLower(const size_t idxParameter) const { return _parametersLimitedLower[idxParameter]; }
			double getParameterLimitUpper(const size_t idxParameter) const { return _parametersLimitUpper[idxParameter]; }
			bool getParameterLimitedUpper(const size_t idxParameter) const { return _parametersLimitedUpper[idxParameter]; }
			double getParameterStep(const size_t idxParameter) const { return _parametersStep[idxParameter]; }

			double val(const rpwa::massDepFit::parameters& fitParameters,
			           const double mass,
			           const size_t idxMass = std::numeric_limits<size_t>::max()) const;

			std::ostream& print(std::ostream& out = std::cout) const;

		private:

			const size_t _id;

			TFormula* _function;

			bool _functionFixed;

			size_t _nrParameters;

			std::vector<bool> _parametersFixed;
			std::vector<double> _parametersLimitLower;
			std::vector<bool> _parametersLimitedLower;
			std::vector<double> _parametersLimitUpper;
			std::vector<bool> _parametersLimitedUpper;
			std::vector<std::string> _parametersName;
			std::vector<double> _parametersStep;

			std::vector<double> _values;
		};

	} // end namespace massDepFit

} // end namespace rpwa


inline
std::ostream&
operator<< (std::ostream& out, const rpwa::massDepFit::fsmd& fsmd)
{
	return fsmd.print(out);
}


#endif // MASSDEPFITFSMD_HH
