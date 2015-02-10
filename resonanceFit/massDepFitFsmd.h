//-----------------------------------------------------------
//
// Description:
//    Final-state mass-dependence of mass dependent fit
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

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
