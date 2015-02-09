//-----------------------------------------------------------
//
// Description:
//      mass dependent fit likelihood rank 1!
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef MASSDEPFITPARAMETERS_HH
#define MASSDEPFITPARAMETERS_HH

#include <complex>

#include <boost/multi_array.hpp>

namespace rpwa {

	namespace massDepFit {

		class parameters {

		public:

			parameters();
			parameters(const size_t maxComponents,
			           const size_t maxChannels,
			           const size_t maxParameters,
			           const size_t maxBins);
			virtual ~parameters() {}

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

			std::ostream& print(std::ostream& out = std::cout) const;

		private:

			const bool _fixed;

			boost::multi_array<std::complex<double>, 2> _branchings;
			boost::multi_array<std::complex<double>, 3> _couplings;
			boost::multi_array<double, 2> _parameters;

		};

	} // end namespace massDepFit

} // end namespace rpwa


inline
std::ostream&
operator<< (std::ostream& out, const rpwa::massDepFit::parameters& parameters)
{
	return parameters.print(out);
}


#endif
