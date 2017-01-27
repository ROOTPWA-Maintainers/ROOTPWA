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
//      wrapper around the ROOT minimizers
//
//-------------------------------------------------------------------------


#ifndef RESONANCEFIT_MINIMIZERROOT_HH
#define RESONANCEFIT_MINIMIZERROOT_HH

#include <memory>
#include <string>
#include <vector>

#include <Math/IFunction.h>

#include "forward.h"
#include "minimizer.h"

namespace ROOT {
	namespace Math {
		class Minimizer;
	}
}

namespace rpwa {

	namespace resonanceFit {

		class minimizerRoot : public rpwa::resonanceFit::minimizer {

		private:

			class functionAdaptor : public ROOT::Math::IBaseFunctionMultiDim {

			public:

				functionAdaptor(const rpwa::resonanceFit::functionConstPtr& fitFunction);
				virtual ~functionAdaptor() {}

				virtual rpwa::resonanceFit::minimizerRoot::functionAdaptor* Clone() const;

				virtual unsigned int NDim() const;

				virtual unsigned int NPoint() const;

			private:

				virtual double DoEval(const double* par) const;

				const rpwa::resonanceFit::functionConstPtr _fitFunction;

			};

		public:

			minimizerRoot(const rpwa::resonanceFit::modelConstPtr& fitModel,
			              const rpwa::resonanceFit::functionConstPtr& fitFunction,
			              const unsigned int maxNmbOfFunctionCalls,
			              const std::string minimizerType[],
			              const int minimizerStrategy,
			              const double minimizerTolerance,
			              const bool quiet);
			virtual ~minimizerRoot();

			std::map<std::string, double> minimize(std::vector<std::string>& freeParameters,
			                                       rpwa::resonanceFit::parameters& fitParameters,
			                                       rpwa::resonanceFit::parameters& fitParametersError,
			                                       TMatrixT<double>& covarianceMatrix,
			                                       rpwa::resonanceFit::cache& cache);

		private:

			bool initParameters(const rpwa::resonanceFit::parameters& fitParameters,
			                    const std::string& freeParameters) const;

			std::unique_ptr<ROOT::Math::Minimizer> _minimizer;

			const rpwa::resonanceFit::modelConstPtr _fitModel;

			rpwa::resonanceFit::minimizerRoot::functionAdaptor _functionAdaptor;

			unsigned int _maxNmbOfIterations;
			unsigned int _maxNmbOfFunctionCalls;

		};

	} // end namespace resonanceFit

} // end namespace rpwa


#endif // RESONANCEFIT_MINIMIZERROOT_HH
