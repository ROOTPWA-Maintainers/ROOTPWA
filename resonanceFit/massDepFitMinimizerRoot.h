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


#ifndef MASSDEPFITMINIMIZERROOT_HH
#define MASSDEPFITMINIMIZERROOT_HH

#include <memory>
#include <string>
#include <vector>

#include <Math/IFunction.h>

#include "massDepFitForward.h"
#include "massDepFitMinimizer.h"

namespace ROOT {
	namespace Math {
		class Minimizer;
	}
}

namespace rpwa {

	namespace massDepFit {

		class minimizerRoot : public rpwa::massDepFit::minimizer {

		private:

			class functionAdaptor : public ROOT::Math::IBaseFunctionMultiDim {

			public:

				functionAdaptor(const rpwa::massDepFit::functionConstPtr& fitFunction);
				virtual ~functionAdaptor() {}

				virtual rpwa::massDepFit::minimizerRoot::functionAdaptor* Clone() const;

				virtual unsigned int NDim() const;

				virtual double DoEval(const double* par) const;

			private:

				const rpwa::massDepFit::functionConstPtr _fitFunction;

			};

		public:

			minimizerRoot(const rpwa::massDepFit::modelConstPtr& fitModel,
			              const rpwa::massDepFit::functionConstPtr& fitFunction,
			              const std::vector<std::string>& freeParameters,
			              const unsigned int maxNmbOfFunctionCalls,
			              const std::string minimizerType[],
			              const int minimizerStrategy,
			              const double minimizerTolerance,
			              const bool quiet);
			virtual ~minimizerRoot();

			unsigned int getNrFreeParameters() const;

			int minimize(rpwa::massDepFit::parameters& fitParameters,
			             rpwa::massDepFit::parameters& fitParametersError,
			             rpwa::massDepFit::cache& cache);

		private:

			bool initParameters(const rpwa::massDepFit::parameters& fitParameters,
			                    const std::string& freeParameters) const;

			std::unique_ptr<ROOT::Math::Minimizer> _minimizer;

			const rpwa::massDepFit::modelConstPtr _fitModel;

			rpwa::massDepFit::minimizerRoot::functionAdaptor _functionAdaptor;

			const std::vector<std::string> _freeParameters;

			unsigned int _maxNmbOfIterations;
			unsigned int _maxNmbOfFunctionCalls;
			bool _runHesse;

		};

	} // end namespace massDepFit

} // end namespace rpwa


#endif
