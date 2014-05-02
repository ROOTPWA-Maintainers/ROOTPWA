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

#ifndef MASSDEPFITLIKELI_HH
#define MASSDEPFITLIKELI_HH

#include <boost/multi_array.hpp>

#include <Math/IFunction.h>
#include <TMatrixT.h>

namespace rpwa {

	namespace massDepFit {

		class model;

		class likelihood : public ROOT::Math::IBaseFunctionMultiDim {

		public:

			likelihood() {}
			virtual ~likelihood() {}

			virtual likelihood* Clone() const;

			virtual unsigned int NDim() const;

			virtual double DoEval(const double* par) const;

			unsigned int NDataPoints() const;

			bool init(rpwa::massDepFit::model* compset,
			          const std::vector<double>& massBinCenters,
			          const boost::multi_array<std::complex<double>, 3>& productionAmplitudes,
			          const boost::multi_array<double, 6>& productionAmplitudesCovariance,
			          const boost::multi_array<std::complex<double>, 4>& spinDensityMatrices,
			          const boost::multi_array<double, 6>& spinDensityCovarianceMatrices,
			          const boost::multi_array<std::pair<size_t, size_t>, 2>& wavePairMassBinLimits);

			bool fitProductionAmplitudes() const { return _fitProductionAmplitudes; }
			void fitProductionAmplitudes(const bool val) { _fitProductionAmplitudes = val; }

			bool useCovariance() const { return _useCovariance; }
			void useCovariance(const bool val) { _useCovariance = val; }

		private:

			double DoEvalProductionAmplitudes() const;
			double DoEvalSpinDensityMatrix() const;

			rpwa::massDepFit::model* _compset;

			size_t _nrBins;
			size_t _nrMassBins;
			size_t _nrWaves;

			size_t _idxMassMin;
			size_t _idxMassMax;

			std::vector<double> _massBinCenters;

			boost::multi_array<std::complex<double>, 3> _productionAmplitudes;
			boost::multi_array<double, 6> _productionAmplitudesCovariance;
			boost::multi_array<TMatrixT<double>, 2> _productionAmplitudesCovMatInv;

			boost::multi_array<std::complex<double>, 4> _spinDensityMatrices;
			boost::multi_array<double, 6> _spinDensityCovarianceMatrices;

			boost::multi_array<std::pair<size_t, size_t>, 2> _wavePairMassBinLimits;

			size_t _idxAnchorWave;

			bool _fitProductionAmplitudes;
			bool _useCovariance;

		};

	} // end namespace massDepFit

} // end namespace rpwa

#endif
