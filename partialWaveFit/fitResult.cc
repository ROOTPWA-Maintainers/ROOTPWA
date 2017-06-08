///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------
//
// Description:
//      Data storage class for PWA fit result of one kinematic bin
//
// Environment:
//      Software developed for the COMPASS experiment at CERN
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#include <algorithm>
#include <set>

#include <boost/multi_array.hpp>

#include "Math/SpecFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"
#include "TDecompChol.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"

#include "fitResult.h"


using namespace std;
using namespace rpwa;


ClassImp(fitResult);


namespace {


	/// \brief returns matrix representation of complex number
	///
	///   c.re  -c.im
	///   c.im   c.re
	// !!! possible optimization: return TMatrixTLazy
	inline
	TMatrixT<double>
	matrixRepr(const std::complex<double>& c)
	{
		TMatrixT<double> m(2, 2);
		m[0][0] = m[1][1] = c.real();
		m[0][1] = -c.imag();
		m[1][0] =  c.imag();
		return m;
	}


}


fitResult::fitResult()
	: _nmbEvents     (0),
	  _normNmbEvents (0),
	  _logLikelihood (0),
	  _rank          (0),
	  _covMatrixValid(false),
	  _converged     (false),
	  _hasHessian    (false)
{ }


fitResult::fitResult(const fitResult& result, const bool fillCovMatrix, const bool fillIntegralMatrices)
	: TObject()
{
	fill(result, fillCovMatrix, fillIntegralMatrices);
}

fitResult::fitResult(const fitResult&                          result,
                     const TMatrixT<double>*                   fitParCovMatrix,         // covariance matrix of fit parameters
                     const std::vector<std::pair<int, int> >*  fitParCovMatrixIndices,  // indices of fit parameters for real and imaginary part in covariance matrix
                     const rpwa::complexMatrix*                normIntegral,            // normalization integral matrix
                     const rpwa::complexMatrix*                acceptedNormIntegral,    // normalization integral matrix with acceptance
                     const std::vector<double>*                phaseSpaceIntegral)     // normalization integral over full phase space without acceptance
	: TObject()
{
	fill(result, fitParCovMatrix, fitParCovMatrixIndices, normIntegral, acceptedNormIntegral, phaseSpaceIntegral);
}


fitResult::~fitResult()
{ }


fitResult*
fitResult::variedProdAmps() const
{
	// copy complete info
	fitResult* result = new fitResult(*this);

	// Cholesky decomposition
	// hopefully this will not slow us down too much; if it does we have to cache
	printDebug << "starting Cholesky decomposition... " << endl;
	TDecompChol decomp(_fitParCovMatrix);
	decomp.Decompose();
	const TMatrixT<Double_t> C(decomp.GetU());
	printDebug << "... done." << endl;

	const unsigned int npar = C.GetNrows();
	TMatrixD x(npar, 1);
	// generate npar independent random numbers
	for (unsigned int ipar = 0; ipar < npar; ++ipar)
		x(ipar, 0) = gRandom->Gaus();

	// this tells us how to permute the parameters taking into account all
	// correlations in the covariance matrix
	const TMatrixD y(C, TMatrixD::kMult, x);

	// now we need mapping from parameters to production amp re and im
	// loop over production amps
	const unsigned int nt = _prodAmps.size();
	for (unsigned int it = 0; it < nt; ++it) {
		const Int_t jre = _fitParCovMatrixIndices[it].first;
		const Int_t jim = _fitParCovMatrixIndices[it].second;
		double      re  = result->_prodAmps[it].Re();
		double      im  = result->_prodAmps[it].Im();
		if(jre > -1)
			re += y(jre, 0);
		if(jim > -1)
			im += y(jim, 0);
		result->_prodAmps[it] = TComplex(re, im);
		printDebug << "varying production amplitude [" << it << "] by ("
		           << ((jre > -1) ? y(jre, 0) : 0) << ", "
		           << ((jim > -1) ? y(jim, 0) : 0) << ")" << endl;
	}
	return result;
}


void
fitResult::reset()
{
	_nmbEvents     = 0;
	_normNmbEvents = 0;
	_multibinBoundaries.clear();
	_logLikelihood = 0;
	_rank          = 0;
	_prodAmps.clear();
	_prodAmpNames.clear();
	_waveNames.clear();
	_covMatrixValid = false;
	_fitParCovMatrix.ResizeTo(0, 0);
	_fitParCovMatrixIndices.clear();
	_normIntegral.resizeTo(0, 0);
	_acceptedNormIntegral.resizeTo(0, 0);
	_normIntIndexMap.clear();
	_phaseSpaceIntegral.clear();
	_converged  = false;
	_hasHessian = false;
}


void
fitResult::fill(const unsigned int              nmbEvents,               // number of events in bin
                const unsigned int              normNmbEvents,           // number of events to normalize to
                const multibinBoundariesType&   multibinBoundaries,      // multibin boundaries
                const double                    logLikelihood,           // log(likelihood) at maximum
                const int                       rank,                    // rank of fit
                const vector<TComplex>&         prodAmps,                // production amplitudes
                const vector<string>&           prodAmpNames,            // names of production amplitudes used in fit
                const TMatrixT<double>*         fitParCovMatrix,         // covariance matrix of fit parameters
                const vector<pair<int, int> >&  fitParCovMatrixIndices,  // indices of fit parameters for real and imaginary part in covariance matrix matrix
                const complexMatrix*            normIntegral,            // normalization integral matrix
                const complexMatrix*            acceptedNormIntegral,    // normalization integral matrix with acceptance
                const vector<double>*           phaseSpaceIntegral,      // normalization integral over full phase space without acceptance
                const bool                      converged,
                const bool                      hasHessian)
{
	_converged              = converged;
	_nmbEvents              = nmbEvents;
	_normNmbEvents          = normNmbEvents;
	_multibinBoundaries     = multibinBoundaries;
	_logLikelihood          = logLikelihood;
	_rank                   = rank;
	_prodAmps               = prodAmps;
	_prodAmpNames           = prodAmpNames;
	_fitParCovMatrixIndices = fitParCovMatrixIndices;

	if (fitParCovMatrix != nullptr) {
		_fitParCovMatrix.ResizeTo(fitParCovMatrix->GetNrows(), fitParCovMatrix->GetNcols());
		_fitParCovMatrix = *fitParCovMatrix;
		_hasHessian         = hasHessian;
		// check whether there really is an error matrix
		if (not (fitParCovMatrix->GetNrows() == 0) and not (fitParCovMatrix->GetNcols() == 0))
			_covMatrixValid = true;
		else
			_covMatrixValid = false;
		if (_covMatrixValid) _hasHessian = true; // there has to be a Hessian if the covariance matrix is valid
	} else {
		_fitParCovMatrix.ResizeTo(0, 0);
		_covMatrixValid = false;
		_hasHessian     = false; // set to false, does not matter what is given externally
	}

	if (normIntegral != nullptr) {
		_normIntegral.resizeTo(normIntegral->nRows(), normIntegral->nCols());
		_normIntegral = *normIntegral;
	} else {
		_normIntegral.resizeTo(0, 0);
	}

	if (acceptedNormIntegral != nullptr) {
		_acceptedNormIntegral.resizeTo(acceptedNormIntegral->nRows(), acceptedNormIntegral->nCols());
		_acceptedNormIntegral = *acceptedNormIntegral;
	} else {
		_acceptedNormIntegral.resizeTo(0, 0);
	}

	if (phaseSpaceIntegral != nullptr) {
		_phaseSpaceIntegral = *phaseSpaceIntegral;
	} else {
		_phaseSpaceIntegral.resize(0);
	}

	for (unsigned int i = 0; i < _prodAmpNames.size(); ++i) {
		const string waveName = waveNameForProdAmp(i);
		const size_t waveIdx  = find(_waveNames.begin(), _waveNames.end(), waveName) - _waveNames.begin();
		if (waveIdx == _waveNames.size())
			_waveNames.push_back(waveName);
		_normIntIndexMap[i] = waveIdx;
	}
	if (not (normIntegral != nullptr or acceptedNormIntegral != nullptr)){
		_normIntIndexMap.clear();
	}

	// check consistency
	if (_prodAmps.size() != _prodAmpNames.size())
		printWarn << "number of production amplitudes (" << _prodAmps.size() << ") "
		          << "does not match number of production amplitude names "
		          << "(" << _prodAmpNames.size() << ")." << endl;
	if (_waveNames.size() > _prodAmpNames.size())
		printWarn << "number of wave names (" << _waveNames.size() << ") "
		          << "larger than number of production amplitude names "
		          << "(" << _prodAmpNames.size() << ")." << endl;
	if (_prodAmps.size() != _fitParCovMatrixIndices.size())
		printWarn << "number of production amplitudes (" << _prodAmps.size() << ") "
				<< "does not match number of covariance matrix indices "
				<< "(" << _fitParCovMatrixIndices.size() << ")." << endl;
	if (fitParCovMatrix != nullptr){
		if (_fitParCovMatrix.GetNrows() != _fitParCovMatrix.GetNcols())
			printWarn << "covariance matrix is not a square matrix "
					<< "(" << _fitParCovMatrix.GetNrows() << ", " << _fitParCovMatrix.GetNcols() << ")." << endl;
		if (_fitParCovMatrix.GetNrows() == 0) {
			// covariance matrix does not exist
			printWarn << "covariance matrix is empty." << endl;
		} else {
			// covariance matrix exists
			for (unsigned int i=0; i<_fitParCovMatrixIndices.size(); ++i) {
				if (_fitParCovMatrixIndices[i].first < 0)
					printWarn << "entry in covariance matrix for real part of production "
							<< "amplitude " << i << " does not exist." << endl;
				if (_fitParCovMatrixIndices[i].first >= _fitParCovMatrix.GetNrows())
					printWarn << "real part of production amplitude " << i << " is mapped to "
							<< "entry in covariance matrix outside covariance matrix size "
							<< "(" << _fitParCovMatrixIndices.size() << ")." << endl;
				if (_fitParCovMatrixIndices[i].second >= _fitParCovMatrix.GetNrows())
					printWarn << "imaginary part of production amplitude " << i << " is mapped to "
							<< "entry in covariance matrix outside covariance matrix size "
							<< "(" << _fitParCovMatrixIndices.size() << ")." << endl;
			}
		}
	}
	if (normIntegral != nullptr){
		if (_normIntegral.nRows() != _normIntegral.nCols())
			printWarn << "normalization integral is not a square matrix "
					<< "(" << _normIntegral.nRows() << ", " << _normIntegral.nCols() << ")." << endl;
		if (_normIntegral.nRows() == 0) {
			// normalization integral matrix does not exist
			printWarn << "normalization integral matrix is empty." << endl;
		} else {
			if (_waveNames.size() != _normIntegral.nRows())
				printWarn << "number of waves (" << _waveNames.size() << ") "
						<< "does not match size of normalization integral "
						<< "(" << _normIntegral.nRows() << ", " << _normIntegral.nCols() << ")." << endl;
		}
	}
	if (acceptedNormIntegral != nullptr){
		if (_acceptedNormIntegral.nRows() != _acceptedNormIntegral.nCols())
			printWarn << "normalization integral is not a square matrix "
					<< "(" << _acceptedNormIntegral.nRows() << ", " << _acceptedNormIntegral.nCols() << ")." << endl;
		if (_acceptedNormIntegral.nRows() == 0) {
			// acceptance integral matrix does not exist
			printWarn << "acceptance integral matrix is empty." << endl;
		} else {
			if (_waveNames.size() != _acceptedNormIntegral.nRows())
				printWarn << "number of waves (" << _waveNames.size() << ") "
						<< "does not match size of acceptance integral "
						<< "(" << _acceptedNormIntegral.nRows() << ", " << _acceptedNormIntegral.nCols() << ")." << endl;
		}
	}
	if (phaseSpaceIntegral != nullptr){
		if (_phaseSpaceIntegral.size() == 0) {
			// phase-space integral does not exist
			printWarn << "phase-space integral is empty." << endl;
		} else {
			if (_waveNames.size() != _phaseSpaceIntegral.size())
				printWarn << "number of waves (" << _waveNames.size() << ") "
						<< "does not match size of phase-space integral "
						<< "(" << _phaseSpaceIntegral.size() << ")." << endl;
		}
	}
	if (normIntegral != nullptr or acceptedNormIntegral != nullptr){
		if (_prodAmps.size() != _normIntIndexMap.size())
			printWarn << "number of production amplitudes (" << _prodAmps.size() << ") "
					<< "does not match number of mappings from production amplitudes to indices in integrals "
					<< "(" << _normIntIndexMap.size() << ")." << endl;
		for (unsigned int i=0; i<_prodAmps.size(); ++i) {
			if (_normIntIndexMap.count(i) != 1)
				printWarn << "production amplitude at index " << i << " is not mapped to any index in integrals." << endl;
			if (_normIntIndexMap[i] < 0 or _normIntIndexMap[i] >= (int)_waveNames.size())
				printWarn << "production amplitude at index " << i << " is mapped to integrals "
							"at index " << _normIntIndexMap[i] <<", which does not exist." << endl;
		}
	}
}

void
fitResult::fill(const unsigned int              nmbEvents,               // number of events in bin
                const unsigned int              normNmbEvents,           // number of events to normalize to
                const multibinBoundariesType&   multibinBoundaries,      // multibin boundaries
                const double                    logLikelihood,           // log(likelihood) at maximum
                const int                       rank,                    // rank of fit
                const vector<complex<double> >& prodAmps,                // production amplitudes
                const vector<string>&           prodAmpNames,            // names of production amplitudes used in fit
                const TMatrixT<double>*         fitParCovMatrix,         // covariance matrix of fit parameters
                const vector<pair<int, int> >&  fitParCovMatrixIndices,  // indices of fit parameters for real and imaginary part in covariance matrix matrix
                const complexMatrix*            normIntegral,            // normalization integral matrix
                const complexMatrix*            acceptedNormIntegral,    // normalization integral matrix with acceptance
                const vector<double>*           phaseSpaceIntegral,      // normalization integral over full phase space without acceptance
                const bool                      converged,
                const bool                      hasHessian)
{
	// In principle the production amplitudes could be written directly to _prodAmps and an empty prodAmps vector could be given to fill.
	// This would remove the additional copy of the prodAmps vector.
	// However, this leads to a warning in the fill function.
	vector<TComplex> prodAmpsT;
	prodAmpsT.resize(prodAmps.size());
	for (unsigned int i = 0; i < prodAmps.size(); ++i)
		prodAmpsT[i] = TComplex(prodAmps[i].real(), prodAmps[i].imag());
	fill(nmbEvents, normNmbEvents, multibinBoundaries, logLikelihood, rank, prodAmpsT, prodAmpNames, fitParCovMatrix, fitParCovMatrixIndices, normIntegral,
	     acceptedNormIntegral, phaseSpaceIntegral, converged, hasHessian);
}

void
fitResult::fill(const unsigned int              nmbEvents,               // number of events in bin
                const unsigned int              normNmbEvents,           // number of events to normalize to
                const multibinBoundariesType&   multibinBoundaries,      // multibin boundaries
                const double                    logLikelihood,           // log(likelihood) at maximum
                const int                       rank,                    // rank of fit
                const vector<complex<double> >& prodAmps,                // production amplitudes
                const vector<string>&           prodAmpNames,            // names of production amplitudes used in fit
                const TMatrixT<double>&         fitParCovMatrix,         // covariance matrix of fit parameters
                const vector<pair<int, int> >&  fitParCovMatrixIndices,  // indices of fit parameters for real and imaginary part in covariance matrix matrix
                const complexMatrix&            normIntegral,            // normalization integral matrix
                const complexMatrix&            acceptedNormIntegral,    // normalization integral matrix with acceptance
                const vector<double>&           phaseSpaceIntegral,      // normalization integral over full phase space without acceptance
                const bool                      converged,
                const bool                      hasHessian)
{
	fill(nmbEvents, normNmbEvents, multibinBoundaries, logLikelihood, rank, prodAmps, prodAmpNames,
	     &fitParCovMatrix, fitParCovMatrixIndices, &normIntegral, &acceptedNormIntegral, &phaseSpaceIntegral,
	     converged, hasHessian);
}


void
fitResult::fill(const fitResult& result,
                const TMatrixT<double>*         fitParCovMatrix,         // covariance matrix of fit parameters
                const vector<pair<int, int> >*  fitParCovMatrixIndices,  // indices of fit parameters for real and imaginary part in covariance matrix matrix
                const complexMatrix*            normIntegral,            // normalization integral matrix
                const complexMatrix*            acceptedNormIntegral,    // normalization integral matrix with acceptance
                const vector<double>*           phaseSpaceIntegral)      // normalization integral over full phase space without acceptance
{
	if (fitParCovMatrix == nullptr)
		fitParCovMatrix = (not result.covMatrixValid()) ? nullptr : &result.fitParCovMatrix();
	if (fitParCovMatrixIndices == nullptr)
		fitParCovMatrixIndices = &result.fitParCovIndices();
	if (normIntegral == nullptr)
		normIntegral = (result.normIntegralMatrix().nRows() == 0 and result.normIntegralMatrix().nCols() == 0) ? nullptr : &result.normIntegralMatrix();
	if (acceptedNormIntegral == nullptr)
		acceptedNormIntegral = (result.acceptedNormIntegralMatrix().nRows() == 0 and result.acceptedNormIntegralMatrix().nCols() == 0) ?
		                       nullptr : &result.acceptedNormIntegralMatrix();
	if (phaseSpaceIntegral == nullptr)
		phaseSpaceIntegral = (result.phaseSpaceIntegralVector().size() == 0) ? nullptr : &result.phaseSpaceIntegralVector();
	fill(result.nmbEvents(), result.normNmbEvents(), result.multibinBoundaries(), result.logLikelihood(),
	     result.rank(), result.prodAmps(), result.prodAmpNames(),
	     fitParCovMatrix,
	     *fitParCovMatrixIndices,
	     normIntegral,
	     acceptedNormIntegral,
	     phaseSpaceIntegral,
	     result.converged(), result.hasHessian());
}


void
fitResult::fill(const fitResult& result, bool fillCovMatrix, bool fillIntegralMatrices)
{
	if (fillCovMatrix)
		fillCovMatrix = result.covMatrixValid();
	bool fillNormIntMatrix = false;
	bool fillAccNormIntMatrix = false;
	bool fillPhaseSpace = false;
	if (fillIntegralMatrices){
		fillNormIntMatrix    = result.normIntegralMatrix().nRows() != 0         or result.normIntegralMatrix().nCols() != 0;
		fillAccNormIntMatrix = result.acceptedNormIntegralMatrix().nRows() != 0 or result.acceptedNormIntegralMatrix().nCols() != 0;
		fillPhaseSpace       = result.phaseSpaceIntegralVector().size() != 0;
	}
	fill(result.nmbEvents(), result.normNmbEvents(), result.multibinBoundaries(), result.logLikelihood(),
	     result.rank(), result.prodAmps(), result.prodAmpNames(),
	     (fillCovMatrix)? &result.fitParCovMatrix() : nullptr,
	     result.fitParCovIndices(),
	     (fillNormIntMatrix)? &result.normIntegralMatrix() : nullptr,
	     (fillAccNormIntMatrix)? &result.acceptedNormIntegralMatrix() : nullptr,
	     (fillPhaseSpace)? &result.phaseSpaceIntegralVector() : nullptr,
	     result.converged(), result.hasHessian());
}


multibinCenterType
fitResult::multibinCenter() const
{
	multibinCenterType multibinCenter;
	for (multibinBoundariesType::const_iterator it = _multibinBoundaries.begin(); it != _multibinBoundaries.end(); ++it) {
		multibinCenter[it->first] = 0.5 * (it->second.first + it->second.second);
	}
	return multibinCenter;
}


/// \brief calculates the model evidence using Rafters occam factor method
///
/// \f[ \ln P(\mathrm{Data}|M_k)\approx \ln P(\mathrm{Data}|A_{\mathrm{ML}}^k,M_k) + \ln{\frac{\sqrt{(2\pi)^d|\mathbf{C}_{A|D}|}}{\sum_i\frac{\pi}{\Psi_{ii}}}} \f]
/// This uses a flat prior in each parameter (see note on Model Selection) such
/// that no single wave can have more than the total intensity measured
double
fitResult::evidence() const
{
	double retval = 0.;
	const std::vector<double> summands = evidenceComponents();
	for(unsigned int i = 0; i < summands.size(); ++i) {
		retval += summands[i];
	}
	return retval;
}


std::vector<double>
fitResult::evidenceComponents() const
{
	std::vector<double> components(4, -numeric_limits<double>::infinity());

	// find the thresholded production amplitudes, assume those are the
	// ones with imaginary and real part equal to zero
	set<unsigned int> thrProdAmpIndices;
	for (unsigned int i = 0; i < nmbProdAmps(); ++i) {
		if (prodAmp(i) == 0.) {
			thrProdAmpIndices.insert(i);
		}
	}

	// from the list of thresholded production amplitudes get the list of
	// thresholded waves, this assumes that if a production amplitude is
	// thresholded, it is thresholded in the same way for all ranks
	set<unsigned int> thrWaveIndices;
	for (set<unsigned int>::const_iterator it = thrProdAmpIndices.begin();
	     it != thrProdAmpIndices.end(); ++it) {
		const int index = waveIndex(waveNameForProdAmp(*it));
		assert(index > 0);
		thrWaveIndices.insert(index);
	}

	// REMOVE CONSTRAINT TO NUMBER OF EVENTS!
	const double l = -logLikelihood();// - intensity(".*");
	components[0]  = l;

	if (covMatrixValid()) {
		// for the list of thresholded production amplitudes get the list of
		// columns and rows of the covariance matrix that are thresholded
		set<Int_t> thrColAndRowIndices;
		for (set<unsigned int>::const_iterator it = thrProdAmpIndices.begin();
		     it != thrProdAmpIndices.end(); ++it) {
			thrColAndRowIndices.insert(fitParCovIndices()[*it].first);
			if (fitParCovIndices()[*it].second != -1)
				thrColAndRowIndices.insert(fitParCovIndices()[*it].second);
		}

		// create a new covariance matrix with the thresholded entries removed
		TMatrixT<Double_t> thrCovMatrix(fitParCovMatrix().GetNrows()-thrColAndRowIndices.size(),
		                                fitParCovMatrix().GetNcols()-thrColAndRowIndices.size());
		Int_t row = -1;
		for (Int_t i = 0; i < fitParCovMatrix().GetNrows(); ++i) {
			if (thrColAndRowIndices.count(i) > 0) {
				continue;
			}
			row++;

			Int_t col = -1;
			for (Int_t j = 0; j < fitParCovMatrix().GetNcols(); ++j) {
				if (thrColAndRowIndices.count(j) > 0) {
					continue;
				}
				col++;

				thrCovMatrix(row, col) = fitParCovMatrix()(i, j);
			}
		}

		double       logDet = thrCovMatrix.Determinant();
		const double d      = (double)thrCovMatrix.GetNcols();

		if(std::isinf(logDet)) {
			printWarn << "found infinite determinant of covariance matrix, trying to calculate log(det) directly..." << endl;
			logDet = 0.;
			TVectorT<double> eigenvalues;
			thrCovMatrix.EigenVectors(eigenvalues);
			for(int i = 0; i < eigenvalues.GetNrows(); ++i) {
				logDet += TMath::Log(eigenvalues[i]);
			}
			if(std::isinf(logDet) or std::isnan(logDet)) {
				printWarn << "calculation failed, log(det) = " << logDet << "." << endl;
			} else {
				printSucc << "calculation succeeded, log(det) = " << logDet << "." << endl;
			}
		} else {
			logDet = TMath::Log(logDet);
		}

		// parameter-volume after observing data
		const double lvad = 0.5 * (d * 1.837877066 + logDet);
		components[1] = lvad;

		// parameter volume prior to observing the data
		if (acceptedNormIntegralMatrix().nCols() != 0 && acceptedNormIntegralMatrix().nRows() != 0) {
			// keep list of required waves for each reflectivity and rank
			boost::multi_array<vector<int>, 2> waves(boost::extents[2][rank()]);

			// loop over production amplitudes and extract combinations of
			// rank and reflectivity
			for (unsigned int i = 0; i < nmbProdAmps(); ++i) {
				// skip thresholded production amplitudes
				if (thrProdAmpIndices.count(i) > 0)
					continue;

				const string waveName = string(waveNameForProdAmp(i));
				const int waveI = waveIndex(waveName);
				assert(waveI >= 0);
				assert(thrWaveIndices.count(waveI) == 0);

				// skip flat wave (handled below)
				if (waveName == "flat")
					continue;

				const int rank = rankOfProdAmp(i);
				assert(rank >= 0 && (unsigned int)rank < this->rank());

				const int refl = partialWaveFitHelper::getReflectivity(waveName);
				assert(refl == -1 || refl == +1);

				if (refl > 0)
					waves[1][rank].push_back(waveI);
				else
					waves[0][rank].push_back(waveI);
			}

			double logDetAcc = 0;
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {
				for (unsigned int iRank = 0; iRank < rank(); ++iRank) {
					const unsigned int dim = waves[iRefl][iRank].size();
					complexMatrix sub(dim, dim);

					for(unsigned int i = 0; i < dim; ++i)
						for(unsigned int j = 0; j < dim; ++j)
							sub.set(i, j, acceptedNormIntegral(waves[iRefl][iRank][i], waves[iRefl][iRank][j]));

					logDetAcc += TMath::Log(sub.determinant().real());
				}
			}

			// parameter volume for flat wave correlates with overall acceptance
			const int idxFlat = waveIndex("flat");
			assert(idxFlat != -1);
			logDetAcc += TMath::Log(acceptedNormIntegral(idxFlat, idxFlat).real());

			// n-Sphere:
			const double lva = TMath::Log(d) + 0.5 * (d * 1.144729886 + (d - 1) * TMath::Log(nmbEvents()))
			                   - ROOT::Math::lgamma(0.5 * d + 1) - 0.5 * logDetAcc;
			components[2] = -lva;
		} else {
			printWarn << "fitResult does not have a acceptance integral matrix. Returning negative infinity for corresponding component of evidence." << endl;
		}

		// finally we calculate the probability of single waves being negligible and
		// take these reults into account
		const unsigned int nwaves = nmbWaves();
		double logprob = 0;
		for (unsigned int iwaves = 0; iwaves < nwaves; ++iwaves) {
			// exclude flat wave
			if (waveName(iwaves) == "flat")
				continue;

			// check that this wave is not amongst the thresholded waves
			if (thrWaveIndices.count(iwaves) > 0)
				continue;

			const double val = intensity   (iwaves);
			const double err = intensityErr(iwaves);
			// P(val>0); (assuming gaussian...) dirty!
			// require 3simga significance!
			const double prob = ROOT::Math::normal_cdf_c(5. * err, err, val);
			logprob += TMath::Log(prob);
		}
		components[3] = logprob;
	} else {
		printWarn << "fitResult does not have a valid error matrix. Returning negative infinity for corresponding components of evidence." << endl;
	}

	return components;
}


int
fitResult::waveIndex(const string& waveName) const
{
	int index = -1;
	for (unsigned int i = 0; i < nmbWaves(); ++i)
		if (waveName == this->waveName(i)) {
			index = i;
			break;  // assumes that waves have unique names
		}
	if (index == -1)
		printWarn << "could not find any wave named '" << waveName << "'." << endl;
	return index;
}


int
fitResult::prodAmpIndex(const string& prodAmpName) const
{
	int index = -1;
	for (unsigned int i = 0; i < nmbProdAmps(); ++i)
		if (prodAmpName == _prodAmpNames[i]) {
			index = i;
			break;  // assumes that production amplitudes have unique names
		}
	if (index == -1)
		printWarn << "could not find any production amplitude named '" << prodAmpName << "'." << endl;
	return index;
}


/// returns fit parameter value by parameter name
double
fitResult::fitParameter(const string& parName) const
{
	// check if parameter corresponds to real or imaginary part of production amplitude
	TString    name(parName);
	const bool realPart = (name.Contains("RE") or name.Contains("flat"));
	// find corresponding production amplitude
	if (realPart)
		name.ReplaceAll("_RE", "");
	else
		name.ReplaceAll("_IM", "");
	const int index = prodAmpIndex(name.Data());
	if (index >= 0) {
		if (realPart)
			return prodAmp(index).real();
		else
			return prodAmp(index).imag();
	}
	return 0;  // not found
}


/// \brief constructs 2n x 2n covariance matrix of production amplitudes specified by index list
// where n is the number of amplitudes
// layout:
//         cov(A0.re, A0.re)        cov(A0.re, A0.im)        ...  cov(A0.re, A(n - 1).re)        cov(A0.re, A(n - 1).im)
//         cov(A0.im, A0.re)        cov(A0.im, A0.im)        ...  cov(A0.im, A(n - 1).re)        cov(A0.im, A(n - 1).im)
//                  .                        .                             .                              .
//                  .                        .                             .                              .
//                  .                        .                             .                              .
//         cov(A(n - 1).re, A0.re)  cov(A(n - 1).re, A0.im)  ...  cov(A(n - 1).re, A(n - 1).re)  cov(A(n - 1).re, A(n - 1).im)
//         cov(A(n - 1).im, A0.re)  cov(A(n - 1).im, A0.im)  ...  cov(A(n - 1).im, A(n - 1).re)  cov(A(n - 1).im, A(n - 1).im)
// !!! possible optimization: exploit symmetry of cov matrix
TMatrixT<double>
fitResult::prodAmpCov(const vector<unsigned int>& prodAmpIndices) const
{
	const unsigned int dim = 2 * prodAmpIndices.size();
	TMatrixT<double>   prodAmpCov(dim, dim);
	if (not covMatrixValid()) {
		printWarn << "fitResult does not have a valid error matrix. Returning zero covariance matrix." << endl;
		return prodAmpCov;
	}

	// get corresponding indices for parameter covariances
	vector<int> parCovIndices;
	for (unsigned int i = 0; i < prodAmpIndices.size(); ++i) {
		parCovIndices.push_back(_fitParCovMatrixIndices[prodAmpIndices[i]].first);   // real part
		parCovIndices.push_back(_fitParCovMatrixIndices[prodAmpIndices[i]].second);  // imaginary part
	}

	// build covariance matrix
	for (unsigned int row = 0; row < dim; ++row)
		for (unsigned int col = 0; col < dim; ++col) {
			const int i = parCovIndices[row];
			const int j = parCovIndices[col];
			if ((i >= 0) and (j >= 0))
				prodAmpCov[row][col] = _fitParCovMatrix(i, j);
		}
	return prodAmpCov;
}


/// finds wave indices for production amplitues A and B and returns the normalization integral of the two waves
complex<double>
fitResult::normIntegralForProdAmp(const unsigned int prodAmpIndexA,
                                  const unsigned int prodAmpIndexB) const
{
	// treat special case of flat wave which has no normalization integral
	const bool flatWaveA = (prodAmpName(prodAmpIndexA).find("flat") != string::npos);
	const bool flatWaveB = (prodAmpName(prodAmpIndexB).find("flat") != string::npos);
	if (flatWaveA and flatWaveB)
		return 1;
	else if (flatWaveA or flatWaveB)
		return 0;
	else {
		map<int, int>::const_iterator indexA = _normIntIndexMap.find(prodAmpIndexA);
		map<int, int>::const_iterator indexB = _normIntIndexMap.find(prodAmpIndexB);
		if ((indexA == _normIntIndexMap.end()) or (indexB == _normIntIndexMap.end())) {
			printWarn << "amplitude index " << prodAmpIndexA << " or " << prodAmpIndexB
			          << " is out of bound. returning 0." << endl;
			return 0;
		}
		return normIntegral(indexA->second, indexB->second);
	}
}


/// \brief calculates spin density matrix element for waves A and B
///
/// \f[ \rho_{AB} = \sum_r V_{Ar} V_{Br}^* \f]
complex<double>
fitResult::spinDensityMatrixElem(const unsigned int waveIndexA,
                                 const unsigned int waveIndexB) const
{
	// spin density matrix element is 0 if the two waves have different
	// reflectivities
	if (partialWaveFitHelper::getReflectivity(waveName(waveIndexA))
	    != partialWaveFitHelper::getReflectivity(waveName(waveIndexB))) {
		return 0;
	}

	// get pairs of amplitude indices with the same rank for waves A and B
	const vector<pair<unsigned int, unsigned int> > prodAmpIndexPairs
		= prodAmpIndexPairsForWaves(waveIndexA, waveIndexB);
	if (prodAmpIndexPairs.size() == 0)
		return 0;

	// sum up amplitude products
	complex<double> spinDens = 0;
	for (unsigned int i = 0; i < prodAmpIndexPairs.size(); ++i)
		spinDens += prodAmp(prodAmpIndexPairs[i].first) * conj(prodAmp(prodAmpIndexPairs[i].second));
	return (double)normNmbEvents() * spinDens;
}


/// \brief calculates covariance matrix of spin density matrix element for waves A and B
///
/// \f[ \rho_{AB} = \sum_r V_{Ar} V_{Br}^* \f]
// !!! possible optimization: make special case for waveIndexA == waveIndexB
TMatrixT<double>
fitResult::spinDensityMatrixElemCov(const unsigned int waveIndexA,
                                    const unsigned int waveIndexB) const
{
	// covariance matrix has to be present
	if (not covMatrixValid()) {
		printWarn << "fitResult does not have a valid error matrix. Returning zero covariance matrix for spin-density matrix element." << endl;
		const TMatrixT<double> spinDensCov(2, 2);
		return spinDensCov;
	}

	// spin density matrix element is 0 if the two waves have different
	// reflectivities
	if (partialWaveFitHelper::getReflectivity(waveName(waveIndexA))
	    != partialWaveFitHelper::getReflectivity(waveName(waveIndexB))) {
		const TMatrixT<double> spinDensCov(2, 2);
		return spinDensCov;
	}

	// get pairs of amplitude indices with the same rank for waves A and B
	const vector<pair<unsigned int, unsigned int> > prodAmpIndexPairs
		= prodAmpIndexPairsForWaves(waveIndexA, waveIndexB);
	if (prodAmpIndexPairs.size() == 0) {
		const TMatrixT<double> spinDensCov(2, 2);
		return spinDensCov;
	}

	// build covariance matrix for amplitudes
	const TMatrixT<double> prodAmpCov = this->prodAmpCov(prodAmpIndexPairs);

	// build Jacobian for rho_AB, which is a 2 x 4m matrix composed of 2m sub-Jacobians:
	// J = (JA0, ..., JA(m - 1), JB0, ..., JB(m - 1))
	// m is the number of production amplitudes for waves A and B that have the same rank
	const unsigned int dim = prodAmpCov.GetNcols();
	TMatrixT<double>   jacobian(2, dim);
	// build m sub-Jacobians for d rho_AB / d V_Ar = M(V_Br^*)
	for (unsigned int i = 0; i < prodAmpIndexPairs.size(); ++i) {
		const unsigned int     ampIndexB    = prodAmpIndexPairs[i].second;
		const TMatrixT<double> subJacobianA = matrixRepr(conj(prodAmp(ampIndexB)));
		jacobian.SetSub(0, 2 * i, subJacobianA);
	}
	// build m sub-Jacobians for d rho_AB / d V_Br = M(V_Ar) {{1,  0},
	//                                                        {0, -1}}
	TMatrixT<double> M(2, 2);  // complex conjugation of V_Br is non-analytic operation
	M[0][0] =  1;
	M[1][1] = -1;
	const unsigned int colOffset = 2 * prodAmpIndexPairs.size();
	for (unsigned int i = 0; i < prodAmpIndexPairs.size(); ++i) {
		const unsigned int ampIndexA    = prodAmpIndexPairs[i].first;
		TMatrixT<double>   subJacobianB = matrixRepr(prodAmp(ampIndexA));
		subJacobianB *= M;
		jacobian.SetSub(0, colOffset + 2 * i, subJacobianB);
	}

	// calculate spin density covariance matrix
	// cov(rho_AB) = J cov(V_A0, ..., V_A(m - 1), V_B0, ..., V_B(m - 1)) J^T
	const TMatrixT<double> jacobianT(TMatrixT<double>::kTransposed, jacobian);
	// binary operations are unavoidable, since matrices are not squared
	// !!! possible optimaztion: use special TMatrixT constructors to perform the multiplication
	const TMatrixT<double> prodAmpCovJT = prodAmpCov * jacobianT;
	const TMatrixT<double> spinDensCov  = jacobian   * prodAmpCovJT;
	return std::pow((double)normNmbEvents(), 2.0) * spinDensCov;
}


/// calculates phase difference between wave A and wave B
double
fitResult::phase(const unsigned int waveIndexA,
                 const unsigned int waveIndexB) const
{
	if (waveIndexA == waveIndexB)
		return 0;
	return arg(spinDensityMatrixElem(waveIndexA, waveIndexB)) * TMath::RadToDeg();
}


/// calculates error of phase difference between wave A and wave B
double
fitResult::phaseErr(const unsigned int waveIndexA,
                    const unsigned int waveIndexB) const
{
	if (not covMatrixValid()) {
		printWarn << "fitResult does not have a valid error matrix. Returning zero error for phase." << endl;
		return 0;
	}

	if (waveIndexA == waveIndexB)
		return 0;

	// construct Jacobian for phi_AB = +- arctan(Im[rho_AB] / Re[rho_AB])
	const complex<double> spinDens = spinDensityMatrixElem(waveIndexA, waveIndexB);
	TMatrixT<double>      jacobian(1, 2);  // phase is real valued function, so J has only one row
	if (norm(spinDens) != 0) {
		jacobian[0][0] = -spinDens.imag() / norm(spinDens);
		jacobian[0][1] =  spinDens.real() / norm(spinDens);
	}

	// calculate variance
	const double phaseVariance = realValVariance(waveIndexA, waveIndexB, jacobian);
	return sqrt(phaseVariance) * TMath::RadToDeg();
}


/// calculates coherence of wave A and wave B
double
fitResult::coherence(const unsigned int waveIndexA,
                     const unsigned int waveIndexB) const
{
	const double          rhoAA = spinDensityMatrixElem(waveIndexA, waveIndexA).real();  // rho_AA is real by definition
	const double          rhoBB = spinDensityMatrixElem(waveIndexB, waveIndexB).real();  // rho_BB is real by definition
	const complex<double> rhoAB = spinDensityMatrixElem(waveIndexA, waveIndexB);
	return sqrt(std::norm(rhoAB) / (rhoAA * rhoBB));
}


/// calculates error of coherence of wave A and wave B
double
fitResult::coherenceErr(const unsigned int waveIndexA,
                        const unsigned int waveIndexB) const
{
	if (not covMatrixValid()) {
		printWarn << "fitResult does not have a valid error matrix. Returning zero error for coherence." << endl;
		return 0;
	}

	// get amplitude indices for waves A and B
	const vector<unsigned int> prodAmpIndicesA = prodAmpIndicesForWave(waveIndexA);
	const vector<unsigned int> prodAmpIndicesB = prodAmpIndicesForWave(waveIndexB);
	if ((prodAmpIndicesA.size() == 0) or (prodAmpIndicesB.size() == 0))
		return 0;

	// build Jacobian for coherence, which is a 1 x 2(n + m) matrix composed of (n + m) sub-Jacobians:
	// J = (JA_0, ..., JA_{n - 1}, JB_0, ..., JB_{m - 1})
	const unsigned int nmbAmps = prodAmpIndicesA.size() + prodAmpIndicesB.size();
	TMatrixT<double>   jacobian(1, 2 * nmbAmps);
	// precalculate some variables
	const double          rhoAA     = spinDensityMatrixElem(waveIndexA, waveIndexA).real();  // rho_AA is real by definition
	const double          rhoBB     = spinDensityMatrixElem(waveIndexB, waveIndexB).real();  // rho_BB is real by definition
	const complex<double> rhoAB     = spinDensityMatrixElem(waveIndexA, waveIndexB);
	const double          rhoABRe   = rhoAB.real();
	const double          rhoABIm   = rhoAB.imag();
	const double          rhoABNorm = std::norm(rhoAB);
	const double          coh       = sqrt(rhoABNorm / (rhoAA * rhoBB));
	if (coh == 0)
		return 0;
	// build m sub-Jacobians for JA_r = coh_AB / d V_Ar
	for (unsigned int i = 0; i < prodAmpIndicesA.size(); ++i) {
		const unsigned int    prodAmpIndexA = prodAmpIndicesA[i];
		const complex<double> prodAmpA      = prodAmp(prodAmpIndexA);
		const int             prodAmpRankA  = rankOfProdAmp(prodAmpIndexA);
		// find production amplitude of wave B with same rank
		complex<double> prodAmpB = 0;
		for (unsigned int j = 0; j < prodAmpIndicesB.size(); ++j) {
			const unsigned int prodAmpIndexB = prodAmpIndicesB[j];
			if (rankOfProdAmp(prodAmpIndexB) == prodAmpRankA) {
				prodAmpB = prodAmp(prodAmpIndexB);
				break;
			}
		}
		jacobian[0][2 * i    ] = rhoABRe * prodAmpB.real() - rhoABIm * prodAmpB.imag() - (rhoABNorm / rhoAA) * prodAmpA.real();
		jacobian[0][2 * i + 1] = rhoABRe * prodAmpB.imag() + rhoABIm * prodAmpB.real() - (rhoABNorm / rhoAA) * prodAmpA.imag();
	}
	// !!! possible optimization: join the loops for JA_r and JB_r
	// build m sub-Jacobians for JB_r = d coh_AB / d V_Br
	const unsigned int colOffset = 2 * prodAmpIndicesA.size();
	for (unsigned int i = 0; i < prodAmpIndicesB.size(); ++i) {
		const unsigned int    prodAmpIndexB = prodAmpIndicesB[i];
		const complex<double> prodAmpB      = prodAmp(prodAmpIndexB);
		const int             prodAmpRankB  = rankOfProdAmp(prodAmpIndexB);
		// find production amplitude of wave A with same rank
		complex<double> prodAmpA = 0;
		for (unsigned int j = 0; j < prodAmpIndicesA.size(); ++j) {
			const unsigned int prodAmpIndexA = prodAmpIndicesA[j];
			if (rankOfProdAmp(prodAmpIndexA) == prodAmpRankB) {
				prodAmpA = prodAmp(prodAmpIndexA);
				break;
			}
		}
		jacobian[0][colOffset + 2 * i    ] = rhoABRe * prodAmpA.real() + rhoABIm * prodAmpA.imag() - (rhoABNorm / rhoBB) * prodAmpB.real();
		jacobian[0][colOffset + 2 * i + 1] = rhoABRe * prodAmpA.imag() - rhoABIm * prodAmpA.real() - (rhoABNorm / rhoBB) * prodAmpB.imag();
	}
	jacobian *= 1 / (coh * rhoAA * rhoBB);

	// build covariance matrix for amplitudes and calculate coherence covariance matrix
	const TMatrixT<double> prodAmpCov   = this->prodAmpCov(prodAmpIndicesA, prodAmpIndicesB);  // 2(n + m) x 2(n + m) matrix
	const TMatrixT<double> jacobianT(TMatrixT<double>::kTransposed, jacobian);                 // 2(n + m) x        1 matrix
	const TMatrixT<double> prodAmpCovJT = prodAmpCov * jacobianT;                              // 2(n + m) x        1 matrix
	const TMatrixT<double> cohCov       = jacobian   * prodAmpCovJT;                           //        1 x        1 matrix
	return sqrt(cohCov[0][0]);
}


/// calculates overlap of wave A and wave B
double
fitResult::overlap(const unsigned int waveIndexA,
                   const unsigned int waveIndexB) const
{
	const complex<double> spinDens = spinDensityMatrixElem(waveIndexA, waveIndexB);
	const complex<double> normInt  = normIntegral         (waveIndexA, waveIndexB);
	return 2 * (spinDens * normInt).real();
}


/// calculates error of overlap of wave A and wave B
double
fitResult::overlapErr(const unsigned int waveIndexA,
                      const unsigned int waveIndexB) const
{
	if (not covMatrixValid()) {
		printWarn << "fitResult does not have a valid error matrix. Returning zero error for overlap." << endl;
		return 0;
	}

	const complex<double> normInt = normIntegral(waveIndexA, waveIndexB);
	TMatrixT<double> jacobian(1, 2);  // overlap is real valued function, so J has only one row
	jacobian[0][0] =  2 * normInt.real();
	jacobian[0][1] = -2 * normInt.imag();
	const double overlapVariance = realValVariance(waveIndexA, waveIndexB, jacobian);
	return sqrt(overlapVariance);
}


/// \brief calculates intensity for set of waves matching name pattern
///
/// int = sum_i int(i) + sum_i sum_{j < i} overlap(i, j)
double
fitResult::intensity(const string& waveNamePattern) const
{
	const vector<unsigned int> waveIndices = waveIndicesMatchingPattern(waveNamePattern);
	const unsigned int         nmbWaves    = waveIndices.size();
	if (nmbWaves == 0) {
		return 0;
	}
	if (nmbWaves == 1) {
		return intensity(waveIndices[0]);
	}

	double intensity = 0;
	for (unsigned int i = 0; i < nmbWaves; ++i) {
		intensity += this->intensity(waveIndices[i]);
		for (unsigned int j = 0; j < i; ++j) {
			intensity += overlap(waveIndices[i], waveIndices[j]);
		}
	}
	return intensity;
}


/// \brief calculates error of intensity of a set of waves matching name pattern
///
/// error calculation is performed on amplitude level using: int = sum_ij Norm_ij sum_r A_ir A_jr*
double
fitResult::intensityErr(const string& waveNamePattern) const
{
	if (not covMatrixValid()) {
		printWarn << "fitResult does not have a valid error matrix. Returning zero error for intensity." << endl;
		return 0.;
	}

	// get amplitudes that correspond to wave name pattern
	const vector<unsigned int> waveIndices = waveIndicesMatchingPattern(waveNamePattern);
	const unsigned int         nmbWaves    = waveIndices.size();
	if (nmbWaves == 0) {
		return 0;
	}
	if (nmbWaves == 1) {
		return intensityErr(waveIndices[0]);
	}

	vector<unsigned int> prodAmpIndices;
	for (unsigned int i = 0; i < nmbWaves; ++i) {
		const vector<unsigned int> prodAmpIndicesWave = prodAmpIndicesForWave(waveIndices[i]);
		prodAmpIndices.insert(prodAmpIndices.end(), prodAmpIndicesWave.begin(), prodAmpIndicesWave.end());
	}
	const unsigned int nmbAmps = prodAmpIndices.size();

	// build Jacobian for intensity, which is a 1 x 2n matrix composed of n sub-Jacobians:
	// J = (JA_0, ..., JA_{n - 1}), where n is the number of production amplitudes
	TMatrixT<double> jacobian(1, 2 * nmbAmps);
	for (unsigned int i = 0; i < nmbAmps; ++i) {
		// build sub-Jacobian for each amplitude; intensity is real valued function, so J has only one row
		// JA_ir = 2 * sum_j (A_jr Norm_ji)
		complex<double>  ampNorm     = 0;  // sum_j (A_jr Norm_ji)
		const int        currentRefl = partialWaveFitHelper::getReflectivity(waveNameForProdAmp(prodAmpIndices[i]));
		const int        currentRank = rankOfProdAmp(prodAmpIndices[i]);
		for (unsigned int j = 0; j < nmbAmps; ++j) {
			if (partialWaveFitHelper::getReflectivity(waveNameForProdAmp(prodAmpIndices[j])) != currentRefl)
				continue;
			if (rankOfProdAmp(prodAmpIndices[j]) != currentRank)
				continue;
			ampNorm += prodAmp(prodAmpIndices[j]) * normIntegralForProdAmp(prodAmpIndices[j], prodAmpIndices[i]);  // order of indices is essential
		}
		jacobian[0][2 * i    ] = ampNorm.real();
		jacobian[0][2 * i + 1] = ampNorm.imag();
	}
	jacobian *= 2;

	const TMatrixT<double> prodAmpCov   = this->prodAmpCov(prodAmpIndices);     // 2n x 2n matrix
	const TMatrixT<double> jacobianT(TMatrixT<double>::kTransposed, jacobian);  // 2n x  1 matrix
	const TMatrixT<double> prodAmpCovJT = prodAmpCov * jacobianT;               // 2n x  1 matrix
	const TMatrixT<double> intensityCov = jacobian * prodAmpCovJT;              //  1 x  1 matrix
	return (double)normNmbEvents() * sqrt(intensityCov[0][0]);
}


std::map<rpwa::multibinBoundariesType, std::list<rpwa::fitResult> >
rpwa::getFitResultsFromFilesInMultibins(
                                        const std::vector<std::string>& fileNames,
                                        const std::string& treeName,
                                        const std::string& branchName,
                                        const bool onlyBestResultInMultibin,
                                        const bool stripMatricesFromNotBestResults,
                                        const bool onlyConvergedResults) {
	std::map<rpwa::multibinBoundariesType, std::list<rpwa::fitResult> > fitResultsInMultibins;
	std::unique_ptr<fitResult> inputFitResult(new fitResult);
	fitResult* inputFitResultPtr = inputFitResult.get();

	// get trees
	std::vector<TFile*> files;
	std::vector<TTree*> trees;
	for (const auto& fileName : fileNames) {
		TFile* file = TFile::Open(fileName.c_str());
		if (file == nullptr or file->IsZombie()) {
			printErr<< "Cannot open rootfile '" << fileName << "'" << std::endl;
			return fitResultsInMultibins;
		}
		TTree* tree;
		file->GetObject(treeName.c_str(), tree);
		if (tree == nullptr) {
			printErr<< "Cannot find tree '" << treeName << "' in rootfile '" << fileName << "'" << std::endl;
			return fitResultsInMultibins;
		}
		tree->SetBranchAddress(branchName.c_str(), &inputFitResultPtr);
		trees.push_back(tree);
		files.push_back(file);

	}

	// find the best result in each multibin and store all others if necessary
	std::set<rpwa::multibinBoundariesType> allMultibinBoundaries;
	std::map<rpwa::multibinBoundariesType, double> negLogLikeOfBestResult;
	std::map<rpwa::multibinBoundariesType, std::pair<size_t, int> > fileTreeEntryOfBestResult;
	std::map<rpwa::multibinBoundariesType, std::list<rpwa::fitResult>::iterator > iterOfBestResult;
	std::map<rpwa::multibinBoundariesType, double> negLogLikeOfBestConvergedResult;
	std::map<rpwa::multibinBoundariesType, std::pair<size_t, int> > fileTreeEntryOfBestConvergedResult;
	std::map<rpwa::multibinBoundariesType, std::list<rpwa::fitResult>::iterator > iterOfBestConvergedResult;
	for (size_t i = 0; i < trees.size(); ++i) {
		printInfo << "load fit results from file '" << fileNames[i] << "'" << std::endl;

		TTree* tree = trees[i];
		for (int j = 0; j < tree->GetEntries(); ++j) {
			tree->GetEntry(j);
			if( onlyConvergedResults and not inputFitResult->converged())
				continue;

			const rpwa::multibinBoundariesType& multibinBoundaries = inputFitResult->multibinBoundaries();
			allMultibinBoundaries.insert(multibinBoundaries);
			std::list<rpwa::fitResult>::iterator it;

			if (not onlyBestResultInMultibin) {
				it = fitResultsInMultibins[multibinBoundaries].insert(fitResultsInMultibins[multibinBoundaries].end(), fitResult());
				it->fill(*inputFitResult, not stripMatricesFromNotBestResults, not stripMatricesFromNotBestResults);

				if(inputFitResult->converged()){
					if (negLogLikeOfBestConvergedResult.find(multibinBoundaries) == negLogLikeOfBestConvergedResult.end() or // first result of this bin
					    inputFitResult->logLikelihood() < negLogLikeOfBestConvergedResult[multibinBoundaries]) { // found new best result
						negLogLikeOfBestConvergedResult[multibinBoundaries] = inputFitResult->logLikelihood();
						fileTreeEntryOfBestConvergedResult[multibinBoundaries] = std::make_pair(i, j);
						iterOfBestConvergedResult[multibinBoundaries] = it;
					}
				}
			}
			if (negLogLikeOfBestResult.find(multibinBoundaries) == negLogLikeOfBestResult.end() or // first result of this bin
			    inputFitResult->logLikelihood() < negLogLikeOfBestResult[multibinBoundaries]) { // found new best result
				negLogLikeOfBestResult[multibinBoundaries] = inputFitResult->logLikelihood();
				fileTreeEntryOfBestResult[multibinBoundaries] = std::make_pair(i, j);
				iterOfBestResult[multibinBoundaries] = it;
			}
		}
	}


	// finally add the result with the best likelihood, if not yet stored
	if (onlyBestResultInMultibin) {
		for (const auto& multibinBoundaries : allMultibinBoundaries) {
			// load best result
			trees[fileTreeEntryOfBestResult[multibinBoundaries].first]->GetEntry(fileTreeEntryOfBestResult[multibinBoundaries].second);
			std::list<fitResult>& fitResults = fitResultsInMultibins[multibinBoundaries];
			assert(fitResults.size() == 0);
			fitResults.push_back(*inputFitResult);
		}
	} else if (stripMatricesFromNotBestResults) {
		for (const auto& multibinBoundaries : allMultibinBoundaries) {
			// load best result
			trees[fileTreeEntryOfBestResult[multibinBoundaries].first]->GetEntry(fileTreeEntryOfBestResult[multibinBoundaries].second);
			assert(negLogLikeOfBestResult[multibinBoundaries] == inputFitResult->logLikelihood());
			assert(iterOfBestResult[multibinBoundaries]->logLikelihood() == inputFitResult->logLikelihood());
			// remove the best result without matrices and add the full result
			iterOfBestResult[multibinBoundaries]->fill(*inputFitResult);

			// if the best one is not converged and we want to store all and we want to strip the matrices from all except the best ones (converged and not converged)
			if( negLogLikeOfBestConvergedResult.find(multibinBoundaries) != negLogLikeOfBestConvergedResult.end()){ // found at least one converged result
				if (not iterOfBestResult[multibinBoundaries]->converged()){
					// load best converged result
					trees[fileTreeEntryOfBestConvergedResult[multibinBoundaries].first]->GetEntry(fileTreeEntryOfBestConvergedResult[multibinBoundaries].second);
					iterOfBestConvergedResult[multibinBoundaries]->fill(*inputFitResult);
				} else {
					assert(negLogLikeOfBestResult[multibinBoundaries] == negLogLikeOfBestConvergedResult[multibinBoundaries]);
				}
			}
		}
	}

	// sort fit results by neg log-likelihood
	if (not onlyBestResultInMultibin) {
		for (auto& elem : fitResultsInMultibins) {
			const multibinBoundariesType& multibinBoundaries = elem.first;
			std::list<fitResult>& fitResults = elem.second;
			// sort preserves the order of equal elements and we set a new best result only if the result < as the old best one,
			// thus the best result is the first best result. Therefore, in the sorted list of the first result is the one without striped matrices
			fitResults.sort([] (const rpwa::fitResult& a, const rpwa::fitResult& b) -> bool {return a.logLikelihood() < b.logLikelihood();});

			assert(fitResults.front().logLikelihood() == negLogLikeOfBestResult[multibinBoundaries]);
		}
	}

	// close files
	for (const auto& file : files) {
		file->Close();
		delete file;
	}
	return fitResultsInMultibins;
}
