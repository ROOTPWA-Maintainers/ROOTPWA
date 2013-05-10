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

#include "Math/SpecFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"
#include "TDecompChol.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TRandom3.h"

// PWA2000 classes
#include "integral.h"

#include "fitResult.h"


using namespace std;
using namespace rpwa;


ClassImp(fitResult);


fitResult::fitResult()
	: _nmbEvents     (0),
	  _normNmbEvents (0),
	  _massBinCenter (0),
	  _logLikelihood (0),
	  _rank          (0),
	  _covMatrixValid(false),
	  _converged     (false),
	  _hasHessian    (false)
{ }


fitResult::fitResult(const TFitBin& fitBin)
	: _nmbEvents             (fitBin.nmbEvents()),
	  _normNmbEvents         (fitBin.normNmbEvents()),
	  _massBinCenter         (fitBin.massBinCenter()),
	  _logLikelihood         (fitBin.logLikelihood()),
	  _rank                  (fitBin.rank()),
	  _covMatrixValid        (fitBin.fitParCovMatrixValid()),
	  _fitParCovMatrix       (fitBin.fitParCovMatrix()),
	  _fitParCovMatrixIndices(fitBin.fitParCovIndices()),
	  _normIntegral          (fitBin.normIntegral()),
	  _normIntIndexMap       (fitBin.prodAmpIndexMap()),
	  _converged             (true),
	  _hasHessian            (false)
{
	_prodAmps = fitBin.prodAmps();
	{
		const unsigned int nmbAmpNames = fitBin.prodAmpNames().size();
		_prodAmpNames.resize(nmbAmpNames, "");
		for (unsigned int i = 0; i < nmbAmpNames; ++i)
			_prodAmpNames[i] = fitBin.prodAmpNames()[i].Data();
	}
	{
		const unsigned int nmbWaveNames = fitBin.waveNames().size();
		_waveNames.resize(nmbWaveNames, "");
		for (unsigned int i = 0; i < nmbWaveNames; ++i)
			_waveNames[i] = fitBin.waveNames()[i].Data();
	}
}


fitResult::fitResult(const fitResult& result)
	: _nmbEvents             (result.nmbEvents()),
	  _normNmbEvents         (result.normNmbEvents()),
	  _massBinCenter         (result.massBinCenter()),
	  _logLikelihood         (result.logLikelihood()),
	  _rank                  (result.rank()),
	  _prodAmps              (result.prodAmps()),
	  _prodAmpNames          (result.prodAmpNames()),
	  _waveNames             (result.waveNames()),
	  _covMatrixValid        (result.covMatrixValid()),
	  _fitParCovMatrix       (result.fitParCovMatrix()),
	  _fitParCovMatrixIndices(result.fitParCovIndices()),
	  _normIntegral          (result.normIntegralMatrix()),
	  _normIntIndexMap       (result.normIntIndexMap()),
	  _phaseSpaceIntegral    (result._phaseSpaceIntegral),
	  _converged             (result._converged),
	  _hasHessian            (result._hasHessian)
{ }


// enable copying from TFitResult for older ROOT versions
#ifdef USE_TFITRESULT
fitResult::fitResult(const TFitResult& result)
	: _nmbEvents             (result.nmbEvents()),
	  _normNmbEvents         (result.normNmbEvents()),
	  _massBinCenter         (result.massBinCenter()),
	  _logLikelihood         (result.logLikelihood()),
	  _rank                  (result.rank()),
	  _prodAmps              (result.prodAmps()),
	  _prodAmpNames          (result.prodAmpNames()),
	  _waveNames             (result.waveNames()),
	  _covMatrixValid        (result.covMatrixValid()),
	  _fitParCovMatrix       (result.fitParCovMatrix()),
	  _fitParCovMatrixIndices(result.fitParCovIndices()),
	  _normIntegral          (result.normIntegralMatrix()),
	  _normIntIndexMap       (result.normIntIndexMap())
{ }
#endif


fitResult::~fitResult()
{ }


fitResult*
fitResult::cloneVar(){
  // copy complete info
  fitResult* result = new fitResult(*this);

  // Cholesky decomposition (hoepfully this will not slow us down too much
  // if it does we have to cache
  cerr << "Starting Cholesky Decomposition" << endl;
  TDecompChol decomp(_fitParCovMatrix);
  decomp.Decompose();
  TMatrixT<Double_t> C(decomp.GetU());
  cerr << "...Done" << endl;

  unsigned int npar=C.GetNrows();
  TMatrixD x(npar,1);
  // generate npar independent random numbers
  for(unsigned int ipar=0;ipar<npar;++ipar){
    x(ipar,0)=gRandom->Gaus();
  }

  // this tells us how to permute the parameters taking into account all
  // correlations in the covariance matrix
  TMatrixD y(C,TMatrixD::kMult,x);

  // now we need mapping from parameters to production amp re and im
  // loop over production amps
  unsigned int nt=_prodAmps.size();
  for(unsigned int it=0;it<nt;++it){
    Int_t jre=_fitParCovMatrixIndices[it].first;
    Int_t jim=_fitParCovMatrixIndices[it].second;
    double re=result->_prodAmps[it].Re(); double im=result->_prodAmps[it].Im();
    if(jre>-1)re+=y(jre,0);
    if(jim>-1)im+=y(jim,0);
    result->_prodAmps[it]=TComplex(re,im);
    cerr << "Modifying amp " << it << " by (" << (jre>-1 ? y(jre,0) : 0)
         << "," << (jim>-1 ? y(jim,0) : 0) <<")" << endl;
  }
  return result;
}


/// \brief calculates the model evidence using Rafters occam factor method
///
/// \f[  \ln P(\mathrm{Data}|M_k)\approx \ln P(\mathrm{Data}|A_{\mathrm{ML}}^k,M_k) + \ln{\frac{\sqrt{(2\pi)^d|\mathbf{C}_{A|D}|}}{\sum_i\frac{\pi}{\Psi_{ii}}}} \f]
/// This uses a flat prior in each parameter (see note on Model Selection) such
/// that no single wave can have more than the total intensity measured
double
fitResult::evidence() const
{
	// find the thresholded production amplitudes, assume those are the
	// ones with imaginary and real part equal to zero
	set<unsigned int> thrProdAmpIndices;
	for (unsigned int i = 0; i < nmbProdAmps(); ++i) {
		// printDebug << "production amplitude[" << i << "] = '" << _prodAmpNames[i]
		//            << "' = " << prodAmp(i) << endl;
		if (prodAmp(i) == 0.) {
			// printDebug << "IS ZERO!!!" << endl;
			thrProdAmpIndices.insert(i);
		}
	}

	// from the list of thresholded production amplitudes get the list of
	// thresholded waves, this assumes that if a production amplitude is
	// thresholded, it is thresholded in the same way for all ranks
	set<unsigned int> thrWaveIndices;
	for (set<unsigned int>::const_iterator it = thrProdAmpIndices.begin();
	     it != thrProdAmpIndices.end(); ++it) {
		const int index = waveIndex(waveNameForProdAmp(*it).Data());
		// printDebug << " wave[" << index << "] = '" << waveNameForProdAmp(*it).Data()
		//            << "' has zero production amp[" << *it << "] = '"
		//            << _prodAmpNames[*it] << "'" << endl;
		assert(index > 0);
		thrWaveIndices.insert(index);
	}

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
			// printDebug << "disregarding covariance row " << i << endl;
			continue;
		}
		row++;

		Int_t col = -1;
		for (Int_t j = 0; j < fitParCovMatrix().GetNcols(); ++j) {
			if (thrColAndRowIndices.count(j) > 0) {
				// printDebug << "disregarding covariance column " << j << endl;
				continue;
			}
			col++;

			thrCovMatrix(row, col) = fitParCovMatrix()(i, j);
		}
	}

	// REMOVE CONSTRAINT TO NUMBER OF EVENTS!
	const double l   = -logLikelihood();// - intensity(".*");
	const double det = thrCovMatrix.Determinant();
	const double d   = (double)thrCovMatrix.GetNcols();

	// parameter-volume after observing data
	const double lvad = 0.5 * (d * 1.837877066 + TMath::Log(det));

	// parameter volume prior to observing the data
	// n-Sphere:
	const double lva = TMath::Log(d) + 0.5 * (d * 1.144729886 + (d - 1) * TMath::Log(_nmbEvents))
		- ROOT::Math::lgamma(0.5 * d + 1);

	// finally we calculate the probability of single waves being negligible and
	// take these reults into account
	const unsigned int nwaves = nmbWaves();
	double logprob = 0;
	for (unsigned int iwaves = 0; iwaves < nwaves; ++iwaves) {
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

	return l + lvad - lva + logprob;
}


/// \brief calculates spin density matrix element for waves A and B
///
/// \f[ \rho_{AB} = \sum_r V_{Ar} V_{Br}^* \f]
complex<double>
fitResult::spinDensityMatrixElem(const unsigned int waveIndexA,
                                 const unsigned int waveIndexB) const
{
	// get pairs of amplitude indices with the same rank for waves A and B
	const vector<pair<unsigned int, unsigned int> > prodAmpIndexPairs
		= prodAmpIndexPairsForWaves(waveIndexA, waveIndexB);
	if (prodAmpIndexPairs.size() == 0)
		return 0;
	// sum up amplitude products
	complex<double> spinDens = 0;
	for (unsigned int i = 0; i < prodAmpIndexPairs.size(); ++i)
		spinDens += prodAmp(prodAmpIndexPairs[i].first) * conj(prodAmp(prodAmpIndexPairs[i].second));
	return spinDens;
}


/// returns fit parameter value by parameter name
double
fitResult::fitParameter(const string& parName) const
{
	// check if parameter corresponds to real or imaginary part of production amplitude
	TString    name(parName);
	const bool realPart = (name.Contains("RE") || name.Contains("flat"));
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
	if (!_covMatrixValid) {
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
			if ((i >= 0) && (j >= 0))
				prodAmpCov[row][col] = fitParameterCov(i, j);
		}
	return prodAmpCov;
}


/// \brief calculates covariance matrix of spin density matrix element for waves A and B
///
/// rho_AB = sum_r V_Ar V_Br^*
// !!! possible optimization: make special case for waveIndexA == waveIndexB
TMatrixT<double>
fitResult::spinDensityMatrixElemCov(const unsigned int waveIndexA,
                                    const unsigned int waveIndexB) const
{
	// get pairs of amplitude indices with the same rank for waves A and B
	const vector<pair<unsigned int, unsigned int> > prodAmpIndexPairs
		= prodAmpIndexPairsForWaves(waveIndexA, waveIndexB);
	if (!_covMatrixValid || (prodAmpIndexPairs.size() == 0)) {
		TMatrixT<double> spinDensCov(2, 2);
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
	// build m sub-Jacobian for d rho_AB / d V_Br = M(V_Ar) {{1,  0},
	//                                                       {0, -1}}
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
	return spinDensCov;
}


/// \brief calculates intensity for set of waves matching name pattern
///
/// int = sum_i int(i) + sum_i sum_{j < i} overlap(i, j)
double
fitResult::intensity(const vector<unsigned int>& waveIndices) const
{
	double               intensity   = 0;
	for (unsigned int i = 0; i < waveIndices.size(); ++i) {
		intensity += this->intensity(waveIndices[i]);
		// cout << "    contribution from " << _waveNames[waveIndices[i]]
		//      << " = " << this->intensity(waveIndices[i]) << endl;
		for (unsigned int j = 0; j < i; ++j) {
			intensity += overlap(waveIndices[i], waveIndices[j]);
			// cout << "        overlap with " << _waveNames[waveIndices[j]]
			//      << " = " << overlap(waveIndices[i], waveIndices[j]) << endl;
		}
	}
	return intensity;
}


/// \brief calculates intensity for set of waves with indices from waveIndices
///
/// int = sum_i int(i) + sum_i sum_{j < i} overlap(i, j)
double
fitResult::intensity(const char* waveNamePattern) const
{
	vector<unsigned int> waveIndices = waveIndicesMatchingPattern(waveNamePattern);

	return intensity(waveIndices);
}


/// finds wave indices for production amplitues A and B and returns the normalization integral of the two waves
complex<double>
fitResult::normIntegralForProdAmp(const unsigned int prodAmpIndexA,
                                  const unsigned int prodAmpIndexB) const
{
	// treat special case of flat wave which has no normalization integral
	const bool flatWaveA = prodAmpName(prodAmpIndexA).Contains("flat");
	const bool flatWaveB = prodAmpName(prodAmpIndexB).Contains("flat");
	if (flatWaveA && flatWaveB)
		return 1;
	else if (flatWaveA || flatWaveB)
		return 0;
	else {
		map<int, int>::const_iterator indexA = _normIntIndexMap.find(prodAmpIndexA);
		map<int, int>::const_iterator indexB = _normIntIndexMap.find(prodAmpIndexB);
		if ((indexA == _normIntIndexMap.end()) || (indexB == _normIntIndexMap.end())) {
			printWarn << "Amplitude index " << prodAmpIndexA << " or " << prodAmpIndexB
			          << " is out of bound." << endl;
			return 0;
		}
		return normIntegral(indexA->second, indexB->second);
	}
}


/// \brief calculates error of intensity of a set of waves with indices from  waveIndices
///
/// error calculation is performed on amplitude level using: int = sum_ij Norm_ij sum_r A_ir A_jr*
double
fitResult::intensityErr(const std::vector<unsigned int>& prodAmpIndices) const
{
	// get amplitudes that correspond to wave name pattern
	const unsigned int         nmbAmps        = prodAmpIndices.size();
	if (!_covMatrixValid || (nmbAmps == 0))
		return 0;
	// build Jacobian for intensity, which is a 1 x 2n matrix composed of n sub-Jacobians:
	// J = (JA_0, ..., JA_{n - 1}), where n is the number of production amplitudes
	TMatrixT<double> jacobian(1, 2 * nmbAmps);
	for (unsigned int i = 0; i < nmbAmps; ++i) {
		// build sub-Jacobian for each amplitude; intensity is real valued function, so J has only one row
		// JA_ir = 2 * sum_j (A_jr Norm_ji)
		complex<double>  ampNorm     = 0;  // sum_j (A_jr Norm_ji)
		const int        currentRank = rankOfProdAmp(prodAmpIndices[i]);
		for (unsigned int j = 0; j < nmbAmps; ++j) {
			if (rankOfProdAmp(prodAmpIndices[j]) != currentRank)
				continue;
			ampNorm += prodAmp(prodAmpIndices[j]) * normIntegralForProdAmp(j, i);  // order of indices is essential
		}
		jacobian[0][2 * i    ] = ampNorm.real();
		jacobian[0][2 * i + 1] = ampNorm.imag();
	}
	jacobian *= 2;
	const TMatrixT<double> prodAmpCov   = this->prodAmpCov(prodAmpIndices);     // 2n x 2n matrix
	const TMatrixT<double> jacobianT(TMatrixT<double>::kTransposed, jacobian);  // 2n x  1 matrix
	const TMatrixT<double> prodAmpCovJT = prodAmpCov * jacobianT;               // 2n x  1 matrix
	const TMatrixT<double> intensityCov = jacobian * prodAmpCovJT;              //  1 x  1 matrix
	return sqrt(intensityCov[0][0]);
}

/// \brief calculates error of intensity of a set of waves matching name pattern
///
/// error calculation is performed on amplitude level using: int = sum_ij Norm_ij sum_r A_ir A_jr*
double
fitResult::intensityErr(const char* waveNamePattern) const
{
	// get amplitudes that correspond to wave name pattern
	const vector<unsigned int> prodAmpIndices = prodAmpIndicesMatchingPattern(waveNamePattern);
	return intensityErr(prodAmpIndices);
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
	if (!_covMatrixValid || (waveIndexA == waveIndexB))
		return 0;
	// construct Jacobian for phi_AB = +- arctan(Im[rho_AB] / Re[rho_AB])
	const complex<double> spinDens = spinDensityMatrixElem(waveIndexA, waveIndexB);
	TMatrixT<double>      jacobian(1, 2);  // phase is real valued function, so J has only one row
	{
		const double x = spinDens.real();
		const double y = spinDens.imag();
		if ((x != 0) || (y != 0)) {
			jacobian[0][0] = 1 / (x + y * y / x);
			jacobian[0][1] = -y / (x * x + y * y);
		}
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
	// get amplitude indices for waves A and B
	const vector<unsigned int> prodAmpIndices[2] = {prodAmpIndicesForWave(waveIndexA),
	                                                prodAmpIndicesForWave(waveIndexB)};
	if (!_covMatrixValid || (prodAmpIndices[0].size() == 0) || (prodAmpIndices[1].size() == 0))
		return 0;
	// build Jacobian for coherence, which is a 1 x 2(n + m) matrix composed of (n + m) sub-Jacobians:
	// J = (JA_0, ..., JA_{n - 1}, JB_0, ..., JB_{m - 1})
	const unsigned int nmbAmps = prodAmpIndices[0].size() + prodAmpIndices[1].size();
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
	for (unsigned int i = 0; i < prodAmpIndices[0].size(); ++i) {
		const unsigned int    prodAmpIndexA = prodAmpIndices[0][i];
		const complex<double> prodAmpA      = prodAmp(prodAmpIndexA);
		const int             prodAmpRankA  = rankOfProdAmp(prodAmpIndexA);
		// find production amplitude of wave B with same rank
		complex<double> prodAmpB = 0;
		for (unsigned int j = 0; j < prodAmpIndices[1].size(); ++j) {
			const unsigned int prodAmpIndexB = prodAmpIndices[1][j];
			if (rankOfProdAmp(prodAmpIndexB) == prodAmpRankA) {
				prodAmpB = prodAmp(prodAmpIndexB);
				break;
			}
		}
		jacobian[0][2 * i    ] = rhoABRe * prodAmpB.real() - rhoABIm * prodAmpB.imag() - (rhoABNorm / rhoAA) * prodAmpA.real();
		jacobian[0][2 * i + 1] = rhoABRe * prodAmpB.imag() + rhoABIm * prodAmpB.real() - (rhoABNorm / rhoAA) * prodAmpA.imag();
	}
	// !!! possible optimization: join the loops for JA_r and JB_r
	// build m sub-Jacobian for JB_r = d coh_AB / d V_Br
	const unsigned int colOffset = 2 * prodAmpIndices[0].size();
	for (unsigned int i = 0; i < prodAmpIndices[1].size(); ++i) {
		const unsigned int    prodAmpIndexB = prodAmpIndices[1][i];
		const complex<double> prodAmpB      = prodAmp(prodAmpIndexB);
		const int             prodAmpRankB  = rankOfProdAmp(prodAmpIndexB);
		// find production amplitude of wave A with same rank
		complex<double> prodAmpA = 0;
		for (unsigned int j = 0; j < prodAmpIndices[0].size(); ++j) {
			const unsigned int prodAmpIndexA = prodAmpIndices[0][j];
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
	const TMatrixT<double> prodAmpCov   = this->prodAmpCov(prodAmpIndices[0], prodAmpIndices[1]);  // 2(n + m) x 2(n + m) matrix
	const TMatrixT<double> jacobianT(TMatrixT<double>::kTransposed, jacobian);                     // 2(n + m) x        1 matrix
	const TMatrixT<double> prodAmpCovJT = prodAmpCov * jacobianT;                                  // 2(n + m) x        1 matrix
	const TMatrixT<double> cohCov       = jacobian   * prodAmpCovJT;                               //        1 x        1 matrix
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
	if (!_covMatrixValid)
		return 0;
	const complex<double> normInt = normIntegral(waveIndexA, waveIndexB);
	TMatrixT<double> jacobian(1, 2);  // overlap is real valued function, so J has only one row
	jacobian[0][0] =  2 * normInt.real();
	jacobian[0][1] = -2 * normInt.imag();
	const double overlapVariance = realValVariance(waveIndexA, waveIndexB, jacobian);
	return sqrt(overlapVariance);
}


void
fitResult::reset()
{
	_nmbEvents     = 0;
	_normNmbEvents = 0;
	_massBinCenter = 0;
	_logLikelihood = 0;
	_rank          = 0;
	_prodAmps.clear();
	_prodAmpNames.clear();
	_waveNames.clear();
	_covMatrixValid = false;
	_fitParCovMatrix.ResizeTo(0, 0);
	_fitParCovMatrixIndices.clear();
	_normIntegral.ResizeTo(0, 0);
	_normIntIndexMap.clear();
	_phaseSpaceIntegral.clear();
	_converged  = false;
	_hasHessian = false;
}


void
fitResult::fill
(const unsigned int              nmbEvents,               // number of events in bin
 const unsigned int              normNmbEvents,	          // number of events to normalize to
 const double                    massBinCenter,	          // center value of mass bin
 const double                    logLikelihood,	          // log(likelihood) at maximum
 const int                       rank,		                // rank of fit
 const vector<complex<double> >& prodAmps,	              // production amplitudes
 const vector<string>&           prodAmpNames,	          // names of production amplitudes used in fit
 const TMatrixT<double>&         fitParCovMatrix,         // covariance matrix of fit parameters
 const vector<pair<int, int> >&  fitParCovMatrixIndices,  // indices of fit parameters for real and imaginary part in covariance matrix matrix
 const TCMatrix&                 normIntegral,            // normalization integral matrix
 const vector<double>&           phaseSpaceIntegral,      // normalization integral over full phase space without acceptance
 const bool                      converged,
 const bool                      hasHessian)
{
	_converged     = converged;
	_hasHessian    = hasHessian;
	_nmbEvents     = nmbEvents;
	_normNmbEvents = normNmbEvents;
	_massBinCenter = massBinCenter;
	_logLikelihood = logLikelihood;
	_rank          = rank;
	_prodAmps.resize(prodAmps.size());
	for (unsigned int i = 0; i < prodAmps.size(); ++i)
		_prodAmps[i] = TComplex(prodAmps[i].real(), prodAmps[i].imag());
	_prodAmpNames  = prodAmpNames;
	_fitParCovMatrix.ResizeTo(fitParCovMatrix.GetNrows(), fitParCovMatrix.GetNcols());
	_fitParCovMatrix        = fitParCovMatrix;
	_fitParCovMatrixIndices = fitParCovMatrixIndices;
	// check whether there really is an error matrix
	if (!(fitParCovMatrix.GetNrows() == 0) && !(fitParCovMatrix.GetNcols() == 0))
		_covMatrixValid = true;
	else
		_covMatrixValid = false;
	_normIntegral.ResizeTo(normIntegral.nrows(), normIntegral.ncols());
	_normIntegral       = normIntegral;
	_phaseSpaceIntegral = phaseSpaceIntegral;

	// get wave list from production amplitudes and fill map for
	// production-amplitude indices to indices in normalization integral
	for (unsigned int i = 0; i < _prodAmpNames.size(); ++i) {
		const TString waveName = waveNameForProdAmp(i);
		if (find(_waveNames.begin(), _waveNames.end(), waveName.Data()) == _waveNames.end())
			_waveNames.push_back(waveName.Data());
		// look for index of first occurence
		unsigned int j;
		for (j = 0; j < _prodAmpNames.size(); ++j)
			if (prodAmpName(j).Contains(waveName))
				break;
		_normIntIndexMap[i] = j;
	}

	// check consistency
	if (_prodAmps.size() != _prodAmpNames.size())
		cout << "fitResult::fill(): warning: number of production amplitudes "
		     << "(" << _prodAmps.size() << ") does not match number of "
		     << "production amplitude names (" << _prodAmpNames.size() << ")." << endl;
	if (_prodAmps.size() != _fitParCovMatrixIndices.size())
		cout << "fitResult::fill(): warning: number of production amplitudes "
		     << "(" << _prodAmps.size() << ") does not match number of "
		     << "covariance matrix indices (" << _fitParCovMatrixIndices.size() << ")." << endl;
	if (   ((int)_waveNames.size() != _normIntegral.nrows())
	       || ((int)_waveNames.size() != _normIntegral.ncols()))
		cout << "fitResult::fill(): warning: number of waves (" << _waveNames.size()
		     << ") does not match size of normalization integral "
		     << "(" << _normIntegral.nrows() << ", " << _normIntegral.ncols() << ")." << endl;

	if (0) {
		// print debug information that allows to check ordering of arrays
		map<TString, complex<double> > V;
		for (unsigned int i = 0; i < nmbProdAmps(); ++i)
			V[prodAmpName(i)] = prodAmp(i);
		printInfo << "production amplitudes:" << endl;
		for (map<TString, complex<double> >::const_iterator i = V.begin(); i != V.end(); ++i)
			cout << "     " << i->first << " = " << i->second << endl;
		printInfo << "normalization integral matrix:" << endl;
		map<TString, unsigned int> waveIndexMap;
		for (unsigned int i = 0; i < nmbWaves(); ++i)
			waveIndexMap[waveName(i)] = i;
		for (map<TString, unsigned int>::const_iterator i = waveIndexMap.begin();
		     i != waveIndexMap.end(); ++i)
			for (map<TString, unsigned int>::const_iterator j = waveIndexMap.begin();
			     j != waveIndexMap.end(); ++j) {
				cout << "    [" << i->first << "][" << j->first << "] = "
				     << this->normIntegral(i->second, j->second) << endl;
			}
		printInfo << "covariance matrix:" << endl;
		map<TString, unsigned int> prodAmpIndexMap;
		for (unsigned int i = 0; i < nmbProdAmps(); ++i)
			prodAmpIndexMap[prodAmpName(i)] = i;
		for (map<TString, unsigned int>::const_iterator i = prodAmpIndexMap.begin();
		     i != prodAmpIndexMap.end(); ++i) {
			const int iIndex[2] = { _fitParCovMatrixIndices[i->second].first,    // real part index
			                        _fitParCovMatrixIndices[i->second].second};  // imaginary part index
			for (map<TString, unsigned int>::const_iterator j = prodAmpIndexMap.begin();
			     j != prodAmpIndexMap.end(); ++j) {
				const int jIndex[2] = { _fitParCovMatrixIndices[j->second].first,    // real part index
				                        _fitParCovMatrixIndices[j->second].second};  // imaginary part index
				cout << "    [" << i->first << "][" << j->first << "]: " << flush;
				if (iIndex[0] >= 0) {
					if (jIndex[0] >= 0)
						cout << "ReRe = " << fitParameterCov(iIndex[0], jIndex[0]) << ", " << flush;
					if (jIndex[1] >= 0)
						cout << "ReIm = " << fitParameterCov(iIndex[0], jIndex[1]) << ", " << flush;
				}
				if (iIndex[1] >= 0) {
					if (jIndex[0] >= 0)
						cout << "ImRe = " << fitParameterCov(iIndex[1], jIndex[0]) << ", " << flush;
					if (jIndex[1] >= 0)
						cout << "ImIm = " << fitParameterCov(iIndex[1], jIndex[1]) << flush;
				}
				cout << endl;
			}
		}
	}
}


int
fitResult::waveIndex(const string& waveName) const
{
	int index = -1;
	for (unsigned int i = 0; i < nmbWaves(); ++i)
		if (waveName == _waveNames[i]) {
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
