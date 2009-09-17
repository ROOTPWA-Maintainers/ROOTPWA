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
// File and Version Information:
// $Id$
//
// Description:
//      Implementation of class TPWALikelihood
//      see TPWALikelihood.hh for details
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


// C/C++ Headers ----------------------
#include <iomanip>
#include <fstream>
#include <cassert>

// Collaborating Class Headers --------
#include "TString.h"
#include "TSystem.h"
#include "TCMatrix.h"
#include "TStopwatch.h"

#include "utilities.h"

// This Class' Header ------------------
#include "TPWALikelihood.h"


#define USE_FDF


using namespace std;


TPWALikelihood::TPWALikelihood()
  : _rank(1),
    _dim(0),
    _nmbEvents(0),
    _nmbWaves(0),
    _nmbWavesPosRefl(0),
    _nmbWavesNegRefl(0),
    _nmbEventsGrad(1000000),
    _Ltime(0),
    _Ntime(0),
    _debug(false),
    _useNorm(true)
{
  for (unsigned int i = 0; i < 4; ++i)
    _nmbCalls[i] = 0;
}


TPWALikelihood::~TPWALikelihood()
{
  clearCache();
}


// likelihood and first derivatives in one go
void
TPWALikelihood::FdF(const double* par,             // parameter array; reduced by rank conditions
		    double&       funcVal,         // function value
		    double*       gradient) const  // array of derivatives
{
  ++(_nmbCalls[FDF]);

  // log consumed time
  TStopwatch timer;
  timer.Start(true);

  // save argument in parameter cache
  for (unsigned int i = 0; i < _dim; ++i)
    _parCache[i] = par[i];

  // build complex production amplitudes from function parameters
  // taking into account rank restrictions
  double                            prodAmpFlat;
  vector<vector<complex<double> > > prodAmps     = copyFromParArray(par, prodAmpFlat);
  const double                      prodAmpFlat2 = prodAmpFlat * prodAmpFlat;

  double L = 0;  // log likelihood
  // derivatives
  // contains derivative with resp to real/imaginary part resp.
  // the dL are NOT well defined complex numbers!
  vector<vector<complex<double> > > dL(_rank, vector<complex<double> >(_nmbWaves, 0));
  vector<vector<complex<double> > > dl(_rank, vector<complex<double> >(_nmbWaves, 0));  // derivative contribution for this event
  double dLdflat = 0;
  // amplitude:
  // two contributions for two reflectivities
  complex<double> ampPos(0, 0);
  complex<double> ampNeg(0, 0); 
  // loop over events and calculate log likelihood and derivatives with respect to parameters
  for (unsigned int iEvt = 0; iEvt < _nmbEvents; ++iEvt) {
    double l = 0;  // likelihood contribution of this event
    // incoherent sum over ranks
    for (unsigned int r = 0; r < _rank; ++r) {
      for (unsigned int iWave = 0; iWave <_nmbWaves; ++iWave) {  // loop over waves
	// contribution to likelihood
	const complex<double> amp = prodAmps[r][iWave] * _decayAmps[iWave][iEvt];
	const int refl = _waveRefl[iWave];
	if (refl == -1)
	  ampNeg += amp;
	else
	  ampPos += amp;
	// contributions to derivatives
	if (iEvt < _nmbEventsGrad)
	  for (unsigned int jWave = 0; jWave < _nmbWaves; ++jWave)  // loop over derivatives 
	    // loop only inside current rank and only for waves with same reflectivity
	    if (refl == _waveRefl[jWave])
	      dl[r][jWave] += amp;
      } // end loop over waves
      l += norm(ampPos);
      l += norm(ampNeg);
      ampPos.real() = 0;
      ampPos.imag() = 0;
      ampNeg.real() = 0;
      ampNeg.imag() = 0;
      assert(l >= 0);
      if (iEvt < _nmbEventsGrad)
	for (unsigned int iWave = 0; iWave <_nmbWaves; ++iWave)  // again loop over derivatives inside current rank
	  dl[r][iWave] *= conj(_decayAmps[iWave][iEvt]);
    }  // end loop over rank
    l += prodAmpFlat2;
    // accumulate log likelihood
    L -= TMath::Log(l);
    // accumulate derivative
    if (iEvt < _nmbEventsGrad) {  //loop over derivatives to incorporate factor 2 / sigma
      const double g = 2. / l;
      for (unsigned int r = 0; r < _rank; ++r)
	for (unsigned int iWave = 0; iWave <_nmbWaves; ++iWave) {
	  dL[r][iWave] -= g * dl[r][iWave];
	  dl[r][iWave].real() = 0;
	  dl[r][iWave].imag() = 0;
	}
      dLdflat -= g * prodAmpFlat;
    }
  }  // end loop over events for calculation of intensity term
 
  // rescale derivatives
  const double scale = (double)_nmbEvents / (double)_nmbEventsGrad;
  for (unsigned int r = 0; r < _rank; ++r)
    for (unsigned int iWave = 0; iWave <_nmbWaves; ++iWave)
      dL[r][iWave] *= scale;
  dLdflat *= scale;
  
  // log consumed time
  const double t1 = timer.RealTime();
  timer.Start(true);

  // add normalization integrals ?? TODO: can we exploit symmetry here?
  complex<double> N(0, 0);
  // event number normalization
  const double nmbEvt      = (_useNorm) ? 1 : (double)_nmbEvents;
  const double twiceNmbEvt = 2 * nmbEvt;
  // normalization contribution of derivative:
  vector<vector<complex<double> > > dLN(_rank, vector<complex<double> >(_nmbWaves, 0));
  for (unsigned int r = 0; r < _rank; ++r) {  // loop over rank
    for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave) {  // outer loop
      for (unsigned int jWave = 0; jWave < _nmbWaves; ++jWave) {  // inner loop
	// check if i and j are of the same reflectivity
	if (_waveRefl[iWave] != _waveRefl[jWave])
	  continue;
	complex<double> p = prodAmps[r][iWave] * conj(prodAmps[r][jWave]);
	//cout << _waveNames[iWave] << "   " << _waveNames[jWave] << endl;
	const complex<double> I = _accMatrix.element(iWave, jWave);
	//cout << "I = " << I << endl;
	//cout << "p = " << p << endl;
	p *= I;
	N += p;
	// calculate contribution to derivative
	dLN[r][iWave] += prodAmps[r][jWave] * conj(I);
      }  // end inner loop over waves
      // account for 2nmbEvents
      dL[r][iWave] += dLN[r][iWave] * twiceNmbEvt;
      //cout << "df(r = " << r << ", i = " << iWave << ") = " << dL[r][iWave] << endl;
    }  // end outer loop over waves
  }  // end loop over rank
  // take care of derivative for fla
  dLdflat  += twiceNmbEvt * prodAmpFlat;
  N.real() += prodAmpFlat2;

  // sort derivative results into output array and cache
  copyToParArray(dL, dLdflat, gradient);
  copyToParArray(dL, dLdflat, &(*(_derivCache.begin())));

  // log consumed time
  const double t2 = timer.RealTime();
  timer.Stop();
  _Ltime += t1;
  _Ntime += t2;
  
  if (_debug)
    printInfo << "Log likelihood =  " << maxPrecisionAlign(L) << ", "
	      << "normalization =  " << maxPrecisionAlign(N.real()) << ", "
	      << "normalized likelihood = " << maxPrecisionAlign(L + nmbEvt * N.real()) << endl
	      << "    Time for likelihood = " << t1 << ", time for normalization = " << t2 << endl
	      << "    Number of events used to calculate gradient = " << _nmbEventsGrad << endl;

  // return likelihood value
  funcVal = L + nmbEvt * N.real();
}


// calculate derivatives with respect to parameters
void
TPWALikelihood::Gradient(const double* par,             // parameter array; reduced by rank conditions
			 double*       gradient) const  // array of derivatives
{
  ++(_nmbCalls[GRADIENT]);
  
#ifdef USE_FDF

  // check whether parameter is in cache
  bool samePar = true;
  for (unsigned int i = 0; i < _dim; ++i)
    if (_parCache[i] != par[i]) {
      samePar = false;
      break;
    }
  if (samePar) {
    //cout << " using cached values." << endl;
    for (unsigned int i = 0; i < _dim ; ++i)
      gradient[i] = _derivCache[i];
    return;
  }
  //cout << endl;
  // call FdF
  double logLikelihood;
  FdF(par, logLikelihood, gradient);

#else

  // log consumed time
  TStopwatch timer;
  timer.Start(true);

  // copy arguments into parameter cache
  for (unsigned int i = 0; i < _dim; ++i)
    _parCache[i] = par[i];

  // build complex production amplitudes from function parameters
  // taking into account rank restrictions
  double                            prodAmpFlat;
  vector<vector<complex<double> > > prodAmps = copyFromParArray(par, prodAmpFlat);

  // array with likelihood derivative with resp to real and imaginary
  // parts of the production amplitudes although stored as complex
  // values, the dL are _not_ well defined complex numbers!
  double                            derivativeFlat = 0;
  vector<vector<complex<double> > > derivatives(_rank, vector<complex<double> >(_nmbWaves, 0));

  // compute derivative for first term of log likelihood
  for (unsigned int iEvt = 0; iEvt < _nmbEvents; ++iEvt) {
    double                            l = 0;                                             // likelihood for this event
    vector<vector<complex<double> > > d(_rank, vector<complex<double> >(_nmbWaves, 0));  // likelihood derivative for this event
    for (unsigned int r = 0; r < _rank; ++r) {  // incoherent sum over ranks
      complex<double> ampPos(0, 0);  // positive reflectivity amplitude for this rank
      complex<double> ampNeg(0, 0);  // negative reflectivity amplitude for this rank
      for (unsigned int iWave = 0; iWave <_nmbWaves; ++iWave) {  // outer loop over waves
	const complex<double> amp  = prodAmps[r][iWave] * _decayAmps[iWave][iEvt];
	const int             refl = _waveRefl[iWave];
	if (refl == -1)
	  ampNeg += amp;
	else
	  ampPos += amp;
	if (iEvt < _nmbEventsGrad)
	  for (unsigned int jWave = 0; jWave < _nmbWaves; ++jWave)  // inner loop over waves
	    // sum up amplitudes for current rank and only for waves with same reflectivity
	    if (refl == _waveRefl[jWave])
	      d[r][jWave] += amp;
      }  // end outer loop over waves
      l += norm(ampPos) + norm(ampNeg);
      assert(l >= 0);
      // loop again over waves for current rank
      if (iEvt < _nmbEventsGrad)
	for (unsigned int iWave = 0; iWave <_nmbWaves; ++iWave)
	  d[r][iWave] *= conj(_decayAmps[iWave][iEvt]);
    }  // end loop over rank
    l += prodAmpFlat * prodAmpFlat;
    // incorporate factor 2 / sigma
    if (iEvt < _nmbEventsGrad) {
      const double factor = 2. / l;
      for (unsigned int r = 0; r < _rank; ++r)
	for (unsigned int iWave = 0; iWave <_nmbWaves; ++iWave)
	  derivatives[r][iWave] -= factor * d[r][iWave];
      derivativeFlat -= factor * prodAmpFlat;
    }
  }  // end loop over events
 
  // rescale derivatives
  if (_nmbEventsGrad != _nmbEvents) {
    const double scale = (double)_nmbEvents / (double)_nmbEventsGrad;
    for (unsigned int r = 0; r < _rank; ++r)
      for (unsigned int iWave = 0; iWave <_nmbWaves; ++iWave)
	derivatives[r][iWave] *= scale;
    derivativeFlat *= scale;
  }
  
  // log consumed time
  const double t1 = timer.RealTime();
  timer.Start(true);

  // normalize derivatives
  const double nmbEvt      = (_useNorm) ? 1 : (double)_nmbEvents;
  const double twiceNmbEvt = 2 * nmbEvt;
  for (unsigned int r = 0; r < _rank; ++r) {  // loop over rank
    for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave) {  // outer loop
      const int refl = _waveRefl[iWave];
      complex<double> normFactor(0, 0);
      for (unsigned int jWave = 0; jWave < _nmbWaves; ++jWave) {  // inner loop
	if (refl != _waveRefl[jWave])  // make sure that waves i and j have the same reflectivity
	  continue;
	normFactor += prodAmps[r][jWave] * conj(_accMatrix.element(iWave, jWave));
      }  // end inner loop over waves
      derivatives[r][iWave] += twiceNmbEvt * normFactor;  // account for 2 * # of events
    }  // end outer loop over waves
  }  // end loop over rank
  derivativeFlat += twiceNmbEvt * prodAmpFlat;

  // return gradient values
  copyToParArray(derivatives, derivativeFlat, gradient);
  copyToParArray(derivatives, derivativeFlat, &(*(_derivCache.begin())));

  // log consumed time
  const double t2 = timer.RealTime();
  timer.Stop();
  _Ltime += t1;
  _Ntime += t2;

  if (0) {  // compare to FdF
    double f;
    double df[_dim];
    FdF(par, f, df);
    double maxDelta = 0;
    for (unsigned int i = 0; i < _dim; ++i) {
      const double delta = df[i] - gradient[i];
      if (delta > maxDelta)
	maxDelta = delta;
    }
    if (maxDelta > 1e-12)
      cout << "!!! max delta df = " << maxPrecision(maxDelta) << endl;
  }

#endif  // USE_FDF
}


unsigned int
TPWALikelihood::NDim() const
{ return _dim; }


// reads wave names and thresholds from wave list file
void 
TPWALikelihood::SetWavelist(const TString& waveListFileName)
{
  printInfo << "Reading amplitude names from wave list file '" << waveListFileName << "'." << endl;
  ifstream waveListFile(waveListFileName.Data());
  if (!waveListFile) {
    printErr << "Cannot open file '" << waveListFileName << "'. Exiting." << endl;
    throw;
  }
  string line;
  int    lineNmb = 0;
  while (getline(waveListFile, line)) {
    if (line[0] == '#')  // comments start with #
      continue;
    stringstream lineStream;
    lineStream.str(line);
    string waveName;
    if (lineStream >> waveName) {
      double threshold;
      // !!! it would be safer to make the threshold value in the wave list file mandatory
      if (!(lineStream >> threshold))
	threshold = 0;
      if (_debug)
	cout << "    reading line " << setw(3) << lineNmb + 1 << ": " << waveName << ",  threshold = " << setw(4) << threshold << " MeV/c^2" << endl;
      if (getReflectivity(waveName) > 0)
	++_nmbWavesPosRefl;
      else
	++_nmbWavesNegRefl;
      _waveNames.push_back(waveName);
      _waveThresholds.push_back(threshold);
    } else
      printWarn << "Cannot parse line '" << line << "' in wave list file '" << waveListFileName << "'." << endl;
    ++lineNmb;
  }
  waveListFile.close();
  printInfo << "Read " << lineNmb << " lines from wave list file '" << waveListFileName << "'." << endl;
  // calculate dimension
  _nmbWaves=_waveNames.size();
  //SetRank(_rank);
}


// prepares data structures and constructs parameter names
void
TPWALikelihood::SetRank(const unsigned int rank)
{
  _rank = rank;
  // calculate dimension
  _dim = 0;
  for (unsigned int r = 0; r < _rank; ++r) {
    int pref  = _nmbWavesPosRefl - r;  // number of positive reflectivity waves in this rank
    int prefc = pref - 1;              // complex valued
    if (prefc < 0)
      prefc = 0;
    int nref  = _nmbWavesNegRefl - r;  // number of negative reflectivity waves in this rank
    int nrefc = nref - 1;              // complex valued
    if (nrefc < 0)
      nrefc = 0;
    _dim += 2 * (prefc + nrefc);  // 2 parameters per complex amplitude
    if (pref > 0)
      ++_dim;  // real amp
    if (nref > 0)
      ++_dim;  // real amp
  }
  _dim += 1;  // additonal flat wave
  // this takes into account rank and the flat wave
  printInfo << "Dimension of likelihood function is " << _dim << "." << endl;
  _parCache.resize(_dim, 0);
  _derivCache.resize(_dim, 0);
  _parNames.resize(_dim, "");
  _waveRefl.resize(_nmbWaves + 1, 0);
  _parThresholds.resize(_dim, 0);
  // build parameter names
  unsigned int outCount = 0;  // output array counter
  for (unsigned int r = 0; r < _rank; ++r) {
    unsigned int posCount = 0; // counter for positive reflectivity waves
    unsigned int negCount = 0; // counter for negative reflectivity waves
    for (unsigned int i = 0; i < _nmbWaves; ++i) {
      int refl = getReflectivity(_waveNames[i]);
      _waveRefl[i] = refl;
      if ((refl > 0) && (posCount < r)) {  // production amplitude is zero
	++posCount;
	continue;
      } else if ((refl < 0) && (negCount < r)){  // production amplitude is zero
	++negCount;
	continue;
      } else if (((refl > 0) && (posCount == r)) || ((refl < 0) && (negCount == r))) {  // production amplitude is real
	ostringstream parName;
	parName << "V" << r << "_" << _waveNames[i] << "_RE";
	_parNames[outCount]      = parName.str();
	_parThresholds[outCount] = _waveThresholds[i];
	++outCount;
      } else {  // production amplitude is complex
	ostringstream parName;
	parName << "V" << r << "_" << _waveNames[i];
	_parNames[outCount]      = parName.str() + "_RE";
	_parThresholds[outCount] = _waveThresholds[i];
	++outCount;
	_parNames[outCount]      = parName.str() + "_IM";
	_parThresholds[outCount] = _waveThresholds[i];
	++outCount;
      }
      if (refl > 0)
	++posCount;
      else
	++negCount;
    }  // end loop over waves
  }  // end loop over rank
  _parNames[outCount]      = "V_flat";
  _parThresholds[outCount] = 0;
  _waveRefl[_nmbWaves]     = 1;
}


void
TPWALikelihood::SetMaxSampDL(const unsigned int nmb)
{
  if ((_nmbEvents > 0) && (_nmbEvents < nmb))
    _nmbEventsGrad = _nmbEvents;
  else
    _nmbEventsGrad = nmb;
}


void
TPWALikelihood::LoadIntegrals(const TString& normIntFileName,
			      const TString& accIntFileName)
{
  printInfo << "Loading normalization integral from '" << normIntFileName << "'." << endl;
  ifstream intFile(normIntFileName.Data());
  if (!intFile) {
    printErr << "Cannot open file '" << normIntFileName << "'. Exiting." << endl;
    throw;
  }
  // !!! integral.scan() performs no error checks!
  _normInt.scan(intFile);
//   if (_debug)
//     _normInt.print(cout);
  intFile.close();
  _normMatrix = reorderedIntegralMatrix(_normInt);

  printInfo << "Loading acceptance integral from '" << accIntFileName << "'." << endl;
  intFile.open(accIntFileName.Data());
  if (!intFile) {
    printErr << "Cannot open file '" << accIntFileName << "'. Exiting." << endl;
    throw;
  }
  // !!! integral.scan() performs no error checks!
  _accInt.scan(intFile);
//   if (_debug)
//     _accInt.print(cout);
  intFile.close();
  //_accInt.events(100000); TODO: add possibility to rescale here!
  _accMatrix = reorderedIntegralMatrix(_accInt);
}


// reads in all amplitude files and stores values in memory
void
TPWALikelihood::LoadAmplitudes(){
  printInfo << "Loading amplitude data." << endl;
  
  // normalization integrals need to be loaded
  if (_normMatrix.nrows() == 0) {
    printErr << "Normalization integrals have to be loaded before loading the amplitudes. Aborting." << endl;
    throw;
  }

  // monitor memory usage
  ProcInfo_t infoBefore;
  if (_debug) {
    gSystem->GetProcInfo(&infoBefore);
    printInfo << "Resident memory usage before loading amplitudes: " << infoBefore.fMemResident << "." << endl;
  }
  
  // clear cache
  clearCache();

  // loop over amplitudes and read in data
  for (unsigned int i = 0; i < _nmbWaves; ++i) {
    if (_debug)
      cout << "    Loading amplitude data for wave '" << _waveNames[i] << "'." << endl;
    ifstream ampFile(_waveNames[i].c_str());  // this assumes that program is executed in amplitude directory
    if (!ampFile) {
      printErr << "Cannot open file '" << _waveNames[i] << "'. Exiting." << endl;
      throw;
    }
    // get integral
    complex<double> normInt = _normMatrix.element(i, i);
    // create data cache and fill it from file
    _decayAmps.push_back(vector<complex<double> >());
    _decayAmps.back().reserve(_nmbEvents);  // _nmbEvents 0 for the first iteration
    complex<double> amp;
    while (ampFile.read((char*) &amp, sizeof(complex<double>))) {
      if (_useNorm)  // normalize data, if option is switched on
	amp /= sqrt(normInt.real());  // rescale decay amplitude
      _decayAmps.back().push_back(amp);
    }
    if (_debug)
      cout << "    Read " << _decayAmps.back().size() << " events from file '" << _waveNames[i] << "'." << endl;
    _nmbEvents = _decayAmps.back().size(); 
  }
  bool allWavesHaveSameNmbEvent = true;
  for (unsigned int i = 1; i < _decayAmps.size(); ++i)
    if (_decayAmps[i - 1].size() != _decayAmps[i].size()) {
      allWavesHaveSameNmbEvent = false;
      break;
    }
  if ((!allWavesHaveSameNmbEvent) || (_decayAmps.back().size() != _nmbEvents))
    printWarn << "Amplitude files do not contain the same number of events for all waves." << endl;
  printInfo << "Loaded " << _nmbEvents << " events into memory." << endl;

  if (_nmbEvents < _nmbEventsGrad)
    _nmbEventsGrad = _nmbEvents;

  // monitor memory usage
  if (_debug) {
    ProcInfo_t infoAfter;
    gSystem->GetProcInfo(&infoAfter);
    printInfo << "Resident memory usage after loading amplitudes: " << infoAfter.fMemResident << ". "
	      << "Memory used for amplitudes: " << infoAfter.fMemResident - infoBefore.fMemResident << "." << endl;
  }

  // rescale integrals, if necessary
  if (_useNorm) {
    // remember: we have two copies of the matrix! -> change both!
    // so far only one is changed!!!!!!!!!!!!!!!
    printInfo << "Rescaling integrals." << endl;
    // normalization integral
    for (unsigned int i = 0; i < _nmbWaves; ++i) {  // outer loop over waves
      double Ni = sqrt(_normMatrix.element(i, i).real());
      for (unsigned int j = 0; j < _nmbWaves; ++j) {  // inner loop over waves
	if (i == j)
	  continue;  // do diagonal terms later
	double Nj = sqrt(_normMatrix.element(j, j).real());
	_normMatrix.element(i, j) /= Ni * Nj;
      }
    }
    // acceptance integral
    // Remember that the matrices _normMatrix (??) and _accMatrix are already normalized
    // to the number of Monte Carlo events!
    for (unsigned int i = 0; i < _nmbWaves; ++i) {   // outer loop over waves
      double Ni = sqrt(_normMatrix.element(i, i).real());
      for (unsigned int j = 0; j < _nmbWaves; ++j) {  // inner loop over waves
	double Nj = sqrt(_normMatrix.element(j, j).real());
	_accMatrix.element(i, j) /= Ni * Nj;
// !!! magic do-nothing line; without it fit result changes
// this is probably due to the -O3 used in compilation which might break some things
// has no influence when using -O0
	cout << "";
      }
    }
    // calculate diagonal terms of normalization matrices
    //!!! why?
    for (unsigned int i = 0; i < _nmbWaves; ++i)
      _normMatrix.element(i, i) = 1;
    if (_debug)
      for (unsigned int i = 0; i < _nmbWaves; ++i)
	for (unsigned int j = 0; j < _nmbWaves; ++j)
	  cout << "    normalization matrix [" << setw(3) << i << ", " << setw(3) << j << "] = "
 	       << "("  << maxPrecisionAlign(_normMatrix.element(i, j).real())
 	       << ", " << maxPrecisionAlign(_normMatrix.element(i, j).imag()) << "), "
	       << "acceptance matrix [" << setw(3) << i << ", " << setw(3) << j << "] = "
 	       << "("  << maxPrecisionAlign(_accMatrix.element(i, j).real())
 	       << ", " << maxPrecisionAlign(_accMatrix.element(i, j).imag()) << ")"
	       << endl;
  }  // end if useNorm
}


void
TPWALikelihood::getIntCMatrix(TCMatrix& normMatrix,
			      TCMatrix& accMatrix)
{
  //normMatrix.ResizeTo(_nmbWaves,_nmbWaves);
  for (unsigned int i = 0; i < _nmbWaves; ++i)  // outer loop
    for (unsigned int j = 0; j < _nmbWaves; ++j) {  // inner loop
      normMatrix.set(i, j, _normMatrix.element(i, j));
      accMatrix.set (i, j, _accMatrix.element (i, j));
    }
  // add flat
  normMatrix.set(_nmbWaves, _nmbWaves, 1);
  accMatrix.set (_nmbWaves, _nmbWaves, 1);
}


// Complex valued Amplitudes and
// mapping of real and imaginary part of amplitudes
// in error matrix (needed by TFitBin)
// takes into account real-valued parameters
void
TPWALikelihood::buildCAmps(const double*             x,
			   vector<complex<double> >& V,
			   vector<pair<int,int> >&   indices,
			   vector<string>&           names,
			   const bool                withFlat)
{
  // build complex numbers from parameters
  // remember rank restrictions!
  V.clear();
  indices.clear();
  names.clear();
  //unsigned int namp=_rank*_nmbWaves;
  //if (withFlat)namp+=1;
  unsigned int outCount = 0;  // output array counter
  for (unsigned int r = 0; r < _rank; ++r) {
    unsigned int posCount = 0; // counter for positive reflectivity waves
    unsigned int negCount = 0; // counter for negative reflectivity waves
    for (unsigned int i = 0; i < _nmbWaves; ++i) {
      int refl = getReflectivity(_waveNames[i]);
      double re, im;
      if (((refl > 0) && (posCount == r)) || ((refl < 0) && (negCount == r))) {  // real production amplitude
	indices.push_back(pair<int, int>(outCount, -1));
	re = x[outCount++];
	im = 0;
	V.push_back(complex<double>(re, im));
	stringstream name;
	name << "V" << r << "_" << _waveNames[i];
	names.push_back(name.str());
      } else if (((refl > 0) && (posCount > r)) || ((refl < 0) && (negCount > r))) {  // complex production amplitude
	indices.push_back(pair<int, int>(outCount, outCount + 1));
	re = x[outCount++];
	im = x[outCount++];
	V.push_back(complex<double>(re, im));
	stringstream name;
	name << "V" << r << "_" << _waveNames[i];
	names.push_back(name.str());
      }
      if (refl > 0)
	++posCount;
      else
	++negCount;
      //cout << "Wave" << i << "=" << V[r * _nmbWaves + i] << endl;
    }
  } // end loop over rank
  if (withFlat){
    V.push_back(complex<double>(x[outCount], 0));
    indices.push_back(pair<int, int>(outCount, -1));
    names.push_back("V_flat");
  }
}


double
TPWALikelihood::DoEval(const double* par) const
{
  ++(_nmbCalls[DOEVAL]);

#ifdef USE_FDF

  // call FdF
  double logLikelihood;
  double gradient[_dim];
  FdF(par, logLikelihood, gradient);
  return logLikelihood;

#else

  // log consumed time
  TStopwatch timer;
  timer.Start(true);

  // build complex production amplitudes from function parameters taking into account rank restrictions
  double                            prodAmpFlat;
  vector<vector<complex<double> > > prodAmps = copyFromParArray(par, prodAmpFlat);

  // compute first term of log likelihood: sum over data
  double logLikelihood = 0;
  for (unsigned int iEvt = 0; iEvt < _nmbEvents; ++iEvt) {
    double l = 0;  // likelihood for this event
    for (unsigned int r = 0; r < _rank; ++r) {  // incoherent sum over ranks
      complex<double> ampPos(0, 0);
      complex<double> ampNeg(0, 0); 
      for (unsigned int iWave = 0; iWave <_nmbWaves; ++iWave) {  // loop over waves
	const complex<double> amp  = prodAmps[r][iWave] * _decayAmps[iWave][iEvt];
	const int             refl = _waveRefl[iWave];
	if (refl == -1)
	  ampNeg += amp;
	else
	  ampPos += amp;
      }  // end loop over waves
      l += norm(ampPos) + norm(ampNeg);
      assert(l >= 0);
    }  // end loop over rank
    l += prodAmpFlat * prodAmpFlat;
    logLikelihood -= TMath::Log(l);  // accumulate log likelihood
  }  // end loop over events for calculation of intensity term
 
  // log consumed time
  const double t1 = timer.RealTime();
  timer.Start(true);

  // compute second term of log likelihood: normalization
  complex<double> normFactor(0, 0);
  for (unsigned int r = 0; r < _rank; ++r)  // loop over rank
    for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave)  // outer loop
      for (unsigned int jWave = 0; jWave < _nmbWaves; ++jWave) {  // inner loop
	if (_waveRefl[iWave] != _waveRefl[jWave])  // make sure that waves i and j have the same reflectivity
	  continue;
	normFactor += prodAmps[r][iWave] * conj(prodAmps[r][jWave]) * _accMatrix.element(iWave, jWave);
      }
  normFactor.real() += prodAmpFlat * prodAmpFlat;

  // log consumed time
  const double t2 = timer.RealTime();
  timer.Stop();
  _Ltime += t1;
  _Ntime += t2;
  
  const double nmbEvt = (_useNorm) ? 1 : (double)_nmbEvents;
  if (_debug)
    printInfo << "Log likelihood =  " << maxPrecisionAlign(logLikelihood) << ", "
	      << "normalization =  (" << maxPrecisionAlign(normFactor.real()) << ", " << maxPrecisionAlign(normFactor.imag()) << "), "
	      << "normalized log likelihood = " << maxPrecisionAlign(logLikelihood + nmbEvt * normFactor.real()) << endl
	      << "    Time for likelihood = " << t1 << ", time for normalization = " << t2 << endl;

  if (0) {  // compare to FdF
    double f;
    double df[_dim];
    FdF(par, f, df);
    const double funcVal = logLikelihood + nmbEvt * normFactor.real();
    const double delta   = f - funcVal;
    if (delta > 1e-12)
      cout << "!!! delta f = " << maxPrecision(delta) << endl;
  }

  // return log likelihood value
  return logLikelihood + nmbEvt * normFactor.real();

#endif  // USE_FDF
}


double
TPWALikelihood::DoDerivative(const double* par,
			     unsigned int  derivativeIndex) const
{
  ++(_nmbCalls[DODERIVATIVE]);

  // check whether parameter is in cache
  bool samePar = true;
  for (unsigned int i = 0; i < _dim; ++i)
    if (_parCache[i] != par[i]) {
      samePar = false;
      break;
    }
  if (samePar) {
    //cout << "using cached derivative! " << endl;
    return _derivCache[derivativeIndex];
  }
  // call FdF
  double logLikelihood;
  double gradient[_dim];
  FdF(par, logLikelihood, gradient);
  return gradient[derivativeIndex];
}


void
TPWALikelihood::clearCache()
{
  for (unsigned int i = 0; i < _decayAmps.size(); ++i)
    _decayAmps[i].clear();
  _decayAmps.clear();
}


// depends on naming convention for waves!!!
// VR_IGJPCMEIso....
int
TPWALikelihood::getReflectivity(const TString& waveName) const
{
  int refl = 0;
  unsigned int reflIndex = 6;  // position of reflectivity in wave
  // check whether it is parameter or wave name
  if (waveName[0] == 'V')
    reflIndex = 9; 
  if (waveName[reflIndex] == '-')
    refl= -1;
  else if (waveName[reflIndex] == '+')
    refl= +1;
  else {
    printErr << "Cannot parse parameter/wave name '" << waveName << "'. Cannot not determine reflectivity. Aborting." << endl;
    throw;
  }
  if (_debug)
    printInfo << "Extracted reflectivity = " << refl << " from parameter name '" << waveName << "' (char position " << reflIndex << ")" << endl;
  return refl;
}


// returns integral matrix reordered according to _waveNames array
matrix<complex<double> >
TPWALikelihood::reorderedIntegralMatrix(integral& integral) const
{
  // get original matrix and list of wave names
  const matrix<complex<double> > intMatrix    = integral.mat();
  const list<string>             intWaveNames = integral.files();
  // build index lookup-table
  vector<unsigned int> indexLookUp;  // lookup table: wave index -> index in normalization integral
  indexLookUp.resize(_nmbWaves, 0);
  for (unsigned int i = 0; i < _nmbWaves; ++i) {
    if (find(intWaveNames.begin(), intWaveNames.end(), _waveNames[i]) == intWaveNames.end()) {
      printErr << "Wave " << _waveNames[i] << " is not in integral. Aborting." << endl;
      throw;
    }
    indexLookUp[i] = integral.index(_waveNames[i]);  // not a const method!!! should be fixed
    if (_debug)
      cout << "    Mapping wave " << setw(3) << i << " '" << _waveNames[i] << "' "
	   << "to index " << setw(3) << indexLookUp[i] << " in integral." << endl;
  }
  // create reordered matrix
  matrix<complex<double> > reorderedMatrix(intMatrix.nrows(), intMatrix.ncols());
  for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave)
    for (unsigned int jWave = 0; jWave < _nmbWaves; ++jWave)
      reorderedMatrix.element(iWave, jWave) = intMatrix.element(indexLookUp[iWave], indexLookUp[jWave]);
  return reorderedMatrix;
}


// copy values from array that corresponds to the function parameters
// to structure that corresponds to the complex production amplitudes
// taking into account rank restrictions
vector<vector<complex<double> > >
TPWALikelihood::copyFromParArray(const double* inPar,             // input parameter array
				 double&       outFlatVal) const  // output value corresponding to flat wave
{
  vector<vector<complex<double> > > outVal(_rank, vector<complex<double> >(_nmbWaves, 0));
  unsigned int inCount = 0;  // input array counter
  for (unsigned int r = 0; r < _rank; ++r) {
    unsigned int posCount = 0;  // counter for positive reflectivity waves
    unsigned int negCount = 0;  // counter for negative reflectivity waves
    double re, im;
    for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave) {
      const int refl = _waveRefl[iWave];
      if (((refl > 0) && (posCount < r)) || ((refl < 0) && (negCount < r)))  // production amplitude is zero
	re = im = 0;
      else if (((refl > 0) && (posCount == r)) || ((refl < 0) && (negCount == r))) {  // production amplitude is real
	re = inPar[inCount++];
	im = 0;
      } else {  // production amplitude is complex
	re = inPar[inCount++];
	im = inPar[inCount++];
      }
      outVal[r][iWave] = complex<double>(re, im);
      if (refl > 0)
	++posCount;
      else
	++negCount;
    }  // end loop over waves
  }  // end loop over rank
  outFlatVal = inPar[inCount];
  return outVal;
}


// copy values from structure that corresponds to the complex
// production amplitudes to array that corresponds to the function
// parameters taking into account rank restrictions
void
TPWALikelihood::copyToParArray(const vector<vector<complex<double> > >& inVal,         // values corresponding to production amplitudes
			       const double                             inFlatVal,     // value corresponding to flat wave
			       double*                                  outPar) const  // output parameter array
{
  unsigned int outCount = 0;  // output array counter
  for (unsigned int r = 0; r < _rank; ++r) {
    unsigned int posCount = 0;  // counter for positive reflectivity waves
    unsigned int negCount = 0;  // counter for negative reflectivity waves
    for (unsigned int iWave = 0; iWave < _nmbWaves; ++iWave) {
      const int refl = _waveRefl[iWave];
      if (((refl > 0) && (posCount == r)) || ((refl < 0) && (negCount == r)))  // production amplitude is real
	outPar[outCount++] = inVal[r][iWave].real();
      else if (((refl > 0) && (posCount > r)) || ((refl < 0) && (negCount > r))) {  // production amplitude is complex
	outPar[outCount++] = inVal[r][iWave].real();
	outPar[outCount++] = inVal[r][iWave].imag();
      }
      if (refl > 0)
	++posCount;
      else
	++negCount;
    }  // end loop over waves
  }  // end loop over rank
  outPar[outCount] = inFlatVal;
}
