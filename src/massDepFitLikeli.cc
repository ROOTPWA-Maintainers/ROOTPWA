#include "massDepFitLikeli.h"
#include "TTree.h"
#include "TF1.h"
#include "TMatrixD.h"
#include "fitResult.h"
#include "TStopwatch.h"
#include <vector>
#include <complex>
#include <iostream>




using namespace std;

unsigned int
rpwa::massDepFitLikeli::NDim() const {return _compset->numPar();}

unsigned int
rpwa::massDepFitLikeli::NDataPoints() const {
  // loop through input data an check how many bins are in massrange
  unsigned int nbins=_tree->GetEntries();
  cout << nbins << " mass bins in input data." << endl;
  unsigned int nbinsInRange=0;
  // loop over mass-bins
  for(unsigned im=0;im<nbins;++im){
    _tree->GetEntry(im);
    double mass=_rhom->massBinCenter();
    if(mass>=_mmin && mass<=_mmax)++nbinsInRange;
  }
  cout << nbinsInRange << " mass bins in allowed mass range {"
       << _mmin<<","<<_mmax<<"}"<< endl;

  // calculate data points, remember (Re,Im)=> factor 2:
  // a complex, symmetric matrix with real diagonal has n^2 elements:
  return nbinsInRange*_wlist.size()*_wlist.size();
}

bool
rpwa::massDepFitLikeli::init(TTree* sp,
                             pwacompset* compset,
                             const vector<string>& waveNames,
                             const vector<pair<double, double> >& waveMassLimits,
                             bool doCov)
{
  _doCov=doCov;
  _compset=compset;
  _tree=sp;
  _rhom=NULL;//new fitResult();
  if(_tree->SetBranchAddress("fitResult_v2",&_rhom))cerr<<"Branch not found!"<<endl;

  _tree->GetEntry(0);

	// check that the wavelists in compset and here are identically
	// not that they should or could have changed, but better safe than sorry
	_wlist = waveNames;

	bool brokenWaveList = false;
	if(_wlist.size() != compset->getWaveList().size()) {
		            brokenWaveList = true;
	} else {
		for (size_t i=0; i<_wlist.size(); ++i) {
			if (_wlist[i] != compset->getWaveList()[i]) {
				brokenWaveList = true;
				break;
			}
		}
	}
	if(brokenWaveList) {
		printErr << "detected a messed up internal array." << endl;
		return false;
	}

  // build waveindex map;
  _index.clear();

  cerr << "Number of components: "<<_compset->n()<< endl;
  cerr << (*_compset) << endl;

  cerr << "Number of waves in fit: "<< _wlist.size() << endl;
// build map from waves as in _index to components and their channels
 for(unsigned int i=0;i<_wlist.size();++i){
   _index.push_back(_rhom->waveIndex(_wlist[i]));
   cerr <<  _wlist[i] << endl;
 }

 // read data into memory
 unsigned int nbins=_tree->GetEntries();
 unsigned int nwaves=_wlist.size();
 _spindens.clear();
 _cov.clear();
 _mass.clear();
 for(unsigned im=0;im<nbins;++im){
   _tree->GetEntry(im);
   _mass.push_back(_rhom->massBinCenter());

   ccmatrix mycov(nwaves,nwaves);
   cmatrix myrhom(nwaves,nwaves);
   // _vrhom.push_back(*_rhom);
   for(unsigned int i=0; i<nwaves; ++i){
     for(unsigned int j=i; j<nwaves; ++j){

       myrhom(i,j)=_rhom->spinDensityMatrixElem(_index[i],_index[j]);

       TMatrixD acov(2,2);
       acov= _rhom->spinDensityMatrixElemCov(_index[i],_index[j]);
       rmatrix c(2,2);c(0,0)=acov(0,0);c(1,0)=acov(1,0);c(1,1)=acov(1,1);c(0,1)=acov(0,1);
       mycov(i,j)=c;

     }
   }
   _spindens.push_back(myrhom);
   _cov.push_back(mycov);
 } // end loop over mass bins

	// sort mass bins
	sort(_mass.begin(), _mass.end());

	printInfo << "found " << _mass.size() << " mass bins, center of first and last mass bins: " << _mass[0] << " and " << _mass[_mass.size() - 1] << " MeV/c^2." << endl;

	const double massStep = (_mass[_mass.size() - 1] - _mass[0]) / (_mass.size() - 1);
	for(size_t idxMass=1; idxMass<_mass.size(); ++idxMass) {
		if(abs(_mass[idxMass]-_mass[idxMass-1] - massStep) > 1000.*numeric_limits<double>::epsilon()) {
			printErr << "mass distance between bins " << idxMass-1 << " (" << _mass[idxMass-1] << " MeV/c^2) and "
			         << idxMass << " (" << _mass[idxMass] << " MeV/c^2) does not agree with nominal distance "
			         << massStep << " MeV/c^2" << endl;
			return false;
		}
	}
	printInfo << "distance between two mass bins is " << massStep << " MeV/c^2." << endl;

	_mmin=_mass[0] - massStep / 2;
	_mmax=_mass[_mass.size() - 1] + massStep / 2;

	printInfo << "mass bins cover the mass range from " << _mmin << " to " << _mmax << " MeV/c^2." << std::endl;

	// determine bins to be used in the fit
	if(waveNames.size() != waveMassLimits.size()) {
		printErr << "detected a messed up internal array." << endl;
		return false;
	}

	std::vector<std::pair<size_t, size_t> > binLimits;
	for(size_t idxWave=0; idxWave<waveNames.size(); ++idxWave) {
		size_t binFirst = 0;
		size_t binLast = _mass.size()-1;
		for(size_t idxMass=0; idxMass<_mass.size(); ++idxMass) {
			if(_mass[idxMass] < waveMassLimits[idxWave].first) {
				binFirst = idxMass+1;
			}
			if(_mass[idxMass] == waveMassLimits[idxWave].first) {
				binFirst = idxMass;
			}
			if(_mass[idxMass] <= waveMassLimits[idxWave].second) {
				binLast = idxMass;
			}
		}
		if(waveMassLimits[idxWave].first < 0) {
			binFirst = 0;
		}
		if(waveMassLimits[idxWave].second < 0) {
			binLast = _mass.size()-1;
		}
		binLimits.push_back(make_pair(binFirst, binLast));
	}

	_massBinLimits.clear();
	_massBinLimits.resize(binLimits.size());
	for(size_t idxWave=0; idxWave<binLimits.size(); ++idxWave) {
		_massBinLimits[idxWave].resize(binLimits.size());
		for(size_t jdxWave=0; jdxWave<binLimits.size(); ++jdxWave) {
			_massBinLimits[idxWave][jdxWave] = make_pair(max(binLimits[idxWave].first,  binLimits[jdxWave].first),
			                                             min(binLimits[idxWave].second, binLimits[jdxWave].second));
		}
	}

	printInfo << "waves and mass limits:" << endl;
	if(_massBinLimits.size() != _wlist.size()) {
		printErr << "detected a messed up internal array." << endl;
		return false;
	}
	for(size_t idxWave=0; idxWave<binLimits.size(); ++idxWave) {
		if(_massBinLimits[idxWave].size() != _wlist.size()) {
			printErr << "detected a messed up internal array." << endl;
			return false;
		}
		ostringstream output;
		for(size_t jdxWave=0; jdxWave<binLimits.size(); ++jdxWave) {
			if(_massBinLimits[idxWave][jdxWave] != _massBinLimits[jdxWave][idxWave]) {
				printErr << "detected a messed up internal array." << endl;
				return false;
			}
			output << _massBinLimits[idxWave][jdxWave].first << "-" << _massBinLimits[idxWave][jdxWave].second << " ";
		}
		printInfo << _wlist[idxWave] << " " << binLimits[idxWave].first << "-" << binLimits[idxWave].second
		          << " (" << (waveMassLimits[idxWave].first < 0 ? _mmin : waveMassLimits[idxWave].first)
		          << "-" << (waveMassLimits[idxWave].second < 0 ? _mmax : waveMassLimits[idxWave].second)
		          << "): " << output.str() << endl;
	}

	return true;
}


double
rpwa::massDepFitLikeli::DoEval(const double* par) const {
  // rank==1 mass dependent fit

  // input values: m
  unsigned int nwaves=_index.size();

  // set parameters for resonances, background and phase space
  _compset->setPar(par);

  // Breitwigner parameters
  //masses;
  //widths;

  double chi2=0;
  unsigned int nbins=_mass.size();

   TStopwatch dataW;dataW.Stop();dataW.Reset();
  TStopwatch modelW;modelW.Stop();modelW.Reset();
  TStopwatch calcW;calcW.Stop();calcW.Reset();

  // loop over mass-bins
  for(unsigned im=0;im<nbins;++im){
    //_tree->GetEntry(im);
    const ccmatrix& mycov =_cov[im];
    const cmatrix& myrhom =_spindens[im];

    //const fitResult* fitR=&_vrhom[im];
    double mass=_mass[im];
    if(mass<_mmin)continue;
    if(mass>_mmax)continue;
    //if(mass>2000)continue;
    //cout << "Mass=" << mass << endl;
    // inpu values: measured spin density matrix
       
    // sum over the contributions to chi2 -> rho_ij
      for(unsigned int i=0; i<nwaves; ++i){
	for(unsigned int j=i; j<nwaves; ++j){
			// check that this mass bin should be taken into account for this
			// combination of waves
			if(im < _massBinLimits[i][j].first || im > _massBinLimits[i][j].second) {
				continue;
			}

	  // calculate target spin density matrix element
	  complex<double> rho;
	  // loop over all waves
	  complex<double> f1;
	  complex<double> f2;
	  const string w1=_wlist[i];
	  const string w2=_wlist[j];
	  // loop over components and look which contains these waves
	  modelW.Start(0);
	  rho=_compset->overlap(i,j,mass);
	  //complex<double> rho2=_compset->overlap(w1,w2,mass);
	  //if(norm(rho-rho2)>1E-4)cerr<< " ##### " << norm(rho-rho2) << endl;
	  modelW.Stop();
	  // compare to measured spin density matrix element
	  dataW.Start(0);
	  complex<double> rhom=rho-myrhom(i,j);//  fitR->spinDensityMatrixElem(_index[i],_index[j]);

	  const rmatrix& cov= mycov(i,j);//fitR->spinDensityMatrixElemCov(_index[i],_index[j]);
	  dataW.Stop();

	  calcW.Start(0);
	   double dchi=norm(rhom);
	  //cov.Print();
	   if(i==j)dchi/=cov(0,0);
	  else {

	    if(_doCov){
	      dchi=rhom.real()*rhom.real()*cov(1,1)-(2*cov(0,1))*rhom.real()*rhom.imag()+cov(0,0)*rhom.imag()*rhom.imag();
	      dchi/=cov(0,0)*cov(1,1)-cov(0,1)*cov(0,1);
	    }
	    else dchi=rhom.real()*rhom.real()/cov(0,0)+rhom.imag()*rhom.imag()/cov(1,1);
	  }

	  //cerr << "d-rho("<<i<<","<<j<<")=" << dchi<<endl;;
	  //cerr << "sigma ="<< sigma << endl;
	  //if(i==j)
	  chi2+=dchi;
	  calcW.Stop();
	}
      }// end loop over i
  } // end loop over mass-bins

  //cerr << "dataAccessTime: " << maxPrecisionAlign(dataW.CpuTime()) << " sec" << endl;
  //cerr << "modelAccessTime: " << maxPrecisionAlign(modelW.CpuTime()) << " sec" << endl;
  //cerr << "calcTime: " << maxPrecisionAlign(calcW.CpuTime()) << " sec" << endl;


  return chi2;
} // end DoEval
