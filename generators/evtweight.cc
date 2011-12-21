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

/** @brief Calculate weight of event list from amp files
 */


#include <complex>
#include <iostream>
//#include <stringstream>

#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <unistd.h>
#include <stdlib.h>

#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

#include "TPWWeight.h"
#include "TDiffractivePhaseSpace.h"

#include "fitResult.h"
#include "event.h"

//#include "diffractiveDissVertex.h"
#include "../amplitude/particle.h"
#include "particleDataTable.h"
#include "libconfig.h++"

using namespace std;
using namespace libconfig;
using namespace rpwa;

extern ::particleDataTable PDGtable;

void printUsage(char* prog, int errCode=0) {
cerr << "usage:" << endl
     << prog
     << " -e <file> -o <file> -w <file> -i <file> -n samples  -m mass"
     << "    where:" << endl
     << "        -e <file> acc or ps events in .evt or .root format"<< endl
     << "        -o <file>  ROOT output file"<< endl
     << "        -w <file.root>  use TFitBin tree as input"<< endl 
     << "        -i <file>  integral file"<< endl 
     << "        -m mass  center of mass bin"<< endl   
     << "        -n samples  sampling of the model parameters (default=1)"<< endl   
     << "        -b width  width of mass bin"<< endl   
     << endl;
 exit(errCode);
}


int getReflectivity(const TString& waveName)
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
  return refl;
}



void parseWaveList(const string& waveListFileName, 
		   vector<string>& waveNames,
		   vector<double>& waveThresholds,
		   vector<int>& refl){
  unsigned int _nmbWavesPosRefl=0;
  unsigned int _nmbWavesNegRefl=0;
  bool _debug=false;
  printInfo << "Reading amplitude names from wave list file '" << waveListFileName << "'." << endl;
  ifstream waveListFile(waveListFileName.c_str());
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
      if (getReflectivity(waveName) > 0){
	++_nmbWavesPosRefl; refl.push_back(1);
      }
      else {
	++_nmbWavesNegRefl; refl.push_back(-1);
      }
      waveNames.push_back(waveName);
      waveThresholds.push_back(threshold);
    } else
      printWarn << "Cannot parse line '" << line << "' in wave list file '" << waveListFileName << "'." << endl;
    ++lineNmb;
  }
  waveListFile.close();
  printInfo << "Read " << lineNmb << " lines from wave list file '" << waveListFileName << "'." << endl;
}



int main(int argc, char** argv)
{

  if(argc<3)printUsage(argv[0],1);

  string output_file("genpw.root");
  string integrals_file;
  string waveListFileName("wavelist");
  string evtfilename;
  double binCenter = 0;
  double binWidth = 60; // MeV
  unsigned int nsamples = 1;

  int c;
  while ((c = getopt(argc, argv, "e:o:w:i:m:n:h")) != -1)
    switch (c) {
    case 'e':
      evtfilename = optarg;
      break;
    case 'o':
      output_file = optarg;
      break;
   case 'w':
      waveListFileName = optarg;
      break;
   case 'i':
      integrals_file = optarg;
      break;
   case 'h':
      printUsage(argv[0]);
      break;
   case 'm':
     binCenter = atof(optarg);
      break;
   case 'n':
     nsamples = atoi(optarg);
      break;
    case 'b':
      binWidth = atof(optarg);
      break;
    }

  binCenter+=0.5*binWidth;

  TFile* outfile=TFile::Open(output_file.c_str(),"RECREATE");
  TH1D* hWeights=new TH1D("hWeights","PW Weights",100,0,100);
  TTree* outtree=new TTree("pwevents","pwevents");
  double weight;
  TClonesArray* p=new TClonesArray("TLorentzVector");
  TLorentzVector beam;
  double qbeam;
  std::vector<int> q; // array of charges

  outtree->Branch("weight",&weight,"weight/d");
  outtree->Branch("p",&p);
  outtree->Branch("beam",&beam);
  outtree->Branch("q",&q);
  outtree->Branch("qbeam",&qbeam,"qbeam/i");

 
  // load integrals ---------------------------------------------------
  integral normInt; 
  ifstream intFile(integrals_file.c_str());
  if (!intFile) {
    printErr << "Cannot open file '" 
	     << integrals_file << "'. Exiting." << endl;
    throw;
  }
  // !!! integral.scan() performs no error checks!
  normInt.scan(intFile);
  intFile.close();

  
  // load production amplitudes ------------------------------------------
  // read TFitResult is used as input
  TFile* fitresults=TFile::Open(waveListFileName.c_str(),"READ");
  fitResult* Bin     = NULL;
  if (!fitresults || fitresults->IsZombie()){
    cerr << "Cannot open start fit results file " << waveListFileName << endl;
    return 1;
  }
  // get tree with start values
  bool hasfit=true;
  TTree* tree;
  fitresults->GetObject("pwa", tree);
  if (!tree)
    cerr << "Cannot find fitbin tree '"<< "pwa" << "' "<< endl;
  else {
    Bin = new fitResult();
    tree->SetBranchAddress("fitResult_v2", &Bin);
    // find entry which is closest to mass bin center
    // and has the lowest likelihood
    unsigned int iBest = 0;
    double mBest = 0;
    double loglike = 0;
    for (unsigned int i = 0; i < tree->GetEntriesFast(); ++i) {
      tree->GetEntry(i);
      if (fabs(binCenter - Bin->massBinCenter()) <= fabs(binCenter - mBest)) {
	// check also if this bin is more likely in case of many fits per bin
	//if (loglike == 0 || Bin->logLikelihood() < loglike){
	  iBest = i;
	  mBest = Bin->massBinCenter();
	  loglike = Bin->logLikelihood();
	  //}
	  //else cerr << "This is a redundant fit" << endl;
      }
    }  // end loop over TFitBins
    cerr << "mBest= " << mBest << endl;


    if(mBest<binCenter-binWidth/2. || mBest>binCenter+binWidth/2.){
      cerr << "No fit found for Mass bin m=" << binCenter << endl;
      Bin->reset();
      hasfit=false;
      return 1;
    }
    else {
      cerr << "Using data from Mass bin m=" << mBest << " bin " << iBest << endl;
      tree->GetEntry(iBest); 
    }
    // write wavelist files for generator
    // string tmpname("genamps.txt");
    // create slightly modified versions of the model according to covariances
    for(unsigned int isamples=0; isamples<nsamples; ++isamples){
      // enumerate genamps files
      TString tmpname("genamps");tmpname+=isamples;tmpname+=".txt";
      ofstream tmpfile(tmpname.Data());
      if(isamples==0)Bin->printAmpsGenPW(tmpfile);
      else {
	fitResult* clone=Bin->cloneVar();
	clone->printAmpsGenPW(tmpfile);
	delete clone;
      }
      tmpfile.close();
      waveListFileName="genamps0.txt";
    } // end loop over model versions
  }

  vector<string> waveNames;

  
  vector<vector<complex<double> >*> prodAmps(nsamples); //production amplitudes
  for(unsigned int isamples=0; isamples<nsamples; ++isamples){
    prodAmps[isamples]=new vector<complex<double> >();
  }
  
  vector<int> reflectivities;
  vector<int> ms;
  vector<int> ranks;
  int maxrank=0;
  // if we have read in a TFitResult the the name of the file has been changed!
  // so it is ok to use the same variable here. See above!

  // read in a list of waveListFiles
  for(unsigned int isamples=0; isamples<nsamples; ++isamples){
    TString tmpname("genamps");tmpname+=isamples;tmpname+=".txt";
    cerr << "Reading back file " << tmpname << endl;
    ifstream wavefile(tmpname);
  while(wavefile.good()){
    TString wavename;
    double RE, IM;
    int rank=0;
    int refl=0;
    int m=0;
    wavefile >> wavename >> RE >> IM;
    //cerr << wavename << endl;

    
    if(wavename.Length()<2)continue;
    if(RE==0 && IM==0)continue;
    // check if there is rank information
    if(wavename(0)=='V'){
      // we multiply rank by to to make space for refl+- production vectors
      rank=2*atoi(wavename(1,1).Data());
      // check reflecitivity to sort into correct production vector
      refl = wavename(9)=='+' ? 0 : 1;
      m= wavename(8)=='0' ? 0 : 1;
      //cerr << wavename(9) << endl;
      wavename=wavename(3,wavename.Length());
    }
    else if(wavename!="flat") {
      refl = wavename(6)=='+' ? 0 : 1;
      m= wavename(5)=='0' ? 0 : 1;
      //cerr << wavename(6) << endl;
    }
   

    std::complex<double> amp(RE,IM);
    prodAmps[isamples]->push_back(amp);
    cerr << wavename << " " << amp << " r=" << rank/2 
	 << " eps=" << refl 
	 << " m="   << m << endl;
    wavefile.ignore(256,'\n');
    // for first file store info on waves
    if(isamples==0){
      waveNames.push_back(wavename.Data());
      reflectivities.push_back(refl);
      ms.push_back(m);
      ranks.push_back(rank);
      if(maxrank<rank)maxrank=rank;
    }
  } // loop over file
  
    cerr << "Rank of fit was:" << maxrank+1 << endl;
  } // end loop over samples

  unsigned int nmbWaves=waveNames.size();
  // TODO reserve list of wheight branches!

  vector<ifstream*> ampfiles;
  // reserve vector beforehand because Branch
  // will take a pointer onto the elements
  //vector<double> weights((nmbWaves+1)*nmbWaves/2);
  vector<double> weights(nmbWaves);
  //unsigned int wcount=0;
  //create wheight vectors for individual intensities and interference terms
  // for(unsigned int iw=0;iw<nmbWaves;++iw){
  //   for(unsigned int jw=iw;jw<nmbWaves;++jw){
  //     TString weightname("W_");
  //     if(iw==jw)weightname+=waveNames[iw];
  //     else weightname+=waveNames[iw] +"_"+ waveNames[jw];
       
  //     outtree->Branch(weightname.Data(),&weights[wcount++],(weightname+"/d").Data());
  //   }
  // }
  for(unsigned int iw=0;iw<nmbWaves;++iw){
    TString weightname("Wintens_");
    weightname+=waveNames[iw];
    outtree->Branch(weightname.Data(),&weights[iw],(weightname+"/d").Data());
  }
   
   // create branches for the weights of the different model variants
   vector<double> modelweights(nsamples);
   for(unsigned int isamples=0;isamples<nsamples;++isamples){
     TString weightname("W");weightname+=isamples;
     outtree->Branch(weightname.Data(),&modelweights[isamples],(weightname+"/d").Data());
   }

 // open decay amplitude files --------------------------------------------
 for(unsigned int iw=0;iw<nmbWaves;++iw){
    ampfiles.push_back(new ifstream(waveNames[iw].c_str()));
  }


  // event loop ------------------------------------------------------------
  unsigned int counter = 0;
  if (evtfilename.find(".evt") > 0) {
    event e;
    list< ::particle> f_mesons;
    ifstream evtfile(evtfilename.c_str());

    while (!evtfile.eof() && evtfile.good()) {
      // evt2tree
      if (counter % 1000 == 0)
        std::cout << ".";
      if (counter++ % 10000 == 0)
        std::cout << counter;
      evtfile >> e;
      p->Delete(); // clear output arrays
      q.clear();
      f_mesons = e.f_mesons();
      fourVec pX;
      list< ::particle>::iterator it = f_mesons.begin();
      while (it != f_mesons.end()) {
        pX = it->get4P();
        new ((*p)[p->GetEntries()]) TLorentzVector(pX.x(), pX.y(), pX.z(), pX.t());
        q.push_back(it->Charge());
        ++it;
      }
      fourVec evtbeam = e.beam().get4P();
      beam.SetPxPyPzE(evtbeam.x(), evtbeam.y(), evtbeam.z(), evtbeam.t());
      qbeam = e.beam().Charge();
      

      // read decay amps for this event
      vector<complex<double> > decayamps(nmbWaves);
      for (unsigned int iw = 0; iw < nmbWaves; ++iw) {
	complex<double> decayamp;
	 ampfiles[iw]->read((char*) &decayamp, sizeof(complex<double> ));
	 decayamps[iw]=decayamp;
      }

      // weighting - do this for each model-sample
      for(unsigned int isample=0;isample<nsamples;++isample){
 
      vector<complex<double> > posm0amps(maxrank + 1); // positive refl vector m=0
      vector<complex<double> > posm1amps(maxrank + 1); // positive refl vector m=1

      vector<complex<double> > negm0amps(maxrank + 1); // negative refl vector m=0
      vector<complex<double> > negm1amps(maxrank + 1); // negative refl vector m=1

      complex<double> flatamp;

      for (unsigned int iw = 0; iw < nmbWaves; ++iw) {
        complex<double> decayamp=decayamps[iw];
        string w1 = waveNames[iw];
	if(w1=="flat"){
	  flatamp=prodAmps[isample]->at(iw);
	  continue;
	}
        //cerr << w1 << "  " << decayamp << endl;
        double nrm = sqrt(normInt.val(w1, w1).real());
        complex<double> amp = decayamp / nrm * prodAmps[isample]->at(iw);
	if(isample==0){// fill wheights of individual waves
	  weights[iw]=norm(amp);
	}
	
        if (reflectivities[iw] == 0) {
          if (ms[iw] == 0)
            posm0amps[ranks[iw]] += amp;
          else if (ms[iw] == 1)
            posm1amps[ranks[iw]] += amp;
        }
        else if(reflectivities[iw] == 1) {
          if (ms[iw] == 0)
            negm0amps[ranks[iw]] += amp;
          else if (ms[iw] == 1)
            negm1amps[ranks[iw]] += amp;
        }
      }
      // end loop over waves

      // incoherent sum over rank and diffrerent reflectivities
      weight = 0;
      if (hasfit) {
        for (int ir = 0; ir < maxrank + 1; ++ir) {
          weight += norm(posm0amps[ir]+posm1amps[ir]);
          //weight += norm(posm1amps[ir]);
          weight += norm(negm0amps[ir]+negm1amps[ir]);
	  // weight += norm(negm1amps[ir]);
        }
      }
      // add flat
      weight += norm(flatamp);


      if(isample==0)hWeights->Fill(weight);
      modelweights[isample]=weight;

      }// end loop over model samples
      outtree->Fill();

    } // end loop over events
  } // end we have an eventfile with decay amplitudes
  else {
    // Initialize PDGTable
    rpwa::particleDataTable& pdt = rpwa::particleDataTable::instance();
    pdt.readFile();

    TClonesArray prodKinPar("TObjString");
    TClonesArray prodKinMom("TVector3");
    TClonesArray decayKinPar("TObjString");
    TClonesArray decayKinMom("TVector3");
    TClonesArray *pprodKinPar = &prodKinPar;
    TClonesArray *pprodKinMom = &prodKinMom;
    TClonesArray *pdecayKinPar = &decayKinPar;
    TClonesArray *pdecayKinMom = &decayKinMom;

    TFile* infile = TFile::Open(evtfilename.c_str(), "READ");
    TTree *evt_tree;
    infile->GetObject("rootPwaEvtTree;1", evt_tree);

    evt_tree->SetBranchAddress("prodKinParticles", &pprodKinPar);
    evt_tree->SetBranchAddress("prodKinMomenta", &pprodKinMom);
    evt_tree->SetBranchAddress("decayKinParticles", &pdecayKinPar);
    evt_tree->SetBranchAddress("decayKinMomenta", &pdecayKinMom);

    unsigned int counter = 0;
    unsigned long entries = evt_tree->GetEntries();
    for (unsigned long i = 0; i < entries; i++) {
      // evt2tree
      evt_tree->GetEntry(i);
      if (counter % 1000 == 0)
        std::cout << ".";
      if (counter++ % 10000 == 0)
        std::cout << counter;
      p->Delete(); // clear output arrays
      q.clear();

      // we need 4 momentum and charge of beam and decay particles
      // in rootPwaEvtTree the particle name plus the 3 momentum is stored
      // so first get charge and particle mass from pdg table via particle name
      // then construct the 4 momenta


      // beam particle
      rpwa::particle beampar(((TString*)prodKinPar[0])->Data()); // beam particle should be at first place -> 0
      qbeam = beampar.charge();
      beam.SetVectM(*(TVector3*)prodKinMom[0], pdt.entry(beampar.bareName())->mass());

      // decay particles
      for (int dp = 0; dp < decayKinPar.GetEntries(); dp++) {
        rpwa::particle dpar(((TString*)prodKinPar[dp])->Data());
        q.push_back(dpar.charge());
        new ((*p)[p->GetEntries()]) TLorentzVector(*(TVector3*)decayKinMom[dp],
            pdt.entry(dpar.bareName())->mass());
      }

      // weighting

      vector<complex<double> > posm0amps(maxrank + 1); // positive refl vector m=0
      vector<complex<double> > posm1amps(maxrank + 1); // positive refl vector m=1

      vector<complex<double> > negm0amps(maxrank + 1); // negative refl vector m=0
      vector<complex<double> > negm1amps(maxrank + 1); // negative refl vector m=1

      for (unsigned int iw = 0; iw < nmbWaves; ++iw) {
        complex<double> decayamp;
        ampfiles[iw]->read((char*) &decayamp, sizeof(complex<double> ));
        string w1 = waveNames[iw];
        //cerr << w1 << "  " << decayamp << endl;
        double nrm = sqrt(normInt.val(w1, w1).real());
        complex<double> amp = decayamp / nrm * prodAmps[0]->at(iw);
        if (reflectivities[iw] == 1) {
          if (ms[iw] == 0)
            posm0amps[ranks[iw]] += amp;
          else if (ms[iw] == 1)
            posm1amps[ranks[iw]] += amp;
        }
        else {
          if (ms[iw] == 0)
            negm0amps[ranks[iw]] += amp;
          else if (ms[iw] == 1)
            negm1amps[ranks[iw]] += amp;
        }
      } // end loop over waves

      // incoherent sum:
      weight = 0;
      if (hasfit) {
        for (int ir = 0; ir < maxrank + 1; ++ir) {
          weight += norm(posm0amps[ir]+posm1amps[ir]);
          //weight += norm(posm1amps[ir]);
          weight += norm(negm0amps[ir]+negm1amps[ir]);
          //weight += norm(negm1amps[ir]);
        }
      }
      hWeights->Fill(weight);
      outtree->Fill();
    }
  } // end of event loop
  cout << endl << "Processed " << counter << " events" << endl;

  outfile->cd();
  hWeights->Write();
  outtree->Write();
  outfile->Close();

  for(unsigned int iw=0;iw<nmbWaves;++iw){
    ampfiles[iw]->close();
    delete ampfiles[iw];
  }
  ampfiles.clear();
  for(unsigned int isamples=0;isamples<nsamples;++isamples){
    delete prodAmps[isamples];
  }
  prodAmps.clear();

  return 0;

}

