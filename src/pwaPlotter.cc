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

// This Class' Header ------------------
#include "pwaPlotter.h"

// C/C++ Headers ----------------------
#include <iostream>
//#include <strstream>

// Collaborating Class Headers --------
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "fitResult.h"

// Class Member definitions -----------

using namespace std;
using namespace rpwa;

  ClassImp(pwaPlotter);


  pwaPlotter::~pwaPlotter(){
    mWavenames.clear();
    
  }
  
  void 
  pwaPlotter::addFit(const std::string& filename,
		     const std::string& title,
		     const unsigned int colour,
		     const std::string& treename,
		     const std::string& branchname){

    // Open and test file and tree
    TFile* infile = TFile::Open(filename.c_str(),"READ");
    if(infile==NULL || infile->IsZombie()){
      cerr << "Input file "<<filename<<" is not a valid file!" << endl;
      return;
    }
    TTree* intree=(TTree*)infile->FindObjectAny(treename.c_str());
    if(intree==NULL || intree->IsZombie()){
      cerr << "Tree "<<treename<<" not found in file "<<filename<< endl;
      return;
    }
    fitResult* result=0;
    if(intree->FindBranch(branchname.c_str())==NULL){
      cerr << "Invalid branch "<<treename<<"."<<branchname<<" in file "
	   <<filename<<endl;
      return;
    }
    
    intree->SetBranchAddress(branchname.c_str(),&result);
    unsigned int nbins=intree->GetEntries();
    // extract info for this fit
    // loop through bins
    // -> getRange in Mass bins
    // -> collect all used waves
    // -> integrate loglikelihood and evidence
    double mass_min=1E6;
    double mass_max=0;
    double logli=0;
    double logliperevt=0;
    double evi=0;
    double eviperevt=0;
    for(unsigned int i=0;i<nbins;++i){
      intree->GetEntry(i);

      double massBinCenter=result->massBinCenter()*0.001;
      if(massBinCenter>mass_max)mass_max=massBinCenter;
      if(massBinCenter<mass_min)mass_min=massBinCenter;

      // check fitResult for used waves
      // if not already registered -> register wave (will create TMultiGraph)
      const vector<string>& waveNames=result->waveNames();
      unsigned int nwaves=waveNames.size();
      for(unsigned int iw=0;iw<nwaves;++iw){
	registerWave(waveNames[iw]);
      }

      // get loglikelihoods
      logli+=result->logLikelihood();
      evi+=result->evidence();
      logliperevt+=result->logLikelihood()/result->nmbEvents();
      eviperevt+=result->evidence()/result->nmbEvents();
    }
    double binwidth=(mass_max-mass_min)/(double)(nbins-1);
    cerr << "Number of bins: " << nbins 
	 << "   Width: " << binwidth << endl;


    // create intensity plots ----------------------------------------------
    // We have registered all graphs in the step before...
    // This has to be done in a separate step! Try not to merge the following
    // with the loop above! You will loose generality!!! You have been warned!

    // create graphs for this fit
    set<string>::iterator it=mWavenames.begin();
    while(it!=mWavenames.end()){
      TGraphErrors* g = new TGraphErrors(nbins);
      stringstream graphName;
      graphName << "g" << title << "_" << *it;
      g->SetName (graphName.str().c_str());
      g->SetTitle(graphName.str().c_str());
      g->SetMarkerStyle(21);
      g->SetMarkerSize(0.5);
      g->SetMarkerColor(colour);
      g->SetLineColor  (colour);
      mIntensities[*it]->Add(g);
      ++it;
    }
    // loop again over fitResults and extract all info simultaneously
    for(unsigned int i=0;i<nbins;++i){
      intree->GetEntry(i);
      // loop through waves
      map<string,TMultiGraph*>::iterator it=mIntensities.begin();
      while(it!=mIntensities.end()){
	TGraphErrors* g=dynamic_cast<TGraphErrors*>(it->second->GetListOfGraphs()->Last());
	g->SetPoint(i,
		    result->massBinCenter(),
		    result->intensity(it->first.c_str()));
	g->SetPointError(i,
		    binwidth,
		    result->intensityErr(it->first.c_str()));
	++it;
      }
      
    }


    // write MetaInfo
    
    
  }
  
  

bool 
pwaPlotter::registerWave(const std::string& wavename){
  pair<set<string>::iterator,bool> inserted=mWavenames.insert(wavename);
  if(inserted.second){ // we had a true insterion
    cerr << "New wave: " << wavename << endl;
    // create intensity graph:
    mIntensities[wavename]=new TMultiGraph();
    mIntensities[wavename]->SetTitle(wavename.c_str());
  }
  
  return inserted.second;
}


void 
pwaPlotter::writeAllIntensities(std::string filename){
  TFile* outfile=TFile::Open(filename.c_str(),"RECREATE");
  if(outfile!=0 && !outfile->IsZombie()){
    writeAllIntensities(outfile);
    outfile->Close();
  }
  else{
    cerr << "Error opening file " << filename << endl;
  }
}


void 
pwaPlotter::writeAllIntensities(TFile* outfile){
  outfile->cd();
  map<string,TMultiGraph*>::iterator it=mIntensities.begin();
  while(it!=mIntensities.end()){
    it->second->Write();
    ++it;
  }
}






// void
// pwaPlotter::plotIntensity(const std::string& wavename, TTree* tr){
//   stringstream drawExpr;
//   drawExpr << branchName << ".intensity(\"" << waveName << "\"):"
// 	   << branchName << ".intensityErr(\"" << waveName << "\"):"
// 	   << branchName << ".massBinCenter() >> h" << waveName << "_" << i;
//   cout << "    running TTree::Draw() expression '" << drawExpr.str() << "' "
//        << "on tree '" << trees[i]->GetName() << "', '" << trees[i]->GetTitle() << "'" << endl;
  
//   cerr << "Drawing" << endl;
  
//   try{
//     trees[i]->Draw(drawExpr.str().c_str(), selectExpr.c_str(), "goff");
//   }
//   catch(std::exception&){
//     cerr << "Cought Exception" << endl;
//     continue;
//   }
  
  
  
//   // extract data from TTree::Draw() result and build graph
//   const int nmbBins = trees[i]->GetSelectedRows();
//   vector<double> x(nmbBins), xErr(nmbBins);
//   vector<double> y(nmbBins), yErr(nmbBins);
//   for (int j = 0; j < nmbBins; ++j) {
//     x   [j] = trees[i]->GetV3()[j] * 0.001;  // convert mass to GeV
//     xErr[j] = 0;
//     y   [j] = trees[i]->GetV1()[j] * normalization;  // scale intensities
//     yErr[j] = trees[i]->GetV2()[j] * normalization;  // scale intensity errors
//   }
  
  
  
//   TGraphErrors* g = new TGraphErrors(nmbBins,
// 				     &(*(x.begin())),      // mass
// 				     &(*(y.begin())),      // intensity
// 				     &(*(xErr.begin())),   // mass error
// 				     &(*(yErr.begin())));  // intensity error
//   {
//     stringstream graphName;
//     graphName << ((graphTitle == "") ? waveName : graphTitle) << "_" << i;
//     g->SetName (graphName.str().c_str());
//     g->SetTitle(graphName.str().c_str());
//   }
//   g->SetMarkerStyle(21);
//   g->SetMarkerSize(0.5);
//   if (graphColors) {
//     g->SetMarkerColor(graphColors[i]);
//     g->SetLineColor  (graphColors[i]);
//   }
//   graph->Add(g);
  
//   // compute maximum for y-axis
//   for (int j = 0; j < nmbBins; ++j)
//     if (maxYVal < (y[j] + yErr[j]))
//       maxYVal = y[j] + yErr[j];
//   const double yMean = g->GetMean(2);
//   if (maxYMean < yMean)
//     maxYMean = yMean;
// }
