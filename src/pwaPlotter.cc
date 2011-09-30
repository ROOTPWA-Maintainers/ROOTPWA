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
#include <limits>

// Collaborating Class Headers --------
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH2D.h"
#include "TF1.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "fitResult.h"

// Class Member definitions -----------

using namespace std;
using namespace rpwa;

TH2D* drawDensity(TGraphErrors* g, TH2D* h, double weight){
  
  unsigned int ybins=h->GetNbinsY();
  unsigned int xpoints=g->GetN();
  TAxis* ax=h->GetYaxis();
  double ymin=ax->GetXmin();
  double ymax=ax->GetXmax();

  TF1 gaus("gaus","gausn(0)",ymin,ymax);
  TString tit=g->GetTitle();
  
  for(unsigned int ip=0; ip<xpoints; ++ip){
    double x=g->GetX()[ip];
    double err=g->GetEY()[ip];
    double val=g->GetY()[ip];
  
 //    if(tit.Contains("1-1++0+sigma_01_a11269") && val<100 && x>1.6 && x<2.0){
//       //cerr << "x="<<x
//       //	   << "  err="<< err 
//       //	   << "  val="<< val << endl;
//       //return h;
//     }
      
    gaus.SetParameters(1,val,err);
    for(unsigned int ibin=1; ibin<ybins;++ibin){
      double y=ax->GetBinCenter(ibin);
      //if(fabs(y-val)<5.*err){
	double w=gaus.Eval(ax->GetBinCenter(ibin))*weight;
	if(w==w)h->Fill(x,y,w);
	//}
    }
  }
  return h;
}



string
getIGJPCMEps(const std::string& wavename){
  return wavename.substr(0,7);
}



ClassImp(pwaPlotter);

pwaPlotter::pwaPlotter()
  : mMinEvidence(0)
{
  mLogLikelihood= new TMultiGraph();
  mLogLikelihood->SetTitle("LogLikelihood");
  mLogLikelihood->SetName("LogLikelihood");
  mLogLikelihoodPerEvent=new TMultiGraph();
  mLogLikelihoodPerEvent->SetTitle("LogLikelihood/Event");
  mLogLikelihoodPerEvent->SetName("LogLikelihoodPerEvent");
  mEvidence= new TMultiGraph();
  mEvidence->SetTitle("Evidence");
  mEvidence->SetName("Evidence");
  mEvidencePerEvent=new TMultiGraph();
  mEvidencePerEvent->SetTitle("Evidence/Event");
  mEvidencePerEvent->SetName("EvidencePerEvent");
  
  if(0){
  std::vector<string> waves;
  waves.push_back("1-2-+0+pi-_02_f21270=pi-+_1_a11269=pi+-_0_rho770.amp");
  waves.push_back("1-2-+0+pi-_22_f21270=pi-+_11_a11269=pi+-_01_rho770.amp");
  waves.push_back("1-2-+0+rho770_02_a21320=pi-_2_rho770.amp");
  waves.push_back("1-2-+0+rho770_02_a11269=pi-_0_rho770.amp");
  waves.push_back("1-0-+0+pi-_00_f01500=rho770_00_rho770.amp");
  waves.push_back("1-0-+0+pi-_00_f01500=sigma_0_sigma.amp");
  waves.push_back("1-0-+0+rho770_00_a11269=pi-_0_rho770.amp");
  waves.push_back("1-0-+0+rho770_22_a11269=pi-_01_rho770.amp");
  waves.push_back("1-0-+0+pi-_22_f21270=sigma_2_sigma.amp");
  waves.push_back("1-1++0+pi-_01_eta11600=pi-+_10_pi1300=pi+-_00_sigma.amp");
  waves.push_back("1-1++0+pi-_01_eta11600=pi-+_01_a11269=pi+-_01_rho770.amp");

  waves.push_back("1-1++0+sigma_01_a11269=pi-_0_rho770.amp");
  waves.push_back("1-1++0+sigma_01_a11269=pi-_1_sigma.amp");
  waves.push_back("1-1++0+rho770_11_a11269=pi-_0_rho770.amp");
  waves.push_back("1-1++0+rho770_12_a11269=pi-_0_rho770.amp");
  waves.push_back("1-1++0+rho770_01_pi1300=pi-_1_rho770.amp");
 
  waves.push_back("1-1++0+pi-_11_f11285=pi-+_11_a11269=pi+-_0_rho770.amp");
  waves.push_back("1-1++0+pi-_01_rho1600=sigma_01_rho770.amp");
  waves.push_back("1-1++0+pi-_10_f01370=rho770_00_rho770.amp");
  waves.push_back("1-1++0+pi-_10_f01500=sigma_0_sigma.amp");
  waves.push_back("1-1++0+pi-_12_b21800=pi-+_12_a21320=pi+-_21_rho770.amp");

  std::vector<strpair> phlist;
  for(unsigned int i=0;i<waves.size();++i){
    for(unsigned int j=i+1;j<waves.size();++j){
      phlist.push_back(strpair(waves[i],waves[j]));
    }
  }



  for(unsigned int i=0;i<phlist.size();++i){
    mPhases[phlist[i]]=new TMultiGraph();
    std::stringstream name;
    name << "PHI"<<phlist[i].first<<"---"<<"PHI"<<phlist[i].second;
    mPhases[phlist[i]]->SetTitle(name.str().c_str());
    mPhases[phlist[i]]->SetName(name.str().c_str());
  }
  }

}


pwaPlotter::~pwaPlotter(){
  mWavenames.clear();
  
}

void 
pwaPlotter::addFit(const std::string& filename,
		   const std::string& title,
		   const unsigned int colour,
		   const std::string& treename,
		   const std::string& branchname,
		   const unsigned int numb_bins){
  
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
  
  unsigned int ifit=mResultMetaInfo.size();
  cerr << "Adding file "<< filename << endl;
  
  intree->SetBranchAddress(branchname.c_str(),&result);
  unsigned int nbins=intree->GetEntries();
  if(numb_bins!=0 && nbins!=numb_bins){
    cerr << "Wrong number of bins "<<nbins<<" in file "
	 <<filename<<endl;
    return;
  }
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
  unsigned int numwaves=0;
  set<string> wavesinthisfit;
  
  for(unsigned int i=0;i<nbins;++i){
    intree->GetEntry(i);
    
    double massBinCenter=result->massBinCenter()*0.001;
    if(massBinCenter>mass_max)mass_max=massBinCenter;
    if(massBinCenter<mass_min)mass_min=massBinCenter;
    
    registerWave(".*"); // Total intensity
    wavesinthisfit.insert(".*");
    registerWave("^.....0"); // Total M=0
    wavesinthisfit.insert("^.....0");
    registerWave("^.....1"); // Total M=1
    wavesinthisfit.insert("^.....1");

    registerWave         ("^......\\+"); // Total Eps=+
    wavesinthisfit.insert("^......\\+");
    registerWave("^......-"); // Total Eps=-
    wavesinthisfit.insert("^......-");

    // check fitResult for used waves
    // if not already registered -> register wave (will create TMultiGraph)
    const vector<string>& waveNames=result->waveNames();
    numwaves=waveNames.size();
    // cout << "Number of Waves ="<< nwaves << endl;
    for(unsigned int iw=0;iw<numwaves;++iw){
      registerWave(waveNames[iw]);
      wavesinthisfit.insert(waveNames[iw]);
      // spin totals...
      registerWave(getIGJPCMEps(waveNames[iw]));
      wavesinthisfit.insert(getIGJPCMEps(waveNames[iw]));
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
  
  //cout << "creating graphs" << endl;
 

  // create graphs for this fit
  set<string>::iterator it=wavesinthisfit.begin();
  while(it!=wavesinthisfit.end()){
    TPwaFitGraphErrors* g = new TPwaFitGraphErrors(nbins,ifit);
    stringstream graphName;
    if(*it==".*") graphName << "g" << title << "_total";
    else if(*it=="^.....0") graphName << "g" << title << "_M0";
    else if(*it=="^.....1") graphName << "g" << title << "_M1";
    else if(*it=="^......\\+") graphName << "g" << title << "_E+";
    else if(*it=="^......-") graphName << "g" << title << "_E-";
    else graphName << "g" << title << "_" << *it;

    //cout << "creating graph   " << graphName.str() << endl;

    g->SetName (graphName.str().c_str());
    g->SetTitle(graphName.str().c_str());
    //g->SetMarkerStyle(21);
    g->SetMarkerSize(0.5);
    g->SetMarkerColor(colour);
    g->SetLineColor(colour);
    g->GetXaxis()->SetTitle("mass (GeV/c^{2})");
    g->GetYaxis()->SetTitle("intensity");
    mIntensities[*it]->Add(g,"p");
    mWaveEvidence[*it]+=evi;
    ++it;
   
  }

  //cout << "building Likelihood graphs" << endl;

  // evidence and likelihood
  TGraph* gLikeli=new TGraph(nbins);
  stringstream graphName;
  graphName << "g" << title << "_LogLikelihood";
  gLikeli->SetName (graphName.str().c_str());
  gLikeli->SetTitle(graphName.str().c_str());
  gLikeli->SetMarkerStyle(21);
  gLikeli->SetMarkerSize(0.5);
  gLikeli->SetMarkerColor(colour);
  gLikeli->SetLineColor(colour);
  mLogLikelihood->Add(gLikeli,"p");
  TGraph* gLikeliPE=new TGraph(nbins);
  graphName.clear();
  graphName << "g" << title << "_LogLikelihoodPerEvent";
  gLikeliPE->SetName (graphName.str().c_str());
  gLikeliPE->SetTitle(graphName.str().c_str());
  gLikeliPE->SetMarkerStyle(21);
  gLikeliPE->SetMarkerSize(0.5);
  gLikeliPE->SetMarkerColor(colour);
  gLikeliPE->SetLineColor(colour);
  mLogLikelihoodPerEvent->Add(gLikeliPE,"p");
  TGraph* gEvidence=new TGraph(nbins);
  graphName.clear();
  graphName << "g" << title << "_Evidence";
  gEvidence->SetName (graphName.str().c_str());
  gEvidence->SetTitle(graphName.str().c_str());
  gEvidence->SetMarkerStyle(21);
  gEvidence->SetMarkerSize(0.5);
  gEvidence->SetMarkerColor(colour);
  gEvidence->SetLineColor(colour);
  mEvidence->Add(gEvidence,"p");
  TGraph* gEvidencePE=new TGraph(nbins);
  graphName.clear();
  graphName << "g" << title << "_EvidencePerEvent";
  gEvidencePE->SetName (graphName.str().c_str());
  gEvidencePE->SetTitle(graphName.str().c_str());
  gEvidencePE->SetMarkerStyle(21);
  gEvidencePE->SetMarkerSize(0.5);
  gEvidencePE->SetMarkerColor(colour);
  gEvidencePE->SetLineColor(colour);
  mEvidencePerEvent->Add(gEvidencePE,"p");


  // create Phase graphs
  std::map<strpair,TMultiGraph*>::iterator iph=mPhases.begin();
  while(iph!=mPhases.end()){
    // check if both waves have bee used in this fit
    std::string w1=iph->first.first;
    std::string w2=iph->first.second;
    if(wavesinthisfit.find(w1)!=wavesinthisfit.end() &&
       wavesinthisfit.find(w2)!=wavesinthisfit.end()  ){
      TPwaFitGraphErrors* g = new TPwaFitGraphErrors(nbins*3,ifit);
      stringstream graphName;
      graphName << "PHI"<<w1<<"---"<<"PHI"<<w2;

      cout << "creating graph   " << graphName.str() << endl;
      
      g->SetName (graphName.str().c_str());
      g->SetTitle(graphName.str().c_str());
      g->SetMarkerStyle(21);
      g->SetMarkerSize(0.5);
      g->SetMarkerColor(colour);
      g->SetLineColor(colour);
      g->GetXaxis()->SetTitle("mass (GeV/c^{2})");
      g->GetYaxis()->SetTitle("#Delta #Phi");
      iph->second->Add(g,"p");
	
    } // endif both waves available
    ++iph;
  } // end create phase graphs







  //cout << "filling data" << endl;

  // loop again over fitResults and extract all info simultaneously
  
  for(unsigned int i=0;i<nbins;++i){
    intree->GetEntry(i);
    // loop through waves
    it=wavesinthisfit.begin();
    while(it!=wavesinthisfit.end()){
      // check if this is a single wave
      if(it->find("amp")!=it->npos){
	// check if Phase Space is already filled
	bool fillps=mPhaseSpace[*it]->GetN()!=(int)nbins;
	if(fillps){
	  unsigned int waveid=result->waveIndex(*it);
	  mPhaseSpace[*it]->Set(i+1);
	  mPhaseSpace[*it]->SetPoint(i,
				     result->massBinCenter()*0.001,
				     result->normIntegral(waveid,
							  waveid).real());
	}
      } // end ccheck for single wave

      TMultiGraph* mg=mIntensities[*it];
      TGraphErrors* g=dynamic_cast<TGraphErrors*>(mg->GetListOfGraphs()->Last());
      g->SetPoint(i,
		  result->massBinCenter()*0.001,
		  result->intensity(it->c_str()));
      g->SetPointError(i,
		       binwidth*0.5,
		       result->intensityErr(it->c_str()));
      
      ++it;
      
    }// end loop through waves

    // loop through phase plots
    iph=mPhases.begin();
    while(iph!=mPhases.end()){
      std::string w1=iph->first.first;
      std::string w2=iph->first.second;
      if(wavesinthisfit.find(w1)!=wavesinthisfit.end() &&
	 wavesinthisfit.find(w2)!=wavesinthisfit.end()  ){
	TGraphErrors* g=dynamic_cast<TGraphErrors*>(iph->second->GetListOfGraphs()->Last());

	double ph=result->phase(w1,w2);
	// check if we should make a transformation by 2pi
	// this is needed because of cyclical variable phi
	if(i>11){
	  double mpre;
	  double phpre;
	  g->GetPoint(i-3,mpre,phpre);
	  double diff1=fabs(ph-phpre);
	  double diff2=fabs(ph+360-phpre);
	  double diff3=fabs(ph-360-phpre);
	  if(diff2<diff1 && diff2<diff3)ph+=360;
	  else if(diff3<diff1 && diff3<diff2)ph-=360;
	}

	

	g->SetPoint(i*3,
		    result->massBinCenter()*0.001,
		    ph);
      g->SetPointError(i*3,
		       binwidth*0.5,
		       result->phaseErr(w1,w2));

      // add point +- 360 degree
	g->SetPoint(i*3+1,
		    result->massBinCenter()*0.001,
		    ph+360);
      g->SetPointError(i*3+1,
		       binwidth*0.5,
		       result->phaseErr(w1,w2));

	g->SetPoint(i*3+2,
		    result->massBinCenter()*0.001,
		    ph-360);
      g->SetPointError(i*3+2,
		       binwidth*0.5,
		       result->phaseErr(w1,w2));


      
      }
      ++iph;
    } // end loop over phase graphs


    gLikeli->SetPoint(i,
		      result->massBinCenter()*0.001,
		      result->logLikelihood());
    gLikeliPE->SetPoint(i,
			result->massBinCenter()*0.001,
			result->logLikelihood()/result->nmbEvents());
    gEvidence->SetPoint(i,
			result->massBinCenter()*0.001,
			result->evidence());
    gEvidencePE->SetPoint(i,
			result->massBinCenter()*0.001,
			result->evidence()/result->nmbEvents());
  }
  
  //cout << "writing meta" << endl;
  // write MetaInfo
  fitResultMetaInfo meta(filename,
			 title,
			 colour,treename,
			 branchname);
  meta.setLikelihoods(logli,logliperevt,evi,eviperevt);
  meta.setBinRange(mass_min-binwidth*0.5,mass_max+binwidth*0.5,nbins);
  meta.setNWaves(numwaves);
  mResultMetaInfo.push_back(meta);
  mMinEvidence=(mMinEvidence*ifit+evi)/(ifit+1);
    
  cout << "Fit Quality Summary: " << endl;
  cout << "  LogLikelihood:   " << logli << endl;
  cout << "  Evidence:        " << evi << endl;
  cout << "  Number of waves: " << numwaves << endl;

  // cleanup
  infile->Close();
  cerr << endl;
  
  
}



bool 
pwaPlotter::registerWave(const std::string& wavename){
  pair<set<string>::iterator,bool> inserted=mWavenames.insert(wavename);
  if(inserted.second){ // we had a true insterion
    cerr << "New wave ("<<mWavenames.size()<<"): " << wavename << endl;
    // create intensity graph:
    mIntensities[wavename]=new TMultiGraph();
    mIntensities[wavename]->SetTitle(wavename.c_str());
    mIntensities[wavename]->SetName(wavename.c_str());
    
    mPhaseSpace[wavename]=new TGraph();
    string psname=wavename;psname.append("PS");
    mPhaseSpace[wavename]->SetTitle(psname.c_str());
    mPhaseSpace[wavename]->SetName(psname.c_str());
    mWaveEvidence[wavename]=0;
  }
  
  return inserted.second;
}



void
pwaPlotter::printStats(){
  map<string,TMultiGraph*>::iterator it=mIntensities.begin();
  multimap<unsigned int, string> m;
  while(it!=mIntensities.end()){
    // get ranges
    TList* graphs=it->second->GetListOfGraphs();
    unsigned int ng=graphs->GetEntries();
    m.insert(pair<unsigned int,string>(ng,it->first));
    //cout << it->first << "    used:" << ng << endl;
    ++it;
  }
  multimap<unsigned int, string>::iterator it2=m.begin();
  
  while(it2!=m.end()){
    cout << it2->second << "    used " << it2->first << " times" 
	 << " with average evidence " << mWaveEvidence[it2->second]/(double)it2->first<< endl;
    ++it2;
  }

  double numwave=0;
  for(unsigned int i=0;i<mResultMetaInfo.size();++i){
    numwave+=mResultMetaInfo[i].nWaves;
  }
  numwave/=(double)mResultMetaInfo.size();
  cout << "Average number of waves: " << numwave << endl;
  
}





void 
pwaPlotter::produceDensityPlots(){
 map<string,TMultiGraph*>::iterator it=mIntensities.begin();
  while(it!=mIntensities.end()){
    // get ranges
    TList* graphs=it->second->GetListOfGraphs();
    unsigned int ng=graphs->GetEntries();
    double xmin=1E9;double ymin=1E9;
    double xmax=0; double ymax=-1E9;
    unsigned int nbins=0;
    for(unsigned int ig=0;ig<ng;++ig){
      TPwaFitGraphErrors* g=dynamic_cast<TPwaFitGraphErrors*>(graphs->At(ig));
      double xmin1,xmax1,ymin1,ymax1;
      g->ComputeRange(xmin1,ymin1,xmax1,ymax1);

      if(ymin > ymin1)ymin=ymin1;
      if(ymax < ymax1)ymax=ymax1;
      unsigned int ifit=g->fitindex;
      if(xmin > mResultMetaInfo[ifit].mMin)xmin=mResultMetaInfo[ifit].mMin;
      if(xmax < mResultMetaInfo[ifit].mMax)xmax=mResultMetaInfo[ifit].mMax;
      if(nbins< mResultMetaInfo[ifit].mNumBins)
	nbins=mResultMetaInfo[ifit].mNumBins;
    }
    double r=fabs(ymax-ymin)*0.1;
    // create 2D Histogram:
    string name="d";name.append(it->first);
    
    //cerr << ymin << " .. " << ymax << endl;
    TH2D* h=new TH2D(name.c_str(),name.c_str(),
		     nbins,xmin,xmax,
		     400,ymin-r,ymax+r);
    
    mIntensityDensityPlots[it->first]=h;
    
    // fill histo
    for(unsigned int ig=0;ig<ng;++ig){
      TPwaFitGraphErrors* g=dynamic_cast<TPwaFitGraphErrors*>(graphs->At(ig));
      unsigned int ifit=g->fitindex;
      double w=(mResultMetaInfo[ifit].mTotalPerEventEvidence)/(double)mResultMetaInfo[ifit].mNumBins;
      //double likeli=0;
      //cout << "weight: "<<TMath::Exp(w)<<endl;
      if(w==w)h=drawDensity(g,h,TMath::Exp(w));
    }
    ++it;
    // rescale each x-bin
    for(unsigned int ibin=1; ibin<=nbins; ++ibin){
      // get maximum bin in y for this x-bin
      double max=0;
      for(unsigned int iy=0;iy<400;++iy){
	unsigned int bin = h->GetBin(ibin,iy);
	double val=h->GetBinContent(bin);
	if(val>max)max=val;
      }
      if(max!=0 && max==max){
	for(unsigned int iy=0;iy<400;++iy){
	  unsigned int bin = h->GetBin(ibin,iy);
	  h->SetBinContent(bin,h->GetBinContent(bin)/max);
	}
      }
    }
    
    
  }// end loop over waves

}



void 
pwaPlotter::writeAll(std::string filename){
  TFile* outfile=TFile::Open(filename.c_str(),"RECREATE");
  if(outfile!=0 && !outfile->IsZombie()){
    writeAll(outfile);
    outfile->Close();
  }
  else{
    cerr << "Error opening file " << filename << endl;
  }
}

void 
pwaPlotter::writeAll(TFile* outfile){
   outfile->cd();
   // write evidence and loglikelihoods
   mLogLikelihood->Write();
   mLogLikelihoodPerEvent->Write();
   mEvidence->Write();
   mEvidencePerEvent->Write();
   writeAllIntensities(outfile);
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
  map<string,TH2D*>::iterator itd=mIntensityDensityPlots.begin();
  while(itd!=mIntensityDensityPlots.end()){
    itd->second->Write();
    ++itd;
  }
  map<string,TGraph*>::iterator itps=mPhaseSpace.begin();
  while(itps!=mPhaseSpace.end()){
    itps->second->Write();
    ++itps;
  }
  map<strpair,TMultiGraph*>::iterator itph=mPhases.begin();
  while(itph!=mPhases.end()){
    itph->second->Write();
    ++itph;
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
