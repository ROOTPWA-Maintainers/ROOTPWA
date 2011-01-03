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


#ifndef PWAPLOTTER_HH
#define PWAPLOTTER_HH

// Base Class Headers ----------------


// Collaborating Class Headers -------
#include <vector>
#include <map>
#include <set>
#include <string>
#include "Rtypes.h"

#include "TGraphErrors.h"

// Collaborating Class Declarations --

class TFile;
class TH2D;
class TMultiGraph;

namespace rpwa {

  typedef std::pair<std::string,std::string> strpair;


  /// \brief Meta information for one fit
  /// Data class to store information that refers to a complete set of
  /// fitResults, eg a complete mass-independent fit over several bins

  class fitResultMetaInfo {
  public:
    fitResultMetaInfo(){};
    fitResultMetaInfo(const std::string& filename,
			const std::string& title,
			unsigned int colour = 1,
			const std::string& treename="pwa",
			const std::string& branchname="fitResult_v2"):
      mfilename(filename),mtitle(title),mcolour(colour),mtreename(treename),mbranchname(branchname){}

    
    std::string mfilename;
    std::string mtitle;
    unsigned int mcolour;
    std::string mtreename;
    std::string mbranchname;


    void setLikelihoods(double logli,double logliperevt,
			double evi, double eviperevt){
      mTotalLogLikelihood=logli;
      mTotalPerEventLogLikelihood=logliperevt;
      mTotalEvidence=evi;
      mTotalPerEventEvidence=eviperevt;
    }
    double mTotalLogLikelihood;
    double mTotalPerEventLogLikelihood;
    double mTotalEvidence;
    double mTotalPerEventEvidence;

    void setBinRange(double min, double max, unsigned int n){
      mMin=min; mMax=max; mNumBins=n;
    }
    double mMin; // minimum mass bin in fit
    double mMax; // maximum mass bin in fit
    unsigned int mNumBins;
   
    

  };


  /// \brief Decorator for TGraph to carry index to fit 
  /// this is needed for density calculation
  class TPwaFitGraphErrors : public TGraphErrors {
  public:
    TPwaFitGraphErrors() : TGraphErrors(), fitindex(0){}
    TPwaFitGraphErrors(unsigned int n, unsigned int i): TGraphErrors(n),fitindex(i){}
    virtual ~TPwaFitGraphErrors(){};
    unsigned int fitindex;
    
    ClassDef(TPwaFitGraphErrors,1);
  };


  /// \brief Plot generator Class reads fitResult trees 
  ///  and creates graphs even for larger number of fits
class pwaPlotter {
public:

  // Constructors/Destructors ---------
  pwaPlotter();
  virtual ~pwaPlotter();

  // Accessors -----------------------
  const std::set<std::string>& wavesNames(){return mWavenames;} 
  const std::set<std::string>& listJPCME(){return mJPCME;} 


  // Modifiers -----------------------
  /// \brief Main function to add information to plotter
  /// title will be used as a prefix for the graphs from this fit

  void addFit(const std::string& filename,
	      const std::string& title,
	      const unsigned int colour=1,
	      const std::string& treename="pwa",
	      const std::string& branchname="fitResult_v2",
	      const unsigned int numb_bins=0);
	      


  /// \brief Create 2D density plots of the intensities
  /// This will produce a wheighted probability density profile
  /// combining the information of all the fits added
  void produceDensityPlots(); 

  // Operations ----------------------
  void writeAllIntensities(std::string filename);
  void writeAllIntensities(TFile* outfile);

//   void writeSpinTotals(std::string filename, std::string opt="");
//   void writeSpinTotals(TFile* outfile, std::string opt="");
  
  void writeAll(std::string filename);
  void writeAll(TFile* outfile);

  void printStats();
  

private:

  // Private Data Members ------------
  std::set<std::string> mWavenames;  ///< list of wave names
  std::set<std::string> mJPCME;    ///< list of available spin totals
  std::vector<fitResultMetaInfo> mResultMetaInfo; ///< list of available fits
  
  ///< 2D-probability distributions (TH2D)
  std::map<std::string,TH2D*> mIntensityDensityPlots; 
 
  
  ///< TMultiGraphs
  std::map<std::string,TMultiGraph*> mIntensities;
  std::map<strpair,TMultiGraph*> mPhases;

  std::map<std::string,TGraph*> mPhaseSpace;
  std::map<std::string,double> mWaveEvidence;

  TMultiGraph* mLogLikelihood;
  TMultiGraph* mLogLikelihoodPerEvent;
  TMultiGraph* mEvidence;
  TMultiGraph* mEvidencePerEvent;
  

  double mMinEvidence; // minimal evidence for a fit -- used to renormalize weighting

  // Private Methods -----------------
  
  bool registerWave(const std::string& wavename); ///< create wave histograms/graphs

  ClassDef(pwaPlotter,2);

};

} // end namespace

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
