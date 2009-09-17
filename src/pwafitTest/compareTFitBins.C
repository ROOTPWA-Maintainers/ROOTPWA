#include <iostream>
#include <iomanip>
#include <string>
#include <cassert>
#include <algorithm>
#include <map>

#include "TFile.h"
#include "TTree.h"

#include "../TFitBin.h"


using namespace std;


bool
binsEqual(TFitBin&             binA,
	  TFitBin&             binB,
	  const double         threshold,
	  map<string, double>& maxDeviations)
{
  if (binA.mass() != binB.mass()) {
    cout << "bins have different mass: " << binA.mass() << " vs. " << binB.mass() << endl;
    return false;
  }

  if (binA.rawEvents() != binB.rawEvents()) {
    cout << "bins have different number of events: " << binA.rawEvents() << " vs. " << binB.rawEvents() << endl;
    return false;
  }

  if (binA.nwaves() != binB.nwaves()) {
    cout << "bins have different number of waves: " << binA.nwaves() << " vs. " << binB.nwaves() << endl;
    return false;
  }
  const int nmbWaves = binA.nwaves();
  bool binsAreEqual = true;
  for (int i = 0; i < nmbWaves; ++i)
    if (binA.waveDesignator(i) != binB.waveDesignator(i)) {
      cout << "bins have different wave designator at index " << setw(2) << i << ": " << binA.waveDesignator(i) << " vs. " << binB.waveDesignator(i) << endl;
      binsAreEqual = false;
    }
  if (!binsAreEqual)
    return false;

  const double likelihoodDiff = binA.logli() - binB.logli();
  if (likelihoodDiff > threshold) {
    cout << "likelihood difference " << likelihoodDiff << " is above threshold." << endl;
    //return false;
  }

  // compare amplitudes
  if (binA.namps() != binB.namps()) {
    cout << "bins have different number of amplitudes: " << binA.namps() << " vs. " << binB.namps() << endl;
    //return false;
  }
  // const int nmbAmps = binA.namps();
//   maxDeviations["amplitude"] = 0;
//   for (int i = 0; i < nmbAmps; ++i) {
//     const TComplex diff = binA.amp(i) - binB.amp(i);
//     if (diff.Re() > maxDeviations["amplitude"])
//       maxDeviations["amplitude"] = diff.Re();
//     if (diff.Im() > maxDeviations["amplitude"])
//       maxDeviations["amplitude"] = diff.Im();
//     if ((diff.Re() > threshold) || (diff.Im() > threshold)) {
//       cout << " --- "<< binA.wavename(i) << endl;
//       cout << " --- "<< binB.wavename(i) << endl;
//       cout << "amplitude difference " << setw(2) << i << " is above threshold: " << binA.amp(i) << " - " << binB.amp(i) << " = " << diff << endl;
//       binsAreEqual = false;
//     }
//   }

  // compare intensities
  maxDeviations["intensity"] = 0;
  for (int i = 0; i < nmbWaves; ++i) {
    const double diff = binA.intens(i) - binB.intens(binA.wavetitle(i));
    if (diff > maxDeviations["intensity"])
      maxDeviations["intensity"] = diff;
    if (diff > threshold)
    {
      cout << " --- "<< binA.wavetitle(i) << endl;
      cout << " --- "<< binB.wavetitle(i) << endl;
      cout << "intensity difference " << setw(2) << i << " is above threshold: " << binA.intens(i) << " - " << binB.intens(i) << " = " << diff << endl;
      binsAreEqual = false;
    }
  }
  // compare intensity errors
  maxDeviations["intensity error"] = 0;
  for (int i = 0; i < nmbWaves; ++i) {
    const double diff = binA.err(i) - binB.err(binA.wavetitle(i));
    if (diff > maxDeviations["intensity error"])
      maxDeviations["intensity error"] = diff;
    if (diff > threshold) {
      cout << "difference of intensity error " << setw(2) << i << " is above threshold: " << binA.err(i) << " - " << binB.err(binA.wavetitle(i)) << " = " << diff << endl;
      binsAreEqual = false;
    }
  }

  if (!binsAreEqual)
    return false;
  else
    return true;
}


bool
compareTFitBins(const string& fileNameA = "result.2000.root",
		const string& fileNameB = "result.root",
		const double  threshold = 0.5)
{
  const string fileNames[2] = {fileNameA,
			       fileNameB};
  const string treeName     = "pwa";
  const string branchName   = "fitbin";

  cout << "comparing TFitBin trees in " << fileNames[0] << " and " << fileNames[1] << "." << endl;

  // open files
  TFile* f[2];
  for (unsigned int i = 0; i < 2; ++i) {
    f[i] = TFile::Open(fileNames[i].c_str(), "READ");
    if (!f[i] || f[i]->IsZombie()) {
      cerr << "cannot open file '" << fileNames[i] << "'." << endl;
      return false;
    }
  }

  // read in fit bins
  TFitBin* b[2];
  TTree*   t[2];
  for (unsigned int i = 0; i < 2; ++i) {
    f[i]->GetObject(treeName.c_str(), t[i]);
    if (!t[i]) {
      cerr << "cannot find tree '"<< treeName << "' in file '" << fileNames[i] << "'." << endl;
      return false;
    } else {
      b[i] = new TFitBin();
      t[i]->SetBranchAddress(branchName.c_str(), &b[i]);
    }
  }
 
  // compare tree entries
  map<string, double> maxDeviations;
  bool                treesEqual    = true;
  const unsigned int  nmbEntries[2] = {t[0]->GetEntries(),
				       t[1]->GetEntries()};
  if (nmbEntries[0] != nmbEntries[1])
    cerr << "warning: trees have not the same number of entries: " << nmbEntries[0] << " != " << nmbEntries[1] << ". comparing up to smaller index." << endl;
  const unsigned int maxNmbEntries = min(nmbEntries[0], nmbEntries[1]);
  for (unsigned int i = 0; i < maxNmbEntries; ++i) {
    t[0]->GetEntry(i);
    t[1]->GetEntry(i);
    map<string, double> deviations;
    if (!binsEqual(*b[0], *b[1], threshold, deviations))
      treesEqual = false;
    for (map<string, double>::const_iterator i = deviations.begin(); i != deviations.end(); ++i)
      if ((maxDeviations.find(i->first) == maxDeviations.end()) || (i->second > maxDeviations[i->first]))
	maxDeviations[i->first] = i->second;
  }

  cout << "maximum deviations:" << endl;
  for (map<string, double>::const_iterator i = maxDeviations.begin(); i != maxDeviations.end(); ++i)
    cout << "    " << i->first << " = " << i->second << "; threshold = " << threshold << endl;

  // cleanup
  for (unsigned int i = 0; i < 2; ++i)
    if (f[i]) {
      f[i]->Close();
      delete f[i];
    } 
  return treesEqual;
}
