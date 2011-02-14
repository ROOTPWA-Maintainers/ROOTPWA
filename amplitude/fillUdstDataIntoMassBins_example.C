//
// reads uDST tree and generates mass bin directory structure with .root files in ROOTPWA format
//
// $Rev:: 1907                                          $: revision of last commit
// $Author:: bgrube                                     $: author of last commit
// $Date:: 2010-08-12 19:55:14 +0200 (Thu, 12 Aug 2010) $: date of last commit
//


#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>

#include <boost/progress.hpp>

#include "TStopwatch.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

//#include "reportingUtils.hpp"


using namespace std;
using namespace boost;

// creates mass bin directory structure and returns event files for mass bins
bool createMassBinFiles(vector<TFile*>& pwaFiles, vector<TTree*>& pwaTrees, const string& dirBaseName = "/tmp/", const unsigned int nmbMassBins = 50,
    const double massBinWidth = 50, // [MeV/c^2]
    const double massRangeMin = 500, // [MeV/c^2]
    const string& treeName = "rootPwaEvtTree") {
  // cleanup
  for (unsigned int i = 0; i < pwaTrees.size(); ++i)
    if (pwaTrees[i])
      delete pwaTrees[i];
  pwaTrees.clear();
  for (unsigned int i = 0; i < pwaFiles.size(); ++i)
    if (pwaFiles[i]) {
      pwaFiles[i]->Close();
      delete pwaFiles[i];
    }
  pwaFiles.clear();
  pwaFiles.resize(nmbMassBins, 0);
  pwaTrees.resize(nmbMassBins, 0);
  bool success = true;
  for (unsigned int i = 0; i < nmbMassBins; ++i) {
    const double binMin = massRangeMin + i * massBinWidth;
    const double binMax = binMin + massBinWidth;
    // create mass bin directory
    stringstream n;
    n << binMin << "." << binMax;
    const string dirName = dirBaseName + "/" + n.str();
    mkdir(dirName.c_str(), S_IRWXU | S_IRWXG); // create directory read/writable by owner and group
    // create output file
    const string pwaFileName = dirName + "/" + n.str() + ".root";
    TFile* pwaFile = TFile::Open(pwaFileName.c_str(), "RECREATE");
    if (not pwaFile or pwaFile->IsZombie()) {
      cout << "problems creating file '" << pwaFileName << "'" << endl;
      success = false;
    }
    else {
      // pwaFile->cd();
      pwaTrees[i] = new TTree(treeName.c_str(), treeName.c_str());
      if (not pwaTrees[i]) {
        cout << "problems creating tree '" << treeName << "' " << "in file " << "'" << pwaFileName
            << "'" << endl;
        success = false;
      }
      else {
        pwaTrees[i]->SetDirectory(pwaFile);
        cout << "successfully created tree '" << treeName << "' in file " << "'" << pwaFileName
            << "'" << endl;
      }
      pwaFiles[i] = pwaFile;
    }
  }
  return success;
}

// creates mass bin directory structure and returns event files for mass bins
bool writeMassBinFile(vector<TTree*>& pwaTrees, TLorentzVector& beamLv, TLorentzVector& targetLv,
    TLorentzVector& pi01, TLorentzVector& pi02, TLorentzVector& pi_out, TLorentzVector& resonance,
    const unsigned int nmbMassBins = 50,
    const double massBinWidth = 50, // [MeV/c^2]
    const double massRangeMin = 500, // [MeV/c^2]
    const string& beamParticleName = "pi-", const string& targetParticleName = "p+",
    const bool debug = false) {
  const double mass = 1000 * resonance.M(); // convert from GeV/c^2 to MeV/c^2
  // make sure that mass is in range
  if ((mass < massRangeMin) or (mass > (massRangeMin + nmbMassBins * massBinWidth)))
    return false;
  const unsigned int bin = (unsigned int) ((mass - massRangeMin) / massBinWidth);
  if (not pwaTrees[bin]) {
    cout << "null pointer for tree for mass bin [" << massRangeMin + bin * massBinWidth << ", "
        << massRangeMin + (bin + 1) * massBinWidth << "]" << endl;
    return false;
  }
  // fill tree
  if (debug)
    cout << "filling tree for bin " << bin << " = [" << massRangeMin + bin * massBinWidth << ", "
    << massRangeMin + (bin + 1) * massBinWidth << "] MeV/c^2" << endl;

  TTree *outTree = pwaTrees[bin];
  //
  // write leafs to tree
  //
  const std::string prodKinParticlesLeafName = "prodKinParticles";
  const std::string prodKinMomentaLeafName = "prodKinMomenta";
  const std::string decayKinParticlesLeafName = "decayKinParticles";
  const std::string decayKinMomentaLeafName = "decayKinMomenta";
  // static leafs speed up reading by factor of 2
  static map<string, TClonesArray*> leafs;
  static bool firstCall = true;
  if (firstCall) {
    leafs[prodKinParticlesLeafName] = new TClonesArray("TObjString");
    leafs[prodKinMomentaLeafName] = new TClonesArray("TVector3");
    leafs[decayKinParticlesLeafName] = new TClonesArray("TObjString");
    leafs[decayKinMomentaLeafName] = new TClonesArray("TVector3");
    firstCall = false;
  }

  // clear arrays
  for (map<string, TClonesArray*>::const_iterator i = leafs.begin(); i != leafs.end(); ++i)
    i->second->Clear();

  // connect leaf variables to tree branches or create branches, if they don't exist yet
  const int split = 0;
  const int bufSize = 256000;
  for (map<string, TClonesArray*>::iterator i = leafs.begin(); i != leafs.end(); ++i)
    if (outTree->GetListOfBranches()->FindObject(i->first.c_str()))
      outTree->SetBranchAddress(i->first.c_str(), &i->second);
    else
      outTree->Branch(i->first.c_str(), "TClonesArray", &i->second, bufSize, split);

  //!!
  // beam
  new ((*leafs[prodKinParticlesLeafName])[0]) TObjString(beamParticleName.c_str());
  new ((*leafs[prodKinMomentaLeafName])[0]) TVector3(beamLv.Vect());
  // target
  new ((*leafs[prodKinParticlesLeafName])[1]) TObjString(targetParticleName.c_str());
  new ((*leafs[prodKinMomentaLeafName])[1]) TVector3(targetLv.Vect());

  // outgoing hadrons
  new ((*leafs[decayKinParticlesLeafName])[0]) TObjString("pi0");
  new ((*leafs[decayKinMomentaLeafName])[0]) TVector3(pi01.Vect());
  new ((*leafs[decayKinParticlesLeafName])[1]) TObjString("pi0");
  new ((*leafs[decayKinMomentaLeafName])[1]) TVector3(pi02.Vect());
  new ((*leafs[decayKinParticlesLeafName])[2]) TObjString("pi-");
  new ((*leafs[decayKinMomentaLeafName])[2]) TVector3(pi_out.Vect());
  /*unsigned int countHadrons = 0;
 for (int charge = -1; charge <= +1; charge += 2)
 for (unsigned int i = 0; i < nmbOfHadrons(charge); ++i) {
 new((*leafs[decayKinParticlesLeafName])[countHadrons])
 TObjString(("pi" + sign(charge)).c_str());
 new((*leafs[decayKinMomentaLeafName  ])[countHadrons])
 TVector3(hadronLv(charge, i).Vect());
 ++countHadrons;
 }*/

  // write to tree
  outTree->Fill();
  return true;
}

void fillUdstDataIntoMassBins_example(const string& inFileNamePattern =
    "/home/spflueger/eventfilters/slot2_full_selection_4to4Gammas.root",
    const string& dirName = "/home/spflueger/data/pwa/3pi_neutral/steve/",
    const long int maxNmbEvents = -1, const unsigned int nmbMassBins = 50,
    const double massBinWidth = 40, // [MeV/c^2]
    const double massRangeMin = 500, // [MeV/c^2]
    const string& uDstTreeName = "2pi0s_dM20MeV_4to4Gammas/2pi0s_dM20MeV_4to4Gammas", 
    const string& pwaTreeName = "rootPwaEvtTree", const bool debug = false) {
  TStopwatch timer;
  timer.Start();

  cout << inFileNamePattern << endl << dirName << endl << uDstTreeName << endl;
  // create chain and connect tree leaf variables to branches
  TChain uDstChain(uDstTreeName.c_str());
  uDstChain.Add(inFileNamePattern.c_str());
  const long int nmbEventsUdstChain = uDstChain.GetEntries();
  //uDstChain.GetListOfFiles()->ls();

  // connect tree leafs
  // nHadronMuoProdEvent* uDstEvent = 0;
  //uDstChain.SetBranchAddress("event", &uDstEvent);
  double t, t_prime;
  TLorentzVector g1, g2, g3, g4, pi01, pi02, pi_out, pi_in, pi_in_temp, proton, proton_rest, resonance;
  TLorentzVector *pg1 = &g1, *pg2 = &g2, *pg3 = &g3, *pg4 = &g4, *ppi_out = &pi_out, *ppi_in = &pi_in, *pproton = &proton;
  proton_rest.SetXYZM(0,0,0,0.938272);

  uDstChain.SetBranchAddress("gamma1", &pg1);
  uDstChain.SetBranchAddress("gamma2", &pg2);
  uDstChain.SetBranchAddress("gamma3", &pg3);
  uDstChain.SetBranchAddress("gamma4", &pg4);
  uDstChain.SetBranchAddress("pi_out", &ppi_out);
  uDstChain.SetBranchAddress("pi_in", &ppi_in);
  uDstChain.SetBranchAddress("proton", &pproton);

  // create directories and .root files
  cout << "creating mass bin directories ..." << endl;
  vector<TFile*> pwaDataFiles;
  vector<TTree*> pwaDataTrees;
  if (not createMassBinFiles(pwaDataFiles, pwaDataTrees, dirName, nmbMassBins, massBinWidth,
      massRangeMin, pwaTreeName)) {
    cout << "there were problems creating the mass bin files/directories. aborting." << endl;
    return;
  }
  cout << "... done. created " << pwaDataFiles.size() << " directories." << endl;

  // loop over events
  cout << "writing events into mass bin directories ..." << endl;
  const long int nmbEvents = ((maxNmbEvents > 0) ? maxNmbEvents : nmbEventsUdstChain);
  unsigned long int countEvWritten = 0;
  unsigned int cutCount =0;
  progress_display progressIndicator(nmbEvents, cout, "");
  cout << "looping over " << nmbEvents << " entries" << endl;
  for (long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex) {
    ++progressIndicator;
    if (uDstChain.LoadTree(eventIndex) < 0)
      break;
    if(uDstChain.GetEntry(eventIndex) == 0) {
	cout << "error reading from tree" << endl;
	break;
    }
    // now just make some minor calculations before
    // construct two pi0s
    pi01 = g1 + g2;
    pi02 = g3 + g4;
    // calculate t'
    resonance = pi01 + pi02 + pi_out;

    double E_beam = resonance.E() + proton.E() - proton.M();
    double scale = E_beam / pi_in.E();
    pi_in_temp = pi_in;
    pi_in_temp.SetVect(pi_in.Vect() * scale);
    pi_in_temp.SetE(E_beam);

    t = (pi_in - resonance) * (pi_in - resonance);
    t_prime = fabs(t) - fabs((resonance.M()*resonance.M() - pi_in.M()*pi_in.M())/(4.0*(pi_in.Vect()*pi_in.Vect())));
    if (0.1 > t_prime || 1.0 < t_prime) {
      cutCount++;
      continue;
    }
    if (writeMassBinFile(pwaDataTrees, pi_in_temp, proton_rest, pi01, pi02, pi_out, resonance, nmbMassBins,
        massBinWidth, massRangeMin, "pi-", "p+", debug))
      ++countEvWritten;
  }

  // write trees
  long unsigned int countTreeEvents = 0;
  for (unsigned int i = 0; i < nmbMassBins; ++i) {
    pwaDataTrees[i]->GetCurrentFile()->Write();
    long unsigned int nmbEvents = pwaDataTrees[i]->GetEntries();
    cout << "successfully written " << setw(9) << nmbEvents << " events to file " << "'"
        << pwaDataTrees[i]->GetCurrentFile()->GetName() << "'" << endl;
    countTreeEvents += nmbEvents;
    pwaDataTrees[i]->GetCurrentFile()->Close();
  }
  pwaDataFiles.clear();
  pwaDataTrees.clear();

  //printInfo << "... done. wrote " << min(countEvWritten, countTreeEvents)
  cout << "... done. " << countTreeEvents << " out of " << nmbEvents
       << " events were in mass range." << cutCount<<endl;
  timer.Stop();
  cout << "this job consumed: ";
  timer.Print();
}

// int main() {
//   fillUdstDataIntoMassBins_example();
//   return 0;
// }
