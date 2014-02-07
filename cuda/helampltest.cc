// #include <fstream>

#include <boost/progress.hpp>

// #include "TVector3.h"
// #include "TLorentzRotation.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
// #include "TFile.h"
// #include "TH1.h"
// #include "TH2.h"
// #include "TSystem.h"

// #include "mathUtils.hpp"
// #include "reportingUtilsRoot.hpp"
// #include "conversionUtils.hpp"
#include "particleDataTable.h"
// #include "diffractiveDissVertex.h"
// #include "massDependence.h"
#include "waveDescription.h"
#include "isobarAmplitude.h"
#include "isobarHelicityAmplitude.h"
// #include "isobarCanonicalAmplitude.h"
// #include "evtTreeHelper.h"


using namespace std;
using namespace boost;
using namespace rpwa;


int
main(int argc, char** argv)
{
  printCompilerInfo();
  printGitHash();
}
