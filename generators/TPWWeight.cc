// This Class' Header ------------------
#include "TPWWeight.h"

// C/C++ Headers ----------------------
#include <list>
#include <fstream>
# include <string>

// Collaborating Class Headers --------
#include <utilities.h>

// Class Member definitions -----------



void 
TPWWeight::addWave(const std::string& keyfilename, 
		   const std::complex<double>& amp,
		   unsigned int vectori){
  
  if(vectori<=m_waves.size()){
    m_waves.resize(vectori+1);
    m_amps.resize(vectori+1);
    m_gamp.resize(vectori+1);
  }
  m_waves[vectori].push_back(keyfilename);
  m_amps[vectori].push_back(amp);

  m_gamp[vectori].addWave(keyfilename);

}


void 
TPWWeight::loadIntegrals(const std::string& normIntFileName){
  //printInfo << "Loading normalization integral from '" 
  //	    << normIntFileName << "'." << endl;
  ifstream intFile(normIntFileName.c_str());
  if (!intFile) {
    printErr << "Cannot open file '" 
	     << normIntFileName << "'. Exiting." << endl;
    throw;
  }
  // !!! integral.scan() performs no error checks!
  m_normInt.scan(intFile);
  intFile.close();
  // renomalize
  // list<string> waves=m_normInt.files();
//   list<string>::iterator it1=waves.begin();
//   list<string>::iterator it2=waves.begin();
//   while(it1!=waves.end()){
//     string w1=*it1;
//     it2=it1;
//     ++it2;
//     while(it2!=waves.end()){
//       string w2=*it2;
//       double nrm=sqrt(m_normInt.el(w1,w1).real()*m_normInt.el(w2,w2).real())/m_normInt.nevents();
//       //cout << w1 << "  " << w2 << "   " <<nrm << endl;
//       m_normInt.el(w1,w2)=m_normInt.el(w1,w2)/nrm;
//       m_normInt.el(w2,w1)=m_normInt.el(w2,w1)/nrm;
//       ++it2;
//     }
//     m_normInt.el(w1,w1)=m_normInt.nevents();
//     ++it1;
//   }
  
}





double
TPWWeight::weight(event& e){

  std::complex<double> amp(0,0);
  unsigned int nvec=m_waves.size();
  for(unsigned int ivec=0;ivec<nvec;++ivec){ // loop over production vectors
    unsigned int nwaves=m_waves[ivec].size();
    for(unsigned int iwaves=0;iwaves<nwaves;++iwaves){
      string w1=m_waves[ivec][iwaves];
      w1.erase(0,w1.find_last_of("/")+1);
      w1.replace(w1.find(".key"),4,".amp");
      double nrm=sqrt(m_normInt.val(w1,w1).real()*m_normInt.val(w1,w1).real());
      amp+=m_gamp[ivec].Amp(iwaves,e)/nrm ;
    }
  } // end loop

  return std::norm(amp);

}



