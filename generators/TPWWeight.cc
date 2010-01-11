// This Class' Header ------------------
#include "TPWWeight.h"

// C/C++ Headers ----------------------
#include <list>
#include <fstream>
# include <string>

// Collaborating Class Headers --------
#include <utilities.h>
#include <particle.h>

// Class Member definitions -----------



void 
TPWWeight::addWave(const std::string& keyfilename, 
		   TProductionAmp* amp,
		   const std::complex<double>& branching,
		   unsigned int vectori){
  
  if(vectori<=m_waves.size()){
    m_waves.resize(vectori+1);
    m_amps.resize(vectori+1);
    m_branchings.resize(vectori+1);
    m_gamp.resize(vectori+1);
  }

  // check if this is a double amplitude (by isospin symmetry)
  if(keyfilename.find("+-",7)!=string::npos){
    cerr << "Decomposing " << keyfilename << endl;
    cerr << "What is the relative phase (0 or pi)?:";
    //string select;
    //cin >> select;
    if(keyfilename.find("f11285",7)!=string::npos)m_relphase[keyfilename]=-1;
    else m_relphase[keyfilename]=1;

    // deccompose into components:
    TString name(keyfilename.c_str());
    int pos=name.Last('/');
    TString w1a=name(pos+8,name.Length()-pos-8);
   
    w1a.ReplaceAll("+-","+");w1a.ReplaceAll("-+","-");
    w1a.Prepend(name(0,pos+8));

    TString w1b=name(pos+8,name.Length()-pos-8);
    w1b.ReplaceAll("+-","-");w1b.ReplaceAll("-+","+");
    w1b.Prepend(name(0,pos+8));

    cerr << w1a << endl;
    cerr << m_relphase[keyfilename] << endl;
    cerr << w1b << endl;

    m_gamp[vectori].addWave(w1a.Data());
    m_waves[vectori].push_back(keyfilename);
    m_branchings[vectori].push_back(branching);
    m_amps[vectori].push_back(amp);
    m_gamp[vectori].addWave(w1b.Data());
    m_waves[vectori].push_back(keyfilename);
    m_branchings[vectori].push_back(branching);
    m_amps[vectori].push_back(amp);

  }
  else {
    m_gamp[vectori].addWave(keyfilename);
    m_waves[vectori].push_back(keyfilename);
    m_branchings[vectori].push_back(branching);
    m_amps[vectori].push_back(amp);
  }
  

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
  m_hasInt=true;
}



std::complex<double> 
TPWWeight::prodAmp(unsigned int iv, 
		   unsigned int iw,
		   event& e){
  // get mass
  std::list<particle> part=e.f_particles();
  fourVec p;
  std::list<particle>::iterator it=part.begin();
  while(it!=part.end())p+=(it++)->get4P();
  double m=p.len();
  return m_amps[iv][iw]->amp(m)*m_branchings[iv][iv];
}





double
TPWWeight::weight(event& e){
   unsigned int nvec=m_waves.size();
  double w=0;
  for(unsigned int ivec=0;ivec<nvec;++ivec){ // loop over production vectors
    std::complex<double> amp(0,0);
    unsigned int nwaves=m_waves[ivec].size();
    for(unsigned int iwaves=0;iwaves<nwaves;++iwaves){
      string w1=m_waves[ivec][iwaves];
      
      // std::cerr << w1 << std::endl;
      //std::cerr.flush();
      w1.erase(0,w1.find_last_of("/")+1);
      w1.replace(w1.find(".key"),4,".amp");
      std::complex<double> decayamp;
      decayamp=m_gamp[ivec].Amp(iwaves,e);
      //std::cerr << "..first component ..."  << std::endl;
      if(w1.find("+-",7)!=string::npos){ // brute force treatment of composed amplitudes!

	// see addWave to understand bookkeeping!
	// Here an isospin clebsch factor is missing
	++iwaves;
	decayamp+=m_relphase[m_waves[ivec][iwaves]]*m_gamp[ivec].Amp(iwaves,e);
      }
      // std::cerr << "..done"  << std::endl;
      double nrm=1;
      if(m_hasInt)nrm=sqrt(m_normInt.val(w1,w1).real());
      amp+=decayamp/nrm*prodAmp(ivec,iwaves,e);
    }
    w+=std::norm(amp);
  } // end loop over production vectors

  return w;

}



