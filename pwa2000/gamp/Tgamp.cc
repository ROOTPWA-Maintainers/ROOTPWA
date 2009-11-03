#include "Tgamp.h"

#include <cstdlib>
#include <unistd.h>
#include <iostream>

extern int keydebug;
extern particleDataTable PDGtable;


std::complex<double> 
Tgamp::Amp(unsigned int i, event& e){

  if(i>=m_waves.size()){
    std::cerr << "Invalid index! Returning (0,0)" << std::endl;
    return std::complex<double>(0,0);
  }

  keyfile keyf;
  keyf.open(m_waves[i]);
  std::complex<double> result;
  keyf.run(e,result);
  keyf.rewind();
  keyf.close();
  return result;
}
