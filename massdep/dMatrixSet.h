///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2011 Sebastian Neubert
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

//-----------------------------------------------------------
//
// Description:
//    set of dMatrix amplitudes for mass dependent fit
//    describes mixing of resonances in several channels  
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#include "dMatrixAmp.h"
#include "TString.h"
#include <iostream>
#include <string>
#include <vector>
#include <map>

class TF1;

///////////////////////////////////////////////////////////////////////////////
// simple class to hold quantum numbers (integer only)
class qntag {
 public: 
  qntag(const TString& tag);
 qntag(unsigned int I, int G, unsigned int J, int P, int C) :
  _I(I),_G(G),_J(J),_P(P),_C(C){};

  friend bool operator== (const qntag& lhs,const qntag& rhs);
  friend bool operator< (const qntag& lhs,const qntag& rhs);
  friend std::ostream& operator<< (std::ostream& o,const qntag& cs);

  unsigned int I()const {return _I;}
 unsigned int G()const {return _G;}
 unsigned int J()const {return _J;}
 unsigned int P()const {return _P;}
 unsigned int C()const {return _C;}
  
 private:
  unsigned int _I;
  int _G;
  unsigned int _J;
  int _P;
  int _C;
};

///////////////////////////////////////////////////////////////////////////////

class dMatrixSet {
  public:
    dMatrixSet():_numpar(0){}
    ~dMatrixSet(){}

    void add(const qntag& igjpc, dMatrixAmp* amp){_dMatrix[igjpc]=amp;}

    // this is the overall phase space of the final state * production factors
    void setPS(TF1* fPS);

    unsigned int n() const {return _dMatrix.size();}
    unsigned int numPar() const;
    
    std::vector<std::string> wavelist()const;

    void setPar(const double* par); // set parameters
    void getPar(double* par);       // return parameters 
    unsigned int nFreePSPar() const {return _freePSpar.size();}
    double getFreePSPar(unsigned int i);
    void getFreePSLimits(unsigned int i, double& lower, double& upper);


    const dMatrixAmp* operator[](const qntag& tag) {return _dMatrix[tag];}

    friend std::ostream& operator<< (std::ostream& o,const dMatrixSet& set);
    
    double ps(double m);
    double intensity(const std::string& wave, double m);
    double phase(const std::string& wave, double m);
    double phase(const std::string& wave1,
		 const std::string& wave2,
		 double m);
    std::complex<double> overlap(const std::string& wave1,
		 const std::string& wave2,
		 double m);
    
  private:
    std::map<qntag,dMatrixAmp*> _dMatrix;
    unsigned int _numpar;
    TF1* _phasespace;
    std::vector<unsigned int> _freePSpar; // parameters of phase space to keep floating

  };
