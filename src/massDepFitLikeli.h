//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      mass dependent fit likelihood rank 1!
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef MASSDEPFITLIKELI_HH
#define MASSDEPFITLIKELI_HH

// Base Class Headers ----------------
#include "Math/IFunction.h"

// Collaborating Class Headers -------
#include <vector>
#include <map>
#include <string>
#include "pwacomponent.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "TMatrixD.h"
// Collaborating Class Declarations --
class TTree;
class TF1;

namespace rpwa {

  class fitResult;

namespace ublas = boost::numeric::ublas;
typedef std::complex<double> cnum;
typedef ublas::matrix<cnum> cmatrix;
typedef ublas::matrix<double> rmatrix;
typedef ublas::matrix<rmatrix> ccmatrix;


  class massDepFitLikeli : public ROOT::Math::IBaseFunctionMultiDim {
  public:
    
    // Constructors/Destructors ---------
    massDepFitLikeli(){}
    virtual ~massDepFitLikeli(){}
    
    
    // Accessors -----------------------
    virtual unsigned int NDim() const;
    unsigned int NDataPoints() const; /// number of data points in fit

    // Modifiers -----------------------
    void init(TTree* fitresulttree,
	      TF1* finalStatePhaseSpace,
	      pwacompset* compset,
	      double mmin=0, double mmax=5000,
	      bool doCov=true);
    


    // Operations ----------------------
    virtual double DoEval  (const double* par) const;
    
    virtual IBaseFunctionMultiDim* Clone()const {return new massDepFitLikeli(*this);}
    
    
    
  private:
    
    // Private Data Members ------------
    TTree* _tree;
    TF1* _finalStatePS;
    pwacompset* _compset;
    fitResult* _rhom; // measured spindensity matrix;
    //std::vector<fitResult> _vrhom; // memory resident data
    std::vector<cmatrix> _spindens; // mremory resident spindensity
    std::vector<ccmatrix> _cov;
    std::vector<double> _mass;

    std::vector<std::string> _wlist;
    std::vector<unsigned int> _index; // wave indices in fitResult 
    double  _mmin, _mmax; // fitrange
    bool    _doCov; // take covariance between Real and Imag into account?
    // Private Methods -----------------
    
  };
  
} // end namespace

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
