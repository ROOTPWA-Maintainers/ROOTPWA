#ifndef MCPHASESPACE_HH
#define MCPHASESPACE_HH

#include <vector>


// Base Class Headers ----------------


// Collaborating Class Declarations --

class TGraph;

namespace rpwa { 

  class nBodyPhaseSpaceGen;


class mcPhaseSpace {
public:

  // Constructors/Destructors ---------
  mcPhaseSpace(unsigned int n, // decay products
	       double* masses, // decay products
	       double mMin, double mMax, // mass range for n-body system
	       unsigned int nsteps,      // mass steps
	       unsigned int nsamples=10000, // samples per mass point
	       int seed=1234567);

  ~mcPhaseSpace(){}
 

  // Accessors -----------------------
  double rho(double m) const ; /// standard accessor: explicitely calculate
  double thres() const {return _thres;}
  double eval(double m) const ; /// uses cached info
  double Evaluate(double *x, double *p)const ; /// accessor used with TF1

  TGraph* getGraph() {return _graph;}

  // Modifiers -----------------------


  // Operations ----------------------

private:

  // Private Data Members ------------
  nBodyPhaseSpaceGen* _gen;
  double _wmax;
  double _thres;
  unsigned int _nsamples;
  unsigned int _nparticles;

  // Data Points for interpolation
  TGraph* _graph;

  // Private Methods -----------------

};


}

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
