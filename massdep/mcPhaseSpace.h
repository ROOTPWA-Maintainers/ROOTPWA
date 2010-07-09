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

/** @brief Monte Carlo Phase Space
 */

#ifndef MCPHASESPACE_HH
#define MCPHASESPACE_HH

#include <vector>


// Base Class Headers ----------------


// Collaborating Class Declarations --

class TGraph;


namespace rpwa { 

  class nBodyPhaseSpaceGen;
  class absMassDep;

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
  double rho(double m, int l=0) const ; /// standard accessor: explicitely calculate
  
  void doCalc(int l);         /// calculate and store complete phase space
  bool hasCached()const {return _graph!=NULL;}
  double thres() const {return _thres;}
  double eval(double m) ; /// uses cached info
  double Evaluate(double *x, double *p) ; /// accessor used with TF1

  TGraph* getGraph() {return _graph;}

  // Modifiers -----------------------

  void setSubSystems21(rpwa::absMassDep* iso); // 3-body decay (pp)p opt=1
  void setSubSystems22(rpwa::absMassDep* iso1, 
		       rpwa::absMassDep* iso2);  // 4-body decay (pp)(pp) opt=2
  void setSubSystems132(rpwa::absMassDep* iso1, 
			rpwa::absMassDep* iso2); // 4-body decay  p((pp)p) opt=3



  // Operations ----------------------

private:

  // Private Data Members ------------
  nBodyPhaseSpaceGen* _gen;
  double _wmax;
  double _thres;
  double _mMax;
  unsigned int _nsteps;
  unsigned int _nsamples;
  unsigned int _nparticles;

  std::vector<absMassDep*> _isobars;
  int _opt; // which isobar option to be used

  // Data Points for interpolation
  TGraph* _graph;

  // Private Methods -----------------

};


}

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
