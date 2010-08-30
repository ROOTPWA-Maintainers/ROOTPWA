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
  class absDecayChannel;

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
  void rho(double m, std::vector<double>& results) const ; /// standard accessor: explicitely calculate
  
  void doCalc();         /// calculate and store complete phase space
  bool hasCached()const {return _graph.size()!=0;}
  double thres() const {return _thres;}
  double eval(double m, unsigned int i=0) ; /// uses cached info
  //double Evaluate(double *x, double *p) ; /// accessor used with TF1

  unsigned int nChannels()const;
  TGraph* getGraph(unsigned int i) {return _graph[i];}

  // Modifiers -----------------------
  void addDecayChannel(absDecayChannel* ch);
 



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

  std::vector<absDecayChannel*> _channels;

  // Data Points for interpolation
  std::vector<TGraph*> _graph;

  // Private Methods -----------------

};


}

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
