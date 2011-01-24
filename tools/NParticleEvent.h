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
//-----------------------------------------------------------
// Description:
//      Utitlity class providing access to substates of an 
//      N-particle final state and GJ frames
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//-----------------------------------------------------------
#ifndef NPARTICLEEVENT_H
#define NPARTICLEEVENT_H

#include "TLorentzVector.h"
#include "NParticleState.h"
#include "TClonesArray.h"
#include <vector>
#include "FSParticle.h"

class TVector3;

class NParticleEvent {
public:
  NParticleEvent(TClonesArray* fs_momenta, 
		 std::vector<int>* fs_charges,
		 TLorentzVector* beam,
		 int* beam_charge,
		 TVector3* vertex);
  
  ~NParticleEvent(){}



  /*** @brief Refresh event
   * clears all data and builds particles new from 
   * Input arrays (can be used together with a root tree)
   * Calls build() to create substates
   */
  void refresh();
  

   /*** @brief Create all possible substates
   */
  unsigned int build();

  TLorentzVector p(); //< returns total momentum of final state
  double tprime(); //< returns effective momentum transfer


  /*** @brief transform into Gottfried Jackson frame
   */
  void toGJ();

  /*** @brief Dump evt format */
  void writeGAMP(ostream& out); //< dump event to PWA200 input format (txt file)
  unsigned int nStates() const {return _NPStates.size();} //< returns number of states
  unsigned int nParticles() const {return _fsparticles.size();} //< returns number of final state particles
  

  /*** @brief returns NParticle(Sub)State
   */
  const NParticleState& getState(unsigned int i) const {return _NPStates[i];}
  FSParticle& getParticle(unsigned int i) //< returns final state particle
  { return _fsparticles[i];}

private:
  /*** @brief Final state particle momenta
   */
  TClonesArray* _fsmomenta; 
  std::vector<int>* _fs_charges;
  TLorentzVector* _beam;
  int* _qbeam;
  TVector3* _vertex;

  /*** @brief vector to hold all the substates
   *
   * Substates are built in the refresh method and stored in this vector.
   */
  std::vector<NParticleState> _NPStates;
  std::vector<FSParticle> _fsparticles;

  /*** @brief method to create all permutations of fs particles
   * this is a magic recursive method!
   */
  void permute(int n, int k, int* permu, int x=-1, int i=1);
  
};


#endif
