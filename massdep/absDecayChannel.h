
///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010 Sebastian Neubert
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

/** @brief Decay Channel in isobar picture
 */

#ifndef ABSDECAYCHANNEL_HH
#define ABSDECAYCHANNEL_HH

#include <vector>
#include "TLorentzVector.h"

namespace rpwa {

  class absDecayChannel {
  public:
    absDecayChannel(double b,int l=0): _branching(b),_l(l){}
    virtual ~absDecayChannel(){}
  
    virtual void tau(std::vector<TLorentzVector>& particles,
	     double& evtweight, // returns isobar part of phase space for event
	     double& breakup) =0 ;   // returns breakup momentum of top vertex

    virtual double branching() const {return _branching;}
    virtual int l() const {return _l;}
    virtual void set_branching(double& b){_branching=b;}

  protected:

    double _branching;
    int _l;

  };



} // end namespace
  
#endif
