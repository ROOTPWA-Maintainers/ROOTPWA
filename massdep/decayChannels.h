
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

#ifndef DECAYCHANNELS_HH
#define DECAYCHANNELS_HH

#include "absDecayChannel.h"
#include "TLorentzVector.h"

namespace rpwa {
  
  class absMassDep;

  // 3body X->(12)3
  class decay21 : public absDecayChannel {
  public:
    decay21(absMassDep* iso, int l, double b=1) :
      absDecayChannel(b,l),isobar(iso){}
    
    virtual void tau( std::vector<TLorentzVector>& particles,
	     double& evtweight, // returns isobar part of phase space for event
	     double& breakup);   // returns breakup momentum of top vertex
 
  private:
    absMassDep* isobar;


  };

  // 4body X->(12)(34)
  class decay22 : public absDecayChannel {
  public:
    decay22(absMassDep* iso1,
	    absMassDep* iso2,
	    int l,
	    double b=1) :
      absDecayChannel(b,l),isobar1(iso1),isobar2(iso2){}
    
    virtual void tau( std::vector<TLorentzVector>& particles,
	     double& evtweight, // returns isobar part of phase space for event
	     double& breakup);   // returns breakup momentum of top vertex
 
  private:
    absMassDep* isobar1;
    absMassDep* isobar2;
    

  };


  // 4body X->((12)3)4
  class decay23 : public absDecayChannel {
  public:
    decay23(absMassDep* iso2,
	    absMassDep* iso3,
	    int l,
	    double b=1) :
      absDecayChannel(b,l),isobar2(iso2),isobar3(iso3){}
    
    virtual void tau( std::vector<TLorentzVector>& particles,
	     double& evtweight, // returns isobar part of phase space for event
	     double& breakup);   // returns breakup momentum of top vertex
 
  private:
    absMassDep* isobar2;
    absMassDep* isobar3;
    

  };


} // end namespace

#endif
















