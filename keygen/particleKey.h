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
// File and Version Information:
// $Id$
//
// Description:
//      
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#ifndef PARTICLEKEY_HH
#define PARTICLEKEY_HH


#include <vector>
#include <ostream>

#include "TString.h"


namespace rpwa {


  class particleKey {

  public:

    particleKey(const TString&     name,           // particle name
		particleKey*       isobar1 = 0,    // first isobar
		particleKey*       isobar2 = 0,    // second isobar (bachelor)
		int                L       = 0,    // relative orbital angular momentum between isobars
		int                S       = -1,   // total spin of isobars
		const std::string& massDep = "");  // mass dependence of amplitude (i.e. for sigma)
    virtual ~particleKey();

    // accessors
    TString name()     const { return _name; }
    TString waveName() const;
    int     L   ()     const { return _L;    }
    int     S   ()     const { return _S;    }

    void setName (const TString& name) { _name = name; }
    void setL    (const int L)         { _L    = L;    }
    void setS    (const int S)         { _S    = S;    }

    const particleKey* isobar(const unsigned int index) const { return _isobars[index]; }

    unsigned int nmbFsParticles(const int                        index = -1) const;  ///< returns number of final state particles for each isobar or total (index == -1)
    void         fsParticles   (std::vector<const particleKey*>& particles)  const;  ///< builds vector of pointers to final state particles
    void         fsCharges     (std::vector<int>&                charges)    const;  ///< extracts charge pattern of final state particles
    unsigned int countFsCharge (const int                        charge)     const;  ///< counts number of final state particles with particular charge

    void setFsIds(const std::vector<int>& fsIds);

    void write(std::ostream&      out,
	       const unsigned int offset = 0) const;

  private:

    TString      _name;               // name of particle or state
    particleKey* _isobars[2];         // pointers to the two isobars this particle decays into
    int          _L;                  // relative orbital angular momentum between isobars
    int          _S;                  // total spin of isobars
    unsigned int _nmbFsParticles[2];  // number of final state particles in iso
    bool         _isFsParticle;       // indicates whether this particle is a final state particle
    unsigned int _id;                 // ID for final state particles
    TString      _massDep;            // mass dependence of amplitude (i.e. for sigma)
  
  };


}  // namespace rpwa


#endif  // PARTICLEKEY_HH


//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
