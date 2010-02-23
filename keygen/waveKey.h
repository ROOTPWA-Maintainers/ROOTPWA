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


#ifndef WAVEKEY_HH
#define WAVEKEY_HH


#include <vector>
#include <ostream>

#include "TString.h"

#include "particleKey.h"


namespace rpwa {


  class waveKey {

  public:

    waveKey(particleKey* mother,   // decay chain of mother particle
	    const int    J,        // total angular momentum
	    const int    P,        // parity
	    const int    M,        // z-projection of J
	    const int    refl,     // reflectivity
	    const int    I =  1,   // isospin
	    const int    G = -1,   // G-parity
	    const int    C = +1);  // charge conjugation

    // accessors
    int                I           () const { return _I;            }
    int                G           () const { return _G;            }
    int                J           () const { return _J;            }
    int                P           () const { return _P;            }
    int                C           () const { return _C;            }
    int                M           () const { return _M;            }
    int                refl        () const { return _refl;         }
    int                L           () const { return _mother->L();  }
    int                S           () const { return _mother->S();  }
    const particleKey& mother      () const { return *_mother;      }
    std::string        outputFormat() const { return _outputFormat; }

    TString waveName(const bool file = false) const;

    void setOutputFormat(const std::string& outputFormat = "binary") { _outputFormat = outputFormat; }

    // constructs reflectivity eigenstate and writes it out
    void write(std::ostream&      out,
	       const std::string& outputFormat = "") const;
    void write(const std::string& outputFormat = "") const;

  private:

    void setStateName(const int J,
		      const int P,
		      const int M) const;

    // assigns FS particle IDs according to the per-charge IDs
    void assignFsIds(std::vector<int>&       fsIds,
		     const std::vector<int>& fsNegIds,
		     const std::vector<int>& fsPosIds,
		     const std::vector<int>& fsNeutIds,
		     const std::vector<int>& fsCharges) const;

    void writeReflectivityBasisAmp(std::ostream&      out,
				   const unsigned int offset = 0) const;

    // constructs amplitude for a particular M and writes it out
    void writeBoseSymmetrizedAmp(std::ostream&      out,
				 const unsigned int offset = 0) const;

    int          _I;             // isospin
    int          _G;             // G-parity
    int          _J;             // total angular momentum
    int          _P;             // parity
    int          _C;             // charge conjugation
    int          _M;             // z-projection of J; since reflectivity basis is used M >= 0
    int          _refl;          // reflectivity
    particleKey* _mother;        // decay chain of mother particle
    std::string  _outputFormat;  // output format of amplitudes: 'binary', 'ascii', or 'none'
    int          _waveDebug;     // value of 'debug' parameter written to key file

  };


}  // namespace rpwa


#endif  // WAVEKEY_HH


//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
