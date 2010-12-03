///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
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
//    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//      TTree leaf persistency storage class for amplitude information
//      needed by fit program
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef AMPLITUDETREELEAF_H
#define AMPLITUDETREELEAF_H


// std::complex is not supported as Tree leafs in ROOT versions below 5.27.06
#include "RVersion.h"
#if ROOT_VERSION_CODE >= 334598  // make sure ROOT version is at least 5.27.06
#define AMPLITUDETREELEAF_ENABLED 1
#else
#define AMPLITUDETREELEAF_ENABLED 0
#endif


#include <vector>
#include <complex>

#include "TObject.h"


namespace rpwa {


#if AMPLITUDETREELEAF_ENABLED


	class amplitudeTreeLeaf : public TObject {
	
	public:
			
		amplitudeTreeLeaf();
    virtual ~amplitudeTreeLeaf();

		void clear();

		// accessors
		unsigned int         nmbIncohSubAmps()                             const { return _incohSubAmps.size();  }
		std::complex<double> incohSubAmp    (const unsigned int index = 0) const { return _incohSubAmps[index];  }

		void setNmbIncohSubAmps(const unsigned int         nmb)	      { _incohSubAmps.resize(nmb, 0); }
		void setIncohSubAmp    (const std::complex<double> amp,
		                        const unsigned int         index = 0) { _incohSubAmps[index] = amp;   }

		static std::string name;  ///< name of tree leaf

	private:

		std::vector<std::complex<double> > _incohSubAmps;  ///< sub amplitudes to be added incoherently in cross section

		ClassDef(amplitudeTreeLeaf,1)
    

	};


#endif  // AMPLITUDETREELEAF_ENABLED


}  // namespace rpwa


#endif  // AMPLITUDETREELEAF_H
