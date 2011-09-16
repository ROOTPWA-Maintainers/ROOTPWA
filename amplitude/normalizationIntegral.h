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
//      container class for complex normalization integral matrices
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef NORMALIZATIONINTEGRAL_H
#define NORMALIZATIONINTEGRAL_H


#include <vector>
#include <map>
#include <complex>
#include <string>
#include <iostream>
#include <fstream>

#include "TObject.h"


class TTree;
namespace rpwa {
	class amplitudeTreeLeaf;
}


namespace rpwa {


	class normalizationIntegral : public TObject {


		typedef std::vector<std::vector<std::complex<double> > >    integralMatrixType;
		typedef std::map<std::string, unsigned int>::const_iterator waveNameWaveIndexMapIterator;


	public:
        
		normalizationIntegral();
		normalizationIntegral(const normalizationIntegral& integral);
		virtual ~normalizationIntegral();

		normalizationIntegral& operator =(const normalizationIntegral& integral);

		// accessors
		unsigned int  nmbWaves () const { return _nmbWaves;  }  ///< returns number of waves in integral
		unsigned long nmbEvents() const { return _nmbEvents; }  ///< returns number of events in integral

		void setNmbEvents(const unsigned long nmbEvents) { _nmbEvents = nmbEvents; }  ///< sets number of events in integral
		
		unsigned int       waveIndex(const std::string& waveName ) const;  ///< returns wave index for a wave name
		const std::string& waveName (const unsigned int waveIndex) const;  ///< returns wave name for a wave index

		std::complex<double> element(const unsigned int waveIndexI,
		                             const unsigned int waveIndexJ) const;  ///< returns integral matrix element defined by index pair
		std::complex<double> element(const std::string& waveNameI,
		                             const std::string& waveNameJ)  const;  ///< returns integral matrix element defined by pair of wave names

		bool integrate(const std::vector<std::string>& binAmpFileNames,
		               const std::vector<std::string>& rootAmpFileNames,
		               const unsigned long             maxNmbEvents   = 0,
		               const std::string&              weightFileName = "");

		void renormalize(const unsigned long nmbEventsRenorm);

		bool writeAscii(std::ostream& out = std::cout) const;
		bool readAscii (std::istream& in  = std::cin );
		bool writeAscii(const std::string& outFileName) const;
		bool readAscii (const std::string& inFileName );

		std::ostream& print(std::ostream& out,
		                    const bool    printIntegralValues = false) const;  ///< prints integral in human-readable form

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	private:

		unsigned long openBinAmpFiles(std::vector<std::ifstream*>&    ampFiles,
		                              const std::vector<std::string>& ampFileNames,
		                              const unsigned int              waveIndexOffset = 0);
		
		unsigned long openRootAmpFiles(std::vector<TTree*>&                   ampTrees,
		                               std::vector<rpwa::amplitudeTreeLeaf*>& ampTreeLeafs,
		                               const std::vector<std::string>&        ampFileNames,
		                               const unsigned int                     waveIndexOffset = 0,
		                               const std::string&                     ampLeafName     = "amplitude");
		
	  static bool _debug;  ///< if set to true, debug messages are printed

		unsigned int                        _nmbWaves;              ///< number of waves in integral
		std::map<std::string, unsigned int> _waveNameWaveIndexMap;  ///< maps wave names to wave indices
		std::vector<std::string>            _waveIndexWaveNameMap;  ///< maps wave indices to wave names
		unsigned long                       _nmbEvents;             ///< number of events in integral matrix
		integralMatrixType                  _integrals;             ///< integral matrix


#ifdef USE_STD_COMPLEX_TREE_LEAFS
		ClassDef(normalizationIntegral,1)
#endif

	};


	inline
	std::ostream&
	operator <<(std::ostream&                out,
	            const normalizationIntegral& integral)
	{
		return integral.print(out);
	}


}  // namespace rpwa


#endif  // NORMALIZATIONINTEGRAL_H
