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


// std::complex is not supported as Tree leafs in ROOT versions below 5.27.06
#include "RVersion.h"
#if ROOT_VERSION_CODE >= 334598  // make sure ROOT version is at least 5.27.06
#define NORMALIZATIONINTEGRAL_ENABLED 1
#else
#define NORMALIZATIONINTEGRAL_ENABLED 0
#endif


#include <vector>
#include <complex>
#include <map>
#include <string>

#include "TObject.h"


namespace rpwa {


#if NORMALIZATIONINTEGRAL_ENABLED


	class normalizationIntegral : public TObject {


		typedef std::vector<std::vector<std::complex<double> > >    integralMatrixType;
		typedef std::map<std::string, unsigned int>::const_iterator waveNameIndexMapIterator;


	public:
        
		normalizationIntegral();
		// normalizationIntegral(char**          files);
		normalizationIntegral(const normalizationIntegral& integral);
		virtual ~normalizationIntegral();

		normalizationIntegral& operator =(const normalizationIntegral& integral);

		// accessors
		std::string  weightFileName() const { return _weightFileName; }
		unsigned int maxNmbEvents  () const { return _maxNmbEvents;   }
		unsigned int nmbWaves      () const { return _nmbWaves;       }
		unsigned int nmbEvents     () const { return _nmbEvents;      }
		
		void setWeightFileName(const std::string& fileName)     { _weightFileName = fileName;     }
		void setMaxNmbevents  (const unsigned int maxNmbEvents) { _maxNmbEvents   = maxNmbEvents; }

		unsigned int waveIndex(const std::string& waveName) const;

		const std::complex<double>& element(const unsigned int waveIndexI,
		                                    const unsigned int waveIndexJ) const;
		const std::complex<double>& element(const std::string& waveNameI,
		                                    const std::string& waveNameJ)  const;


		// normalizationIntegral& files(char**                        files);
		// normalizationIntegral& files(const std::list<std::string>& files);
		// std::list<std::string> files()       const;
		// char**                 files_c_str() const;

		bool integrate();
		// normalizationIntegral& integrate();
		// normalizationIntegral& renormalize(const int n);
		// normalizationIntegral& events(const int n);

		// std::complex<double> val(const std::string& iName,
		//                          const std::string& jName);

		// normalizationIntegral get(char**                        files);
		// normalizationIntegral get(const std::list<std::string>& files);


		// matrix<std::complex<double> > mat();

		bool writeAscii(std::ostream& out = std::cout) const;
		bool readAscii (std::istream& in  = std::cin );

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	private:

	  static bool _debug;  ///< if set to true, debug messages are printed

		std::string  _weightFileName;  ///< path to file with weights used in MC importance sampling
		unsigned int _maxNmbEvents;    ///< maximum number of events to process

		unsigned int                        _nmbWaves;          ///< number of waves in integral
		std::map<std::string, unsigned int> _waveNameIndexMap;  ///< maps wave names to wave indices
		unsigned int                        _nmbEvents;         ///< number of events in integral matrix
		integralMatrixType                  _integrals;         ///< integral matrix


		ClassDef(normalizationIntegral,1)

	};


#endif  // NORMALIZATIONINTEGRAL_ENABLED


}  // namespace rpwa


#endif  // NORMALIZATIONINTEGRAL_H
