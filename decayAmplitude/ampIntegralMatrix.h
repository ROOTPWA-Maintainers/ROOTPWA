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
//
// Description:
//      container class for complex amplitude integral matrices
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef AMPINTEGRALMATRIX_H
#define AMPINTEGRALMATRIX_H


#include <vector>
#include <map>
#include <complex>
#include <string>
#include <iostream>
#include <fstream>


#ifndef __CINT__
#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
#endif

#include "TObject.h"

#include "waveDescription.h"


class TTree;
namespace rpwa {
	class amplitudeTreeLeaf;
}


namespace rpwa {


	class ampIntegralMatrix : public TObject {


		typedef std::map<std::string, unsigned int>::const_iterator waveNameToIndexMapIterator;


	public:


#ifndef __CINT__
		typedef boost::multi_array<std::complex<double>, 2> integralMatrixType;
	private:
		typedef integralMatrixType::size_type               sizeType;
	public:
#endif


		ampIntegralMatrix();
		ampIntegralMatrix(const ampIntegralMatrix& integral);
		virtual ~ampIntegralMatrix();

		void clear();

		ampIntegralMatrix& operator =(const ampIntegralMatrix& integral);

		friend bool operator ==(const ampIntegralMatrix& lhsInt,
		                        const ampIntegralMatrix& rhsInt);

		// arithmetic operators for integrals
		ampIntegralMatrix& operator +=(const ampIntegralMatrix& integral);
		ampIntegralMatrix& operator -=(const ampIntegralMatrix& integral);
		ampIntegralMatrix& operator *=(const double             factor);
		ampIntegralMatrix& operator /=(const double             factor);

		// accessors
		unsigned int  nmbWaves () const { return _nmbWaves;  }  ///< returns number of waves in integral
		unsigned long nmbEvents() const { return _nmbEvents; }  ///< returns number of events in integral

		void setNmbWaves (const unsigned int  nmbWaves)  { _nmbWaves  = nmbWaves;  }  ///< sets number of waves in integral
		void setNmbEvents(const unsigned long nmbEvents) { _nmbEvents = nmbEvents; }  ///< sets number of events in integral

		bool               containsWave(const std::string& waveName ) const;  ///< returns whether wave is in integral matrix
		unsigned int       waveIndex   (const std::string& waveName ) const;  ///< returns wave index for a wave name
		const std::string& waveName    (const unsigned int waveIndex) const;  ///< returns wave name for a wave index

		const std::vector<rpwa::waveDescription>& waveDescriptions() const { return _waveDescriptions; }  ///< returns array of wave descriptions
		const waveDescription* waveDesc(const unsigned int waveIndex) const;  ///< returns wave descriptuion for wave index, if existent
		const waveDescription* waveDesc(const std::string& waveName ) const { return waveDesc(waveIndex(waveName)); }  ///< returns wave descriptuion for wave name, if existent
		bool allWavesHaveDesc() const;  ///< returns whether all waves in integral have description object

#ifndef __CINT__
		integralMatrixType&       matrix()       { return _integrals; }  ///< returns integral matrix
		const integralMatrixType& matrix() const { return _integrals; }  ///< returns integral matrix
#endif

		std::complex<double> element(const unsigned int waveIndexI,
		                             const unsigned int waveIndexJ) const;  ///< returns integral matrix element divided by number of events defined by index pair
		std::complex<double> element(const std::string& waveNameI,
		                             const std::string& waveNameJ)  const  ///< returns integral matrix element divided by number of events defined by pair of wave names
		{ return element(waveIndex(waveNameI), waveIndex(waveNameJ)); }

		bool integrate(const std::vector<TTree*>&      ampTrees,
		               const std::vector<std::string>& waveNames,
		               const unsigned long             maxNmbEvents   = 0,
		               const std::string&              weightFileName = "");

		void renormalize(const unsigned long nmbEventsRenorm);

		bool writeAscii(std::ostream&      out = std::cout) const;
		bool readAscii (std::istream&      in  = std::cin );
		bool writeAscii(const std::string& outFileName) const;
		bool readAscii (const std::string& inFileName );

		std::ostream& print(std::ostream& out,
		                    const bool    printIntegralValues = false) const;  ///< prints integral in human-readable form

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	private:

		void rebuildWaveNameToIndexMap();  ///< rebuilds the wave name -> index map from _waveNames

		bool hasIdenticalWaveSet(const ampIntegralMatrix& integral) const;  ///< checks whether other integral matrix has exactly the same set of waves


		static bool _debug;  ///< if set to true, debug messages are printed

		unsigned int                        _nmbWaves;            ///< number of waves in integral
		std::map<std::string, unsigned int> _waveNameToIndexMap;  //! ///< maps wave names to wave indices
		std::vector<std::string>            _waveNames;           ///< maps wave indices to wave names
		unsigned long                       _nmbEvents;           ///< number of events in integral matrix

		void storeMultiArray();  ///< copies multiarray into storage variables written to ROOT file
		void readMultiArray ();  ///< rebuilds multiarray from storage variables read from ROOT file
#ifndef __CINT__
		integralMatrixType        _integrals;  //! ///< integral matrix
#endif
		std::vector<unsigned int> _intStorageShape;        ///< array shape
		unsigned int              _intStorageNmbElements;  ///< number of elements in array
		std::complex<double>*     _intStorageData;         //[_intStorageNmbElements]

		std::vector<rpwa::waveDescription> _waveDescriptions;  ///< wave descriptions of all waves


		ClassDef(ampIntegralMatrix,1)

	};


	// comparison operators
	inline
	bool
	operator ==(const ampIntegralMatrix& lhsInt,
	            const ampIntegralMatrix& rhsInt)
	{
		if (not lhsInt.hasIdenticalWaveSet(rhsInt)
		    or (lhsInt.nmbEvents() != rhsInt.nmbEvents()))
			return false;
		for (unsigned int i = 0; i < lhsInt.nmbWaves(); ++i) {
			const std::string waveNameI = lhsInt.waveName(i);
			for (unsigned int j = 0; j < lhsInt.nmbWaves(); ++j) {
				const std::string waveNameJ = lhsInt.waveName(j);
				if (lhsInt.matrix()[i][j]
				    != rhsInt.matrix()[rhsInt.waveIndex(waveNameI)][rhsInt.waveIndex(waveNameJ)])
					return false;
			}
		}
		return true;
	}

	inline
	bool
	operator !=(const ampIntegralMatrix& lhsInt,
	            const ampIntegralMatrix& rhsInt)
	{
		return not(lhsInt == rhsInt);
	}


	// arithmetic operators for integrals
	inline
	ampIntegralMatrix
	operator +(const ampIntegralMatrix& integralA,
	           const ampIntegralMatrix& integralB)
	{
		ampIntegralMatrix result = integralA;
		result += integralB;
		return result;
	}

	inline
	ampIntegralMatrix
	operator -(const ampIntegralMatrix& integralA,
	           const ampIntegralMatrix& integralB)
	{
		ampIntegralMatrix result = integralA;
		result -= integralB;
		return result;
	}

	inline
	ampIntegralMatrix
	operator *(const ampIntegralMatrix& integral,
	           const double             factor)
	{
		ampIntegralMatrix result = integral;
		result *= factor;
		return result;
	}

	inline
	ampIntegralMatrix
	operator *(const double             factor,
	           const ampIntegralMatrix& integral)
	{
		return integral * factor;
	}

	inline
	ampIntegralMatrix
	operator /(const ampIntegralMatrix& integral,
	           const double             factor)
	{
		ampIntegralMatrix result = integral;
		result /= factor;
		return result;
	}


	inline
	std::ostream&
	operator <<(std::ostream&            out,
	            const ampIntegralMatrix& integral)
	{
		return integral.print(out);
	}


}  // namespace rpwa


#endif  // AMPINTEGRALMATRIX_H
