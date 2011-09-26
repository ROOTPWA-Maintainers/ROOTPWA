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


#include <vector>
#include <complex>

#include "TObject.h"


namespace rpwa {


	class amplitudeTreeLeaf : public TObject {

	public:

		amplitudeTreeLeaf();
		virtual ~amplitudeTreeLeaf();

		void clear();

		amplitudeTreeLeaf& operator =(const amplitudeTreeLeaf& amp);

		// friend bool operator ==(const amplitudeTreeLeaf& lhsAmp,
		//                         const amplitudeTreeLeaf& rhsAmp);

		// arithmetic operators for integrals
		amplitudeTreeLeaf& operator +=(const amplitudeTreeLeaf&   amp);
		amplitudeTreeLeaf& operator -=(const amplitudeTreeLeaf&   amp);
		amplitudeTreeLeaf& operator *=(const double               factor);
		amplitudeTreeLeaf& operator /=(const double               factor);
		amplitudeTreeLeaf& operator *=(const std::complex<double> factor);
		amplitudeTreeLeaf& operator /=(const std::complex<double> factor);
		
		// accessors
		unsigned int         nmbIncohSubAmps()                             const { return _incohSubAmps.size();  }
		std::complex<double> incohSubAmp    (const unsigned int index = 0) const { return _incohSubAmps[index];  }

		void setNmbIncohSubAmps(const unsigned int         nmb)	      { _incohSubAmps.resize(nmb, 0); }
		void setIncohSubAmp    (const std::complex<double> amp,
		                        const unsigned int         index = 0) { _incohSubAmps[index] = amp;   }

		std::ostream& print(std::ostream& out) const;  ///< prints amplitudes in human-readable form

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	private:

	  static bool _debug;  ///< if set to true, debug messages are printed

		std::vector<std::complex<double> > _incohSubAmps;  ///< sub amplitudes to be added incoherently in cross section


#ifdef USE_STD_COMPLEX_TREE_LEAFS
		ClassDef(amplitudeTreeLeaf,1)
#endif

	};


	// comparison operators
	inline
	bool
	operator ==(const amplitudeTreeLeaf& lhsAmp,
	            const amplitudeTreeLeaf& rhsAmp)
	{
		const unsigned int nmbSubAmps = lhsAmp.nmbIncohSubAmps();
		if (nmbSubAmps != lhsAmp.nmbIncohSubAmps())
			return false;
		for (unsigned int i = 0; i < nmbSubAmps; ++i)
			if (lhsAmp.incohSubAmp(i) != rhsAmp.incohSubAmp(i))
				return false;
		return true;
	}

	inline
	bool
	operator !=(const amplitudeTreeLeaf& lhsAmp,
	            const amplitudeTreeLeaf& rhsAmp)
	{
		return not(lhsAmp == rhsAmp);
	}


	// arithmetic operators for integrals
	inline
	amplitudeTreeLeaf
	operator +(const amplitudeTreeLeaf& ampA,
	           const amplitudeTreeLeaf& ampB)
	{
		amplitudeTreeLeaf result = ampA;
		result += ampB;
		return result;
	}

	inline
	amplitudeTreeLeaf
	operator -(const amplitudeTreeLeaf& ampA,
	           const amplitudeTreeLeaf& ampB)
	{
		amplitudeTreeLeaf result = ampA;
		result -= ampB;
		return result;
	}

	inline
	amplitudeTreeLeaf
	operator *(const amplitudeTreeLeaf& amp,
	           const double             factor)
	{
		amplitudeTreeLeaf result = amp;
		result *= factor;
		return result;
	}

	inline
	amplitudeTreeLeaf
	operator *(const double             factor,
	           const amplitudeTreeLeaf& amp)
	{
		return amp * factor;
	}

	inline
	amplitudeTreeLeaf
	operator /(const amplitudeTreeLeaf& amp,
	           const double             factor)
	{
		amplitudeTreeLeaf result = amp;
		result /= factor;
		return result;
	}

	inline
	amplitudeTreeLeaf
	operator *(const amplitudeTreeLeaf&   amp,
	           const std::complex<double> factor)
	{
		amplitudeTreeLeaf result = amp;
		result *= factor;
		return result;
	}

	inline
	amplitudeTreeLeaf
	operator *(const std::complex<double> factor,
	           const amplitudeTreeLeaf&   amp)
	{
		return amp * factor;
	}

	inline
	amplitudeTreeLeaf
	operator /(const amplitudeTreeLeaf&   amp,
	           const std::complex<double> factor)
	{
		amplitudeTreeLeaf result = amp;
		result /= factor;
		return result;
	}


	inline
	std::ostream&
	operator <<(std::ostream&            out,
	            const amplitudeTreeLeaf& amp)
	{
		return amp.print(out);
	}


}  // namespace rpwa


#endif  // AMPLITUDETREELEAF_H
