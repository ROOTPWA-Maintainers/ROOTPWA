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
//      TTree leaf storage class for amplitude information
//
//      class allows to define named subamplitudes that are summed
//      incoherently in the intensity formula
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
#include <map>
#include <complex>

#include "TObject.h"


namespace rpwa {


	class amplitudeTreeLeaf : public TObject {


		typedef std::map<std::string, unsigned int>::const_iterator labelToIndexMapIterator;


	public:

		amplitudeTreeLeaf();
		virtual ~amplitudeTreeLeaf();

		void clear();  ///< clears all subamps and subamp labels

		amplitudeTreeLeaf& operator =(const amplitudeTreeLeaf& amp);

		friend bool operator ==(const amplitudeTreeLeaf& lhsAmp,
		                        const amplitudeTreeLeaf& rhsAmp);

		// arithmetic operators for integrals
		amplitudeTreeLeaf& operator +=(const amplitudeTreeLeaf& amp);  ///< adds all subamps of two amplitudeTreeLeafs with same set of subamps
		amplitudeTreeLeaf& operator -=(const amplitudeTreeLeaf& amp);  ///< subtracts all subamps of two amplitudeTreeLeafs with same set of subamps

		template<typename T>
		amplitudeTreeLeaf& operator *=(const T& factor);  ///< muliplies all subamps with factor
		template<typename T>
		amplitudeTreeLeaf& operator /=(const T& factor);  ///< divides all subamps by factor

		// accessors
		unsigned int nmbIncohSubAmps() const { return _incohSubAmps.size(); }  ///< returns number of incoherent subamps

		bool               containsIncohSubAmp(const std::string& subAmpLabel) const;
		unsigned int       incohSubAmpIndex   (const std::string& subAmpLabel) const;
		const std::string& incohSubAmpName    (const unsigned int subAmpIndex) const;

		const std::complex<double>& incohSubAmp(const unsigned int index = 0  ) const
		{ return _incohSubAmps[index];                         }  ///< returns incoherent subamp at index
		const std::complex<double>& incohSubAmp(const std::string& subAmpLabel) const
		{ return _incohSubAmps[incohSubAmpIndex(subAmpLabel)]; }  ///< returns incoherent subamp at index

		const std::complex<double>& amp() const	{ return incohSubAmp(0);       }  ///< returns first incoherent subamp (meant for cases where there is only one amplitude)

		void defineIncohSubAmps(const std::vector<std::string>& subAmpLabels);  ///< defines number of subamps for this amplitude and their labels; vector has to be more than one entry

		void setIncohSubAmp(const std::complex<double>& amp,
		                    const unsigned int          index = 0) { _incohSubAmps[index] = amp; }  ///< sets incoherent subamp defined by index

		void setAmp (const std::complex<double>& amp) { setIncohSubAmp(amp, 0); }  ///< returns first incoherent subamp (meant for cases where there is only one amplitude)

		std::ostream& print(std::ostream& out) const;  ///< prints amplitudes in human-readable form

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag

		void rebuildSubAmpLabelMap();  ///< rebuilds label-to-index map for subamps


	private:

		static bool _debug;  ///< if set to true, debug messages are printed

		std::vector<std::complex<double> >  _incohSubAmps;       ///< sub-amplitudes to be added incoherently in intensity formula
		std::vector<std::string>            _incohSubAmpLabels;  ///< labels for subamps
		std::map<std::string, unsigned int> _labelToIndexMap;    //! ///< maps subamp labels to indices


		ClassDef(amplitudeTreeLeaf,1)

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
	template<typename T>
	inline
	amplitudeTreeLeaf&
	amplitudeTreeLeaf::operator *=(const T& factor)
	{
		for (unsigned int i = 0; i < nmbIncohSubAmps(); ++i)
			_incohSubAmps[i] *= factor;
		return *this;
	}


	template<typename T>
	inline
	amplitudeTreeLeaf&
	amplitudeTreeLeaf::operator /=(const T& factor)
	{
		for (unsigned int i = 0; i < nmbIncohSubAmps(); ++i)
			_incohSubAmps[i] /= factor;
		return *this;
	}


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

	template<typename T>
	inline
	amplitudeTreeLeaf
	operator *(const amplitudeTreeLeaf& amp,
	           const T&                 factor)
	{
		amplitudeTreeLeaf result = amp;
		result *= factor;
		return result;
	}

	template<typename T>
	inline
	amplitudeTreeLeaf
	operator *(const T&                 factor,
	           const amplitudeTreeLeaf& amp)
	{
		return amp * factor;
	}

	template<typename T>
	inline
	amplitudeTreeLeaf
	operator /(const amplitudeTreeLeaf& amp,
	           const T&                 factor)
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
