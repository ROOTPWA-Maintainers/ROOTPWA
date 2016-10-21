///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2016 Sebastian Uhl (TUM)
//
//    This file is part of ROOTPWA
//
//    ROOTPWA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ROOTPWA is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ROOTPWA.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      implementation of the storage for the data required by the
//      resonance fit
//
//-------------------------------------------------------------------------


#include "data.h"

#include <reportingUtils.hpp>


namespace {


	template<typename T>
	void
	checkSize(const std::vector<T>& obj,
	          const size_t dim1,
	          const std::string& msg1)
	{
		if(dim1 != obj.size()) {
			printErr << msg1 << " expected: " << dim1 << ", found: " << obj.size() << ". Aborting..." << std::endl;
			throw;
		}
	}


	template<typename T, const size_t dim = 0>
	void
	checkSize(const boost::multi_array<T, 1+dim>& obj,
	          const size_t dim1,
	          const std::string& msg1)
	{
		if(dim1 != *(obj.shape()+dim)) {
			printErr << msg1 << " expected: " << dim1 << ", found: " << *(obj.shape()+dim) << ". Aborting..." << std::endl;
			throw;
		}
	}


	template<typename T, const size_t dim = 0>
	void
	checkSize(const boost::multi_array<T, 2+dim>& obj,
	          const size_t dim1,
	          const std::string& msg1,
	          const size_t dim2,
	          const std::string& msg2)
	{
		if(dim1 != *(obj.shape()+dim)) {
			printErr << msg1 << " expected: " << dim1 << ", found: " << *(obj.shape()+dim) << ". Aborting..." << std::endl;
			throw;
		}
		checkSize<T, dim+1>(obj, dim2, msg2);
	}


}


rpwa::resonanceFit::data::data(const std::vector<size_t>& nrMassBins,
                               const boost::multi_array<double, 2>& massBinCenters)
	: _nrMassBins(nrMassBins),
	  _massBinCenters(massBinCenters)
{
	// get dimensions from one array and make sure that all other arrays
	// have the same dimensions
	const size_t nrBins = _nrMassBins.size();
	if(nrBins == 0) {
		printErr << "number of bins is zero, cannot perform the fit. Aborting..." << std::endl;
		throw;
	}

	const size_t maxMassBins = *(std::max_element(_nrMassBins.begin(), _nrMassBins.end()));
	if(maxMassBins == 0) {
		printErr << "maximal number of mass bins is zero, cannot perform the fit. Aborting..." << std::endl;
		throw;
	}

	checkSize(_nrMassBins,
	          nrBins, "number of bins is not correct for number of mass bins.");
	checkSize(_massBinCenters,
	          nrBins, "number of bins is not correct for centers of mass bins.",
	          maxMassBins, "maximal number of mass bins is not correct for centers of mass bins.");
}
