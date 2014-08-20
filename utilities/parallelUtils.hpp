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
//      functions for simple often used parallel computations
//
//
// Author List:
//      Armin Gensler          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef PARALLELUTILS_HPP
#define PARALLELUTILS_HPP

#include <complex>
#include <iostream>
#include <vector>
#include <boost/date_time/posix_time/posix_time.hpp>

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "reportingUtils.hpp"
#include "reportingUtilsRoot.hpp"

namespace rpwa {
	
	template<typename T>
	inline void parallelAdd(std::vector<T>& xs, const std::vector<T>& ys)
	{
		std::cout << "EPL: parallelAdd" << std::endl;
		if(xs.size() != ys.size()) {
			printErr << "size of per-event-data vectors does not match. aborting." << std::endl;
			throw;
		}
		// !! EVENT PARALLEL LOOP
		boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
		const unsigned int size = xs.size();
		#pragma omp parallel for
		for(unsigned int i = 0; i < size; ++i) {
			xs[i] += ys[i];
		}
		boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
		uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
		std::cout << "    timediff = " << timeDiff << std::endl;
	}
	
	template<typename T>
	inline void parallelSub(std::vector<T>& xs, const std::vector<T>& ys)
	{
		std::cout << "EPL: parallelSub" << std::endl;
		if(xs.size() != ys.size()) {
			printErr << "size of per-event-data vectors does not match. aborting." << std::endl;
			throw;
		}
		// !! EVENT PARALLEL LOOP
		boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
		const unsigned int size = xs.size();
		#pragma omp parallel for
		for(unsigned int i = 0; i < size; ++i) {
			xs[i] -= ys[i];
		}
		boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
		uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
		std::cout << "    timediff = " << timeDiff << std::endl;
	}
	
	template<typename T>
	inline void parallelMultiply(std::vector<T>& xs, const T& factor)
	{
		std::cout << "EPL: parallelMultiply" << std::endl;
		// !! EVENT PARALLEL LOOP
		boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
		const unsigned int size = xs.size();
		#pragma omp parallel for
		for(unsigned int i = 0; i < size; ++i) {
			xs[i] = xs[i] * factor;
		}
		boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
		uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
		std::cout << "    timediff = " << timeDiff << std::endl;
	}
	
	template<typename T>
	inline void parallelDivide(std::vector<T>& xs, const T& factor)
	{
		std::cout << "EPL: parallelDivide" << std::endl;
		// !! EVENT PARALLEL LOOP
		boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
		const unsigned int size = xs.size();
		#pragma omp parallel for
		for(unsigned int i = 0; i < size; ++i) {
			xs[i] = xs[i] / factor;
		}
		boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
		uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
		std::cout << "    timediff = " << timeDiff << std::endl;
	}

	inline void parallelLorentzRotation(std::vector<TLorentzVector>& lzVec, const std::vector<TLorentzRotation>& trafo)
	{
		std::cout << "EPL: parallelLorentzRotation" << std::endl;
		if(lzVec.size() != trafo.size()) {
			printErr << "size of per-event-data vectors does not match. aborting." << std::endl;
			throw;
		}
		// !! EVENT PARALLEL LOOP
		boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
		const unsigned int size = lzVec.size();
		#pragma omp parallel for
		for(unsigned int i = 0; i < size; ++i) {
			lzVec[i].Transform(trafo[i]);
		}
		boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
		uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
		std::cout << "    timediff = " << timeDiff << std::endl;
	}
	
	inline void parallelLorentzBoost(std::vector<TLorentzVector>& lzVec, const std::vector<TVector3>& trafo)
	{
		std::cout << "EPL: parallelLorentzBoost" << std::endl;
		if(lzVec.size() != trafo.size()) {
			printErr << "size of per-event-data vectors does not match. aborting." << std::endl;
			throw;
		}
		// !! EVENT PARALLEL LOOP
		boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
		const unsigned int size = lzVec.size();
		#pragma omp parallel for
		for(unsigned int i = 0; i < size; ++i) {
			lzVec[i].Boost(trafo[i]);
		}
		boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
		uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
		std::cout << "    timediff = " << timeDiff << std::endl;
	}
	
	/*inline void parallelLorentzVectorMag2(const std::vector<TLorentzVector>& lzVec, std::vector<double>& result)
	{
		std::cout << "EPL: parallelLorentzVectorMag2" << std::endl;
		if(lzVec.size() != result.size()) {
			printErr << "size of per-event-data vectors does not match. aborting." << std::endl;
			throw;
		}
		// !! EVENT PARALLEL LOOP
		boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
		const unsigned int size = lzVec.size();
		#pragma omp parallel for
		for(unsigned int i = 0; i < size; ++i) {
			result[i] = lzVec[i].Mag2();
		}
		boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
		uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
		std::cout << "    timediff = " << timeDiff << std::endl;
	}*/
	
	/*inline void parallelLorentzVectorScalarProduct(const std::vector<TLorentzVector>& xs, 
						       const std::vector<TLorentzVector>& ys, 
						       std::vector<double>& result)
	{
		std::cout << "EPL: parallelLorentzVectorScalarProduct" << std::endl;
		if(xs.size() != ys.size() || xs.size() != result.size()) {
			printErr << "size of per-event-data vectors does not match. aborting." << std::endl;
			throw;
		}
		// !! EVENT PARALLEL LOOP
		boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
		const unsigned int size = xs.size();
		#pragma omp parallel for
		for(unsigned int i = 0; i < size; ++i) {
			result[i] = xs[i] * ys[i];
		}
		boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
		uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
		std::cout << "    timediff = " << timeDiff << std::endl;
	}*/
	
	/*inline void parallelLorentzVectorNegate3(std::vector<TLorentzVector>& lzVec)
	{
		std::cout << "EPL: parallelLorentzVectorNegate3" << std::endl;
		// !! EVENT PARALLEL LOOP
		boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
		const unsigned int size = lzVec.size();
		#pragma omp parallel for
		for(unsigned int i = 0; i < size; ++i) {
			lzVec[i] = TLorentzVector(-lzVec[i].Vect(), lzVec[i].E());
		}
		boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
		uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
		std::cout << "    timediff = " << timeDiff << std::endl;
	}*/
	
	inline void parallelLorentzRotationInvert(std::vector<TLorentzRotation>& lzRot)
	{
		std::cout << "EPL: parallelLorentzRotationInvert" << std::endl;
		// !! EVENT PARALLEL LOOP
		boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
		const unsigned int size = lzRot.size();
		#pragma omp parallel for
		for(unsigned int i = 0; i < size; ++i) {
			lzRot[i].Invert();
		}
		boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
		uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
		std::cout << "    timediff = " << timeDiff << std::endl;
	}

}  // namespace rpwa


#endif  // PARALLELUTILS_HPP
