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
//      various helper functions used througout the complete resonance fit
//
//-------------------------------------------------------------------------


#ifndef RESONANCEFIT_RESONANCEFITHELPER_HH
#define RESONANCEFIT_RESONANCEFITHELPER_HH


#include <reportingUtils.hpp>


namespace rpwa {


	namespace resonanceFit {


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


		template<typename T, const size_t dim = 0>
		void
		checkSize(const boost::multi_array<T, 3+dim>& obj,
		          const size_t dim1,
		          const std::string& msg1,
		          const size_t dim2,
		          const std::string& msg2,
		          const size_t dim3,
		          const std::string& msg3)
		{
			if(dim1 != *(obj.shape()+dim)) {
				printErr << msg1 << " expected: " << dim1 << ", found: " << *(obj.shape()+dim) << ". Aborting..." << std::endl;
				throw;
			}
			checkSize<T, dim+1>(obj, dim2, msg2, dim3, msg3);
		}


		template<typename T, const size_t dim = 0>
		void
		checkSize(const boost::multi_array<T, 4+dim>& obj,
		          const size_t dim1,
		          const std::string& msg1,
		          const size_t dim2,
		          const std::string& msg2,
		          const size_t dim3,
		          const std::string& msg3,
		          const size_t dim4,
		          const std::string& msg4)
		{
			if(dim1 != *(obj.shape()+dim)) {
				printErr << msg1 << " expected: " << dim1 << ", found: " << *(obj.shape()+dim) << ". Aborting..." << std::endl;
				throw;
			}
			checkSize<T, dim+1>(obj, dim2, msg2, dim3, msg3, dim4, msg4);
		}


	}


}


#endif // RESONANCEFIT_RESONANCEFITHELPER_HH
