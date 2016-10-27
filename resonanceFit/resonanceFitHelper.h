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


#include <boost/multi_array.hpp>

#include <TMatrixT.h>

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


		template<typename T>
		void
		copyMultiArrayElement(const T& source,
		                      T& target)
		{
			target = source;
		}


		inline
		void
		copyMultiArrayElement(const TMatrixT<double>& source,
		                      TMatrixT<double>& target)
		{
			target.ResizeTo(source);
			target = source;
		}


		template<typename T, size_t dim1, size_t dim2>
		void
		copyMultiArray(const boost::multi_array<T, dim1>& source,
		               std::vector<size_t> sourceIndices,
		               boost::multi_array<T, dim2>& target,
		               std::vector<size_t> targetIndices)
		{
			assert(sourceIndices.size() <= dim1);
			assert(targetIndices.size() <= dim2);
			assert(dim1-sourceIndices.size() == dim2-targetIndices.size());
			assert(((sourceIndices.size() == dim1) and (targetIndices.size() == dim2))
			       or ((sourceIndices.size() != dim1) and (targetIndices.size() != dim2)));

			if(sourceIndices.size() == dim1) {
				copyMultiArrayElement(source(sourceIndices), target(targetIndices));
			} else {
				const size_t maxIdx = source.shape()[sourceIndices.size()];

				sourceIndices.push_back(0);
				targetIndices.push_back(0);
				for(size_t idx = 0; idx < maxIdx; ++idx) {
					sourceIndices.back() = idx;
					targetIndices.back() = idx;
					copyMultiArray(source, sourceIndices, target, targetIndices);
				}
			}
		}


		// increase the extensions of a boost::multi_array such that one part
		// of it can simply be reset
		template<typename T>
		void
		adjustSize(boost::multi_array<T, 1>& master,
		           const size_t minSizeFirst = 1)
		{
			// initialize the next size with the current size
			std::vector<size_t> newSize(master.shape(), master.shape()+master.num_dimensions());

			// resize if the minimal size required for the first dimension
			// is larger than its current size
			if(newSize[0] < minSizeFirst) {
				newSize[0] = minSizeFirst;
				master.resize(newSize);
			}
		}


		// increase the extensions of a boost::multi_array such that one part
		// of it can simply be reset
		template<typename T, size_t dim>
		void
		adjustSize(boost::multi_array<T, dim>& master,
		           const boost::multi_array<T, dim-1>& part,
		           const size_t minSizeFirst = 1)
		{
			// initialize the next size with the current size
			std::vector<size_t> newSize(master.shape(), master.shape()+master.num_dimensions());

			// resize if the minimal size required for the first dimension
			// is larger than its current size
			bool resize = newSize[0] < minSizeFirst;
			newSize[0] = std::max(newSize[0], minSizeFirst);

			// compare the other dimensions with the dimenstions of the
			// part to be set
			for(size_t i = 1; i < dim; ++i) {
				if(newSize[i] < part.shape()[i-1]) {
					resize = true;
					newSize[i] = part.shape()[i-1];
				}
			}

			if(resize) {
				master.resize(newSize);
			}
		}


		// increase the extensions of a boost::multi_array such that one part
		// of it can simply be reset
		// special case for TMatrixT because TMatrixT::operator= does not
		// adjust the size of the matrices
		template<size_t dim>
		void
		adjustSize(boost::multi_array<TMatrixT<double>, dim>& master,
		           const boost::multi_array<TMatrixT<double>, dim-1>& part,
		           const size_t minSizeFirst = 1)
		{
			// initialize the next size with the current size
			std::vector<size_t> newSize(master.shape(), master.shape()+master.num_dimensions());

			// resize if the minimal size required for the first dimension
			// is larger than its current size
			bool resize = newSize[0] < minSizeFirst;
			newSize[0] = std::max(newSize[0], minSizeFirst);

			// compare the other dimensions with the dimenstions of the
			// part to be set
			for(size_t i = 1; i < dim; ++i) {
				if(newSize[i] < part.shape()[i-1]) {
					resize = true;
					newSize[i] = part.shape()[i-1];
				}
			}

			if(resize) {
				boost::multi_array<TMatrixT<double>, dim> temp(std::vector<size_t>(master.shape(), master.shape()+master.num_dimensions()));

				copyMultiArray(master, std::vector<size_t>(), temp, std::vector<size_t>());

				// clear the current content, in particular make sure
				// that after the next 'resize' all TMatrixT have
				// dimension 0x0
				master.resize(std::vector<size_t>(master.num_dimensions(), 0));

				// resize to new size
				master.resize(newSize);

				copyMultiArray(temp, std::vector<size_t>(), master, std::vector<size_t>());
			}
		}


		// increase the extensions of a std::vector such that a new element can
		// simply be set
		template<typename T>
		void
		adjustSize(std::vector<T>& master,
		           const size_t minSize)
		{
			if(master.size() < minSize) {
				master.resize(minSize);
			}
		}


		// increase the extensions of a boost::multi_array such that one part
		// of it can simply be reset, and also set this part
		template<typename T>
		void
		adjustSizeAndSet(boost::multi_array<T, 1>& master,
		                 const size_t idx,
		                 const T& part)
		{
			adjustSize(master, idx+1);
			master[idx] = part;
		}


		// increase the extensions of a boost::multi_array such that one part
		// of it can simply be reset, and also set this part
		template<typename T, size_t dim>
		void
		adjustSizeAndSet(boost::multi_array<T, dim>& master,
		                 const size_t idx,
		                 const boost::multi_array<T, dim-1>& part)
		{
			adjustSize(master, part, idx+1);
			copyMultiArray(part, std::vector<size_t>(), master, std::vector<size_t>(1, idx));
		}


		// increase the extensions of a std::vector such that a new element can
		// simply be set, and also set this part
		template<typename T>
		void
		adjustSizeAndSet(std::vector<T>& master,
		                 const size_t idx,
		                 const T& part)
		{
			adjustSize(master, idx+1);
			master[idx] = part;
		}


	}


}


#endif // RESONANCEFIT_RESONANCEFITHELPER_HH
