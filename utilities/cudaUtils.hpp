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
//      some CUDA macros for code that is intended to run on host and device
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef CUDAUTILS_HPP
#define CUDAUTILS_HPP


#ifdef __CUDACC__
#define HOST __host__
#else
#define HOST
#endif


#ifdef __CUDACC__
#define DEVICE __device__
#else
#define DEVICE
#endif


#ifdef __CUDACC__
//#warning "info: generating CUDA device functions for nDimArray.hpp"
#define HOST_DEVICE __host__ __device__
#else
#define HOST_DEVICE
#endif


#ifdef __CUDACC__
#define ALIGN(nmbBytes) __align__(nmbBytes)
#else
#define ALIGN(nmbBytes)
#endif


#endif  // CUDAUTILS_HPP
