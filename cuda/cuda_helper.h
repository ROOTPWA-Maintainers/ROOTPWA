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
//      Test CUDA program
//
//
// Author List:
//      Philipp Meyer          TUM            (original author)
//
//
//-------------------------------------------------------------------------
#ifndef CUDA_HELPER
#define CUDA_HELPER

typedef double Scalar;
#define ALIGN 2 * sizeof( Scalar)

// my includes:
#include "CudaComplex.h" // includes my Complex datatype

////////////////////////////////////////////////////////////////////////////////
// declarations:

double SumArrayCUDA2(rpwa::complex<double>* prodAmps, const unsigned int prod_mem_size, const double prodAmpFlat, const unsigned int nmbEvents, const unsigned int rank, const unsigned int nmbWavesRefl[2], rpwa::complex<double>* d_decayAmps, unsigned int num_threads, unsigned int num_blocks);

void PrepareCUDA2(rpwa::complex<double>* decayAmps, const unsigned int decay_mem_size, rpwa::complex<double>** d_decayAmps,unsigned int &num_threads, unsigned int &num_blocks);

#endif // CUDA_HELPER
