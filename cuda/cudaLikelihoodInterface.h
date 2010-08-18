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
//      interface for CUDA likelihood functions
//
//
// Author List:
//      Philipp Meyer          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef CUDALIKELIHOODINTERFACE_H
#define CUDALIKELIHOODINTERFACE_H


#include <iostream>

typedef double Scalar;
#define ALIGN 2 * sizeof(Scalar)
#include "CudaComplex.h"


namespace rpwa {


	template<typename complexT>
	class cudaLikelihoodInterface {

	public:

		typedef typename complexT::value_type value_type;

		static cudaLikelihoodInterface& instance() { return _instance; }  ///< get singleton instance

		static bool         cudaInitialized   () { return _cudaInitialized;    }  ///< returns status of CUDA initialization
		static unsigned int totalDeviceMem    ();                                 ///< returns total memory capacity of used CUDA device
		static unsigned int freeDeviceMem     ();                                 ///< returns available memory capacity of used CUDA device
		static unsigned int nmbBlocks         () { return _nmbBlocks;          }  ///< returns number of CUDA thread blocks 
		static unsigned int nmbThreadsPerBlock() { return _nmbThreadsPerBlock; }  ///< returns number of CUDA threads per block

		static bool init(const complexT*    decayAmps,
		                 const unsigned int nmbDecayAmps,
		                 const unsigned int nmbEvents,
		                 const unsigned int nmbWavesRefl[2],
		                 const bool         reshuffleArray = true);  ///< convenience routine that initializes the CUDA device and copies decay amplitudes into device memory
		static bool initCudaDevice();  ///< initializes CUDA device
		static bool loadDecayAmps(const complexT*    decayAmps,
		                          const unsigned int nmbDecayAmps,
		                          const unsigned int nmbEvents,
		                          const unsigned int nmbWavesRefl[2],
		                          const bool         reshuffleArray = true);  ///< copies decay amplitudes into CUDA device memory

		static value_type sumLogLikelihood(const complexT*    prodAmps,
		                                   const unsigned int nmbProdAmps,
		                                   const value_type   prodAmpFlat,
		                                   const unsigned int rank);  ///< computes log likelihood sum for given production amplitudes
		
		static std::ostream& print(std::ostream& out);  ///< prints properties of used CUDA device

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	private:

		cudaLikelihoodInterface () { }
		~cudaLikelihoodInterface();
		cudaLikelihoodInterface (const cudaLikelihoodInterface&);
		cudaLikelihoodInterface& operator =(const cudaLikelihoodInterface&);


		static cudaLikelihoodInterface _instance;  ///< singleton instance

		static bool                  _cudaInitialized;     ///< indicates whether CUDA environment was initialized correctly
		static int                   _nmbOfCudaDevices;    ///< number of found CUDA capable devices
		static int                   _cudaDeviceId;        ///< device ID of used CUDA device
		static struct cudaDeviceProp _cudaDeviceProp;      ///< properties of used CUDA device
		static unsigned int          _nmbBlocks;           ///< number of CUDA blocks used to run kernels
		static unsigned int          _nmbThreadsPerBlock;  ///< number of CUDA threads that is run per block
		static complexT*             _d_decayAmps;         ///< device pointer to precalculated decay amplitudes
		static unsigned int          _nmbEvents;           ///< number of events to process
		static unsigned int          _nmbWavesRefl[2];     ///< number of waves for each reflectivity

		static bool _debug;  ///< if set to true, debug messages are printed

	};


	template<typename complexT>
	inline
	std::ostream&
	operator <<(std::ostream&                            out,
	            const cudaLikelihoodInterface<complexT>& interface)
	{
		return interface.print(out);
	}


}  // namespace rpwa


#endif // CUDALIKELIHOODINTERFACE_H
