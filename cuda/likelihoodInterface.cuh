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


#ifndef LIKELIHOODINTERFACE_CUH
#define LIKELIHOODINTERFACE_CUH


#include <iostream>


namespace rpwa {

	namespace cuda {


		// some helper classes to wrap kernel calls so that they can be
		// used as template parameters
		template<typename complexT, typename T>
		struct sumKernelCaller {
			typedef T value_type;
			static void call(const T*           d_sumsPrev,
			                 const unsigned int nmbSumsPrev,
			                 T*                 d_sumsNext,
			                 const unsigned int nmbSumsNext,
			                 const bool         debug = false);
		};


		template<typename complexT, typename T>
		struct sumDerivativesKernelCaller {
			typedef T value_type;
			static void call(const T*           d_sumsPrev,
			                 const unsigned int nmbSumsPrev,
			                 T*                 d_sumsNext,
			                 const unsigned int nmbSumsNext,
			                 const bool         debug = false);
		};


		template<typename complexT>
		class likelihoodInterface {

			template<typename U, typename T> friend class sumDerivativesKernelCaller;


		public:

			typedef typename complexT::value_type value_type;

			static likelihoodInterface& instance() { return _instance; }  ///< get singleton instance

			static bool init(const complexT*    decayAmps,
			                 const unsigned int nmbDecayAmps,
			                 const unsigned int nmbEvents,
			                 const unsigned int nmbWavesRefl[2],
			                 const bool         reshuffleArray = true);  ///< convenience routine that initializes the CUDA device and copies decay amplitudes into device memory
			static bool initCudaDevice();   ///< initializes CUDA device
			static void closeCudaDevice();  ///< closes CUDA device
			static bool loadDecayAmps(const complexT*    decayAmps,
			                          const unsigned int nmbDecayAmps,
			                          const unsigned int nmbEvents,
			                          const unsigned int nmbWavesRefl[2],
			                          const bool         reshuffleArray = true);  ///< copies decay amplitudes into CUDA device memory

			static bool         cudaInitialized   () { return _cudaInitialized;    }  ///< returns status of CUDA initialization
			static unsigned int totalDeviceMem    ();                                 ///< returns total memory capacity of used CUDA device
			static unsigned int availableDeviceMem();                                 ///< returns available memory capacity of used CUDA device
			static const struct cudaDeviceProp* deviceProperties();                   ///< returns pointer to properties of used CUDA device

			template<typename kernelT>
			static std::ostream& printKernelAttributes(std::ostream& out,
			                                           kernelT*      kernel);  ///< prints attributes of CUDA kernel
			template<typename kernelT>
			static void estimateOptimumKernelGrid(kernelT*           kernel,
			                                      unsigned int&      nmbBlocks,
			                                      unsigned int&      nmbThreadsPerBlock,
			                                      const unsigned int minNmbThreads = 0);

			static value_type logLikelihood(const complexT*    prodAmps,
			                                const unsigned int nmbProdAmps,
			                                const value_type   prodAmpFlat,
			                                const unsigned int rank);  ///< computes log likelihood for given production amplitudes
		
			static void logLikelihoodDeriv(const complexT*    prodAmps,
			                               const unsigned int nmbProdAmps,
			                               const value_type   prodAmpFlat,
			                               const unsigned int rank,
			                               complexT*          derivatives,
			                               value_type&        derivativeFlat);  ///< computes derivatives of log likelihood for given production amplitudes
		
			static std::ostream& print(std::ostream& out);  ///< prints properties of used CUDA device

			static bool debug() { return _debug; }                             ///< returns debug flag
			static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


		private:

			likelihoodInterface () { }
			~likelihoodInterface();
			likelihoodInterface (const likelihoodInterface&);
			likelihoodInterface& operator =(const likelihoodInterface&);

			template<class kernelCaller>
			static void cascadedKernelSum(const unsigned int                  nmbOfSumsAtEachStage,
			                              typename kernelCaller::value_type*& d_sumsPrev,
			                              unsigned int                        nmbSumsPrev,
			                              const unsigned int                  sumElementSize = 1);


			static likelihoodInterface _instance;  ///< singleton instance

			static bool                  _cudaInitialized;     ///< indicates whether CUDA environment was initialized correctly
			static int                   _nmbOfCudaDevices;    ///< number of found CUDA capable devices
			static int                   _cudaDeviceId;        ///< device ID of used CUDA device
			static struct cudaDeviceProp _cudaDeviceProp;      ///< properties of used CUDA device
			static complexT*             _d_decayAmps;         ///< device pointer to precalculated decay amplitudes
			static unsigned int          _nmbEvents;           ///< number of events to process
			static unsigned int          _nmbWavesRefl[2];     ///< number of waves for each reflectivity
			static unsigned int          _nmbWavesMax;         ///< maximum extent of wave index for production and decay amplitude arrays
			static unsigned int          _rank;                ///< rank of spin-density matrix

			static bool _debug;  ///< if set to true, debug messages are printed

		};


		template<typename complexT>
		inline
		std::ostream&
		operator <<(std::ostream&                        out,
		            const likelihoodInterface<complexT>& interface)
		{
			return interface.print(out);
		}


	}  // namespace cuda

}  // namespace rpwa


#endif // LIKELIHOODINTERFACE_CUH
