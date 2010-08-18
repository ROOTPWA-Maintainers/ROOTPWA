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


#include "cudaLikelihoodInterface.cu"


using namespace rpwa;
using namespace std;


template<typename complexT, typename T>
T
SumArray(complexT*          decayAmps,
         complexT*          prodAmps,
         const T            prodAmpFlat,
         const unsigned int nmbEvents,
         const unsigned int rank,
         const unsigned int nmbWavesRefl[2])
{
	const T            prodAmpFlat2   = prodAmpFlat * prodAmpFlat;
	const unsigned int prodAmpDim [3] = {rank,      2, max(nmbWavesRefl[0], nmbWavesRefl[1])};
	const unsigned int decayAmpDim[3] = {2, max(nmbWavesRefl[0], nmbWavesRefl[1]), nmbEvents};
	// loop over events and calculate first term of log likelihood
	T logLikelihood = 0;
	for (unsigned int iEvt = 0; iEvt < nmbEvents; iEvt++) {
		T l = 0;  // likelihood for this event
		for (unsigned int iRank = 0; iRank < rank; ++iRank) {  // incoherent sum over ranks
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
				complexT ampProdSum = 0;  // amplitude sum for negative/positive reflectivity for this rank
				for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
					// compute likelihood term
					const unsigned int prodAmpIndices [3] = {iRank, iRefl, iWave};
					const unsigned int decayAmpIndices[3] = {iRefl, iWave, iEvt};
					ampProdSum = ampProdSum + prodAmps   [indicesToOffset<unsigned int>(prodAmpIndices,  prodAmpDim,  3)]
						* decayAmps[indicesToOffset<unsigned int>(decayAmpIndices, decayAmpDim, 3)];
				}
				l += norm(ampProdSum);
			}
			//     assert(l >= 0);
		}  // end loop over rank
		l             += prodAmpFlat2;
		logLikelihood -= log(l);  // accumulate log likelihood
	}  // end loop over events
  
	return logLikelihood;
}


void
StartCalc(const unsigned int nmbEvents,
          const unsigned int rank,
          const unsigned int nmbWavesRefl[2],
          const unsigned int iterat)
{
	// set decay amplitudes
	const unsigned int     decayAmpDim[3] = {2, max(nmbWavesRefl[0], nmbWavesRefl[1]), nmbEvents};
	rpwa::complex<Scalar>* decayAmps;
	const unsigned int     decayAmpsSize
		= allocatePseudoNdimArray<rpwa::complex<double>, unsigned int>(decayAmps, decayAmpDim, 3);
	//  printInfo << "size of decay amplitude array is " << decayAmpsSize / (1024. * 1024.) << " MiBytes; "
	//	    << 100 * nmbEvents * (nmbWavesRefl[0] + nmbWavesRefl[1]) * sizeof(rpwa::complex<double>)
	//               / (double)decayAmpsSize << " % used" << endl;
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave)
			for (unsigned int iEvt = 0; iEvt < nmbEvents; ++iEvt) {
				const unsigned int decayAmpIndices[3] = {iRefl, iWave, iEvt};
				const unsigned int decayAmpOffset     = indicesToOffset<unsigned int>(decayAmpIndices,
					                                                                    decayAmpDim, 3);
				const double       val                = iEvt * 1000 + iRefl * 100 + iWave;
				decayAmps[decayAmpOffset] = rpwa::complex<double>(val, val + 0.5);
			}

	// set production amplitudes
	const unsigned int     prodAmpDim[3] = {rank, 2, max(nmbWavesRefl[0], nmbWavesRefl[1])};
	rpwa::complex<double>* prodAmps;
	const unsigned int     prodAmpsSize  = allocatePseudoNdimArray<rpwa::complex<double>, unsigned int>
		(prodAmps, prodAmpDim, 3);
	// printInfo << "size of production amplitude array is " << prodAmpsSize << " bytes" << endl;
	unsigned int parIndex = 1;
	for (unsigned int iRank = 0; iRank < rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < nmbWavesRefl[iRefl]; ++iWave) {
				const unsigned int prodAmpIndices[3] = {iRank, iRefl, iWave};
				const unsigned int prodAmpOffset     = indicesToOffset<unsigned int>(prodAmpIndices,
				                                                                     prodAmpDim, 3);
				const double       val               = iRank * 1000 + iRefl * 100 + iWave;
				prodAmps[prodAmpOffset] = rpwa::complex<double>(val, val + 0.5);
			}
  
	const double prodAmpFlat = parIndex;

	Scalar resultcpu = 0;
	Scalar resultcuda = 0;
	//Scalar resultcuda2 = 0;

	unsigned int timer = 0;
	cutCreateTimer( &timer);
	cutStartTimer( timer);

	for (int i=0; i < iterat; i++)
	{
	  resultcpu = SumArray<rpwa::complex<Scalar>, Scalar>(decayAmps, prodAmps, prodAmpFlat,
	                                                      nmbEvents, rank, nmbWavesRefl);
	}

	cout << "CPU result:" << resultcpu << endl;

	cutStopTimer( timer);
	printf( "time: %f (ms)\n", cutGetTimerValue( timer));
	cutDeleteTimer( timer);
	cout << endl;

	// unsigned int timer2 = 0;
	// cutCreateTimer( &timer2);
	// cutStartTimer( timer2);

	// for (int i=0; i < iterat; i++)
	// {
	//   resultcuda2 = SumArrayCUDATest<rpwa::complex<Scalar>, Scalar>(decayAmps,decayAmpsSize,prodAmps,prodAmpsSize,prodAmpFlat,nmbEvents,rank,nmbWavesRefl);
	// }

	// cout << "GPU Test result:" << resultcuda2 << endl;

	// cutStopTimer( timer2);
	// printf( "time: %f (ms)\n", cutGetTimerValue( timer2));
	// cutDeleteTimer( timer2);
	// cout << endl;

	unsigned int timer3 = 0;
	cutCreateTimer( &timer3);
	cutStartTimer( timer3);

	unsigned int nmbThreadsPerBlock = 0;
	unsigned int nmbBlocks          = 0;
	rpwa::complex<Scalar>* d_decayAmps;

	PrepareCUDA<rpwa::complex<Scalar> >(decayAmps, decayAmpsSize, d_decayAmps, nmbBlocks, nmbThreadsPerBlock);

	printf("nmbThreadsPerBlock: %i \n", nmbThreadsPerBlock);
	printf("nmbBlocks         : %i \n", nmbBlocks);
  
	for (int i=0; i < iterat; i++)
	{
		resultcuda = runCudaLogLikelihoodKernels<rpwa::complex<Scalar> >
		  (prodAmps, prodAmpsSize, prodAmpFlat, d_decayAmps, nmbEvents,
		   rank, nmbWavesRefl, nmbBlocks, nmbThreadsPerBlock);
	}

	cout << "GPU result:" << resultcuda << endl;

	cutStopTimer( timer3);
	printf( "time: %f (ms)\n", cutGetTimerValue( timer3));
	cutDeleteTimer( timer3);
	cout << endl;

	if (abs((resultcuda - resultcpu) / resultcpu) <= 10^(-10)) {printf("PASSED \n");} else {printf("FAILED \n");}

}


int
main(int    argc,
     char** argv) 
{
	const unsigned int nmbEvents = 100000;
	// const unsigned int nmbEvents = 32 * 4 * 10;
	const unsigned int rank = 2;
	const unsigned int nmbWavesRefl[2] = {7,34};
	const unsigned int iterat = 10;
	// const unsigned int rank = 1;
	// const unsigned int nmbWavesRefl[2] = {1,1};
	//  printf("Num Threads: %i \n", nmbThreadsPerBlock);
	//  printf("Num Blocks: %i \n", nmbBlocks);
  
	StartCalc(nmbEvents, rank, nmbWavesRefl, iterat); // change the datatype only here

	//  cutilExit(argc, argv);
}
