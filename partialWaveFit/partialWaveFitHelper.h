///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2015 Sebastian Uhl (TUM)
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
//    along with ROOTPWA. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      helper functions related to the partial-wave fit
//
//-------------------------------------------------------------------------


#ifndef PARTIALWAVEFITHELPER_HH
#define PARTIALWAVEFITHELPER_HH


#include <iostream>


template<typename T> class TMatrixT;
template<typename T> class TVectorT;


namespace rpwa {


	class fitResult;
	template<typename T> class pwaLikelihood;


	namespace partialWaveFitHelper {


		void extractWaveList(const rpwa::fitResult& fitResult, std::ostream& waveList);

		template<typename T1, typename T2>
		void getEigenvectors(const rpwa::pwaLikelihood<T1>& L,
                                     const TMatrixT<T2>&      matrix,
                                     TMatrixT<T2>&            eigenvectors,
                                     TVectorT<T2>&            eigenvalues);

		int getReflectivity(const std::string& name);


	}  // namespace partialWaveFitHelper


}  // namespace rpwa


// Get the Eigenvectors of a matrix taking care of fixed parameters. A reduced
// matrix is created from the original one skipping the entries corresponding
// to fixed parameters. The Eigenvectors and Eigenvalues of this reduced matrix
// are calculated. Finally the matrix containing the Eigenvectors is modified
// by adding rows containing only 0. for the fixed parameters, the number of
// Eigenvalues does not change, so neither the number of columns, nor the size
// of the Eigenvalues vector is changed.
//
// The input matrix can be the Hessian or covariance matrix.
//
// Eigenvectors is a (nmbPars x nmbPars-nmbParsFixed) matrix, Eigenvalues is
// a vector with size (nmbPars-nmbParsFixed).
template<typename T1, typename T2>
inline
void
rpwa::partialWaveFitHelper::getEigenvectors(const pwaLikelihood<T1>& L,
                                            const TMatrixT<T2>&      matrix,
                                            TMatrixT<T2>&            eigenvectors,
                                            TVectorT<T2>&            eigenvalues)
{
	const unsigned int nmbPars      = L.nmbPars();
	const unsigned int nmbParsFixed = L.nmbParsFixed();

	TMatrixT<T2> reducedMatrix(nmbPars - nmbParsFixed, nmbPars - nmbParsFixed);
	{
		unsigned int iSkip = 0;
		for (unsigned int i = 0; i < nmbPars; ++i) {
			if (L.parFixed(i)) {
				iSkip++;
				continue;
			}

			unsigned int jSkip = 0;
			for (unsigned int j = 0; j < nmbPars; ++j) {
				if (L.parFixed(j)) {
					jSkip++;
					continue;
				}

				reducedMatrix[i - iSkip][j - jSkip] = matrix[i][j];
			}
		}
	}

	eigenvalues.ResizeTo(nmbPars - nmbParsFixed);
	TMatrixT<T2> reducedEigenvectors = reducedMatrix.EigenVectors(eigenvalues);

	eigenvectors.ResizeTo(nmbPars, nmbPars - nmbParsFixed);
	{
		unsigned int iSkip = 0;
		for (unsigned int i = 0; i < nmbPars; ++i) {
			if (L.parFixed(i)) {
				for (unsigned int j = 0; j < nmbPars - nmbParsFixed; ++j) {
					eigenvectors[i][j] = 0.;
				}
				iSkip++;
				continue;
			}

			for (unsigned int j = 0; j < nmbPars - nmbParsFixed; ++j) {
				eigenvectors[i][j] = reducedEigenvectors[i - iSkip][j];
			}
		}
	}
}


#endif  // PARTIALWAVEFITHELPER_HH
