///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
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
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////

#include <boost/numeric/ublas/lu.hpp>

#include "complexMatrix.h"
#include "reportingUtils.hpp"

using namespace std;
using namespace rpwa;

namespace bnu = boost::numeric::ublas;

ClassImp(complexMatrix)


complexMatrix complexMatrix::t() const
{
	return complexMatrix(bnu::trans(_matrix));
}
complexMatrix complexMatrix::dagger() const
{
	return complexMatrix(bnu::conj(_matrix));
}


complex<double> complexMatrix::determinant() const
{
	// based on http://programmingexamples.net/wiki/CPP/Boost/Math/uBLAS/determinant
	if(_matrix.size1() != _matrix.size2()) {
		printWarn << "cannot calculate determinant for non-symmetric matrix, returning 0." << endl;
		return 0.;
	}
	bnu::matrix<complex<double> > matrixCopy(_matrix); // otherwise _matrix would be changed
	bnu::permutation_matrix<size_t> permutationMatrix(matrixCopy.size1());
	complex<double> det = 1.;
	if(bnu::lu_factorize(matrixCopy, permutationMatrix)) {
		det = 0.;
	} else {
		for(unsigned int i = 0; i < matrixCopy.size1(); ++i) {
			det *= matrixCopy(i, i); // multiply by elements on diagonal
		}
		int determinantSign = 1;
		size_t size = permutationMatrix.size();
		for(size_t i = 0; i < size; ++i) {
			if (i != permutationMatrix(i)) {
				determinantSign *= -1; // swap_rows would swap a pair of rows here, so we change sign
			}
		}
		det = det * (double)determinantSign;
	}
	return det;
}


complexMatrix rpwa::operator*(const complexMatrix& c1, const complexMatrix& c2)
{
	return complexMatrix(c1._matrix * c2._matrix);
}


complexMatrix rpwa::operator-(const complexMatrix& c1, const complexMatrix& c2)
{
	return complexMatrix(c1._matrix - c2._matrix);
}


complexMatrix rpwa::operator+(const complexMatrix& c1, const complexMatrix& c2)
{
	return complexMatrix(c1._matrix + c2._matrix);
}


ostream& rpwa::operator<<(ostream& out, const complexMatrix& A)
{
	out << A._matrix << endl;
	return out;
}


void complexMatrix::readMatrix()
{
	if(_size1 > 0 and _size2 > 0) {
		_matrix.resize(_size1, _size2);
		unsigned int dataIndex = 0;
		for(unsigned int i = 0; i < _size1; ++i) {
			for(unsigned int j = 0; j < _size2; ++j) {
				_matrix(i, j) = _data[dataIndex++];
			}
		}
	}
	delete [] _data;
	_data = 0;
}


void complexMatrix::storeMatrix()
{
	_size1 = _matrix.size1();
	_size2 = _matrix.size2();
	_nmbDataElements = _size1 * _size2;
	_data = new std::complex<double>[_nmbDataElements];
	unsigned int dataIndex = 0;
	for(unsigned int i = 0; i < _size1; ++i) {
		for(unsigned int j = 0; j < _size2; ++j) {
			_data[dataIndex++] = _matrix(i, j);
		}
	}
}


// copied from ampIntegralMatrix.cc
void
complexMatrix::Streamer(TBuffer& R__b)
{
	if (R__b.IsReading()) {
		R__b.ReadClassBuffer(rpwa::complexMatrix::Class(), this);
		readMatrix();
	} else {
		storeMatrix();
		R__b.WriteClassBuffer(rpwa::complexMatrix::Class(), this);
	}
}
