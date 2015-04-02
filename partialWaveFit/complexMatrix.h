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
// a matrix with complex elements

#ifndef COMPLEXMATRIX_HH
#define COMPLEXMATRIX_HH

#include <complex>

#ifndef __CINT__
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#endif

#include "TObject.h"


namespace rpwa {

	class complexMatrix : public TObject {

	public:

		// BEWARE: the matrix constructor does not set the matrix elements
		//         to zero, neither does resize. So if you initialize the
		//         complexMatrix, make sure to set all the entries to zero.

		complexMatrix()
			: _matrix(),
			  _size1(0),
			  _size2(0),
			  _nmbDataElements(0),
			  _data(0) { }

		complexMatrix(const int i, const int j)
			: _matrix(i, j),
			  _size1(0),
			  _size2(0),
			  _nmbDataElements(0),
			  _data(0) { }

		~complexMatrix();

		void resizeTo(const int i, const int j) { _matrix.resize(i, j); }
		void set(const int i, const int j, const std::complex<double>& c) { _matrix(i, j) = c; }
		std::complex<double> get(const int i, const int j) const { return _matrix(i, j); }
		std::complex<double> operator()(const int i, const int j) const { return this->get(i, j); }
		unsigned int nRows() const { return _matrix.size1(); }
		unsigned int nCols() const { return _matrix.size2(); }

		bool operator==(const complexMatrix& rhs) const;
		bool equalToPrecision(const complexMatrix& rhs, const double& precision = 1e-16, const bool& verbose = false) const;

		virtual void Print(const Option_t* = "") const { std::cout << _matrix << std::endl; }

		complexMatrix t() const; // return transpose matrix
		complexMatrix dagger() const; // return adjoint matrix
		std::complex<double> determinant() const;

		friend complexMatrix operator*(const complexMatrix& c1, const complexMatrix& c2);
		friend complexMatrix operator-(const complexMatrix& c1, const complexMatrix& c2);
		friend complexMatrix operator+(const complexMatrix& c1, const complexMatrix& c2);
		friend std::ostream& operator<<(std::ostream& out, const complexMatrix& A);

	private:

#ifndef __CINT__
		complexMatrix(boost::numeric::ublas::matrix<std::complex<double> > matrix)
			: _matrix(matrix),
			  _size1(0),
			  _size2(0),
			  _nmbDataElements(0),
			  _data(0) { }

		boost::numeric::ublas::matrix<std::complex<double> > _matrix; //!
#endif

	public:

		void readMatrix();
		void storeMatrix();

	private:

		unsigned int _size1;
		unsigned int _size2;
		unsigned int _nmbDataElements;
		std::complex<double>* _data; //[_nmbDataElements]

	public:

		ClassDef(complexMatrix, 2);

	};

	complexMatrix operator*(const complexMatrix& c1, const complexMatrix& c2);
	complexMatrix operator-(const complexMatrix& c1, const complexMatrix& c2);
	complexMatrix operator+(const complexMatrix& c1, const complexMatrix& c2);
	std::ostream& operator<<(std::ostream& out, const complexMatrix& A);

}

#endif  // COMPLEXMATRIX_HH
