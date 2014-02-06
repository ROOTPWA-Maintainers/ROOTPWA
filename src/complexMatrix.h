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

#include "TObject.h"
#include "TComplex.h"
#include "TMatrixD.h"

#include "reportingUtilsRoot.hpp"

namespace rpwa {

	class complexMatrix : public TObject {

	public:

		complexMatrix() { }

		complexMatrix(TMatrixD re, TMatrixD im) :
				_re(re), _im(im) { }

		complexMatrix(const int i, const int j) :
				_re(i, j), _im(i, j) { }

		~complexMatrix() { }

		void ResizeTo(const int i, const int j) { _re.ResizeTo(i, j); _im.ResizeTo(i, j); }
		void set(const int i, const int j, const std::complex<double>& c);
		void set(const int i, const int j, const TComplex& c);
		TComplex get(const int i, const int j) const { return TComplex(_re[i][j], _im[i][j]); }
		TComplex operator()(const int i, const int j) const { return this->get(i, j); }
		int nRows() const { return _re.GetNrows(); }
		int nCols() const { return _re.GetNcols(); }
		virtual void Print(const Option_t* = "") const { _re.Print(); _im.Print(); }

		complexMatrix t() const; // return transpose matrix
		complexMatrix dagger() const; // return adjoint matrix

		friend complexMatrix operator*(const complexMatrix& c1, const complexMatrix& c2);
		friend complexMatrix operator-(const complexMatrix& c1, const complexMatrix& c2);
		friend complexMatrix operator+(const complexMatrix& c1, const complexMatrix& c2);

	private:

		TMatrixD _re;
		TMatrixD _im;

	public:

		ClassDef(complexMatrix, 1);

	};

	complexMatrix operator*(const complexMatrix& c1, const complexMatrix& c2);
	complexMatrix operator-(const complexMatrix& c1, const complexMatrix& c2);
	complexMatrix operator+(const complexMatrix& c1, const complexMatrix& c2);

	inline std::ostream&
	operator <<(std::ostream& out, const complexMatrix& A) {
		for(int row = 0; row < A.nRows(); ++row) {
			out << "row " << row << " = (";
			for(int col = 0; col < A.nCols(); ++col) {
				out << A(row, col);
				if(col < A.nCols() - 1) {
					out << ", ";
				}
			}
			if(row < A.nRows() - 1) {
				out << "), " << std::endl;
			} else {
				out << ")";
			}
		}
		return out;
	}

}

#endif  // TCMATRIX_HH
