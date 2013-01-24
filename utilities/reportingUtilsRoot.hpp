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
//      some stream operators for common ROOT classes
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef REPORTINGUTILSROOT_HPP
#define REPORTINGUTILSROOT_HPP


#include <iostream>

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TComplex.h"
#include "TMatrixT.h"


namespace rpwa {


	//////////////////////////////////////////////////////////////////////////////
	// simple stream operators for some common ROOT classes
	inline
	std::ostream&
	operator << (std::ostream&   out,
	             const TVector3& vec)
	{
		out << "(" << vec.X() << ", " << vec.Y() << ", " << vec.Z() << ")";
		return out;
	}


	inline
	std::ostream&
	operator << (std::ostream&         out,
	             const TLorentzVector& vec)
	{
		out << "(" << vec.X() << ", " << vec.Y() << ", " << vec.Z() << "; " << vec.T() << ")";
		return out;
	}


	inline
	std::ostream&
	operator << (std::ostream&   out,
	             const TComplex& c)
	{
		out << "(" << c.Re() << ", " << c.Im() << ")";
		return out;
	}


	template<typename T>
	std::ostream&
	operator << (std::ostream&      out,
	             const TMatrixT<T>& A)
	{
		for (int row = 0; row < A.GetNrows(); ++row) {
			out << "row " << row << " = (";
			for (int col = 0; col < A.GetNcols(); ++col) {
				out << A[row][col];
				if (col < A.GetNcols() - 1)
					out << ", ";
			}
			if (row < A.GetNrows() - 1)
				out << "), " << std::endl;
			else
				out << ")";
		}
		return out;
	}


}  // namespace rpwa


#endif  // REPORTINGUTILSROOT_HPP
