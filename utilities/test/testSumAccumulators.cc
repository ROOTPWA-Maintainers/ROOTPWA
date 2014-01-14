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
//      simple test for sum accumulators
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "reportingUtils.hpp"
#include "sumAccumulators.hpp"


using namespace std;
using namespace rpwa;
using namespace boost::accumulators;


int
main(int argc,
     char** argv)
{
	const size_t nmbValues = 1000000;
	const float  value     = 1e-6;
	accumulator_set<float, stats<tag::sum> >              naiveSumAcc;
	accumulator_set<float, stats<tag::sum(cascaded)> >    cascadedSumAcc
		(tag::cascadedSum::nmbElements = nmbValues);
	accumulator_set<float, stats<tag::sum(compensated)> > compensatedSumAcc;
	for (size_t i = 0; i < nmbValues; ++i) {
		naiveSumAcc      (value);
		cascadedSumAcc   (value);
		compensatedSumAcc(value);
	}
	printInfo << nmbValues << " * " << value << " = " << nmbValues * value << endl
	          << "    naive sum ......... " << sum(naiveSumAcc)       << endl
	          << "    cascaded sum ...... " << sum(cascadedSumAcc)    << endl
	          << "    compensated sum ... " << sum(compensatedSumAcc) << endl;
}
