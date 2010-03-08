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
//      basic test program for vertex and decay topology
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "TVector3.h"

#include "utilities.h"
#include "particleDataTable.h"
#include "particle.h"
#include "diffractiveDissVertex.h"
//#include "isobarDecayVertex.h"


using namespace std;
using namespace rpwa;


int
main(int argc, char** argv)
{

  if (1) {
    // switch on debug output
    //particleProperties::setDebug(true);
    // particleDataTable::setDebug(true);
    particle::setDebug(true);
    vertex::setDebug(true);
  }

  particleDataTable& pdt = particleDataTable::instance();
  pdt.readFile();

  // test construction of vertices
  if (1) {
    TVector3 mom;
    mom = TVector3(1, 2, 3);
    particle beam("pi", -1, mom);
    mom = TVector3(2, 3, 4);
    particle X("X", -1, mom);
    printInfo << "created particles: " << endl
	      << beam << endl
	      << X    << endl;
    diffractiveDissVertex vert(beam, X);
    printInfo << "created vertex: " << endl
	      << vert;
  }

}
