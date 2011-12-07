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
//      base class that encapsulates naming scheme for partial waves
//
//      in the general case the intensity is parameterized by
//
//      I = sum_j sum_k sum_l |sum_i T_i^{j l} A_i^{j k}|^2
//
//      where
//      T_i^{j l} - transition amplitude
//      A_i^{j k} - decay amplitude
//      i         - coherent sum index common for T and A (e.g. I^G J^PC M [isobar decay chain])
//      j         - incoherent sum index common for T and A (e.g. reflectivity)
//      k         - incoherent sum index for A only (e.g. FS-particle helicity)
//      l         - incoherent sum index for T only (e.g. rank)
//
//      so both T and A are defined by 3 sets of quantum numbers
//      the total amplitude name thus consists of 3 substrings
//
//      the spin-density matrix is given by
//
//      rho^j_{i, i'} = sum_l T_i^{j l} T*_{i'}^{j l}
//
//      so partial waves are characterized by i and j only
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "waveName.h"

  
using namespace std;
using namespace boost;
using namespace rpwa;


ClassImp(waveName);


bool waveName::_debug = false;


waveName::waveName()
	: TObject      (),
	  _cohQnLabel  (""),
		_incohQnLabel("")
{
	//waveName::Class()->IgnoreTObjectStreamer();  // don't store TObject's fBits and fUniqueID
}


waveName::waveName(const isobarAmplitudePtr& amp)
{
	setQnLabel(amp);
}


waveName::~waveName()
{ }


void
waveName::clear()
{
	_cohQnLabel   = "";
	_incohQnLabel = "";
}


waveName&
waveName::operator =(const waveName& ampName)
{
	if (this != &ampName) {
		TObject::operator =(ampName);
		_cohQnLabel   = ampName._cohQnLabel;
		_incohQnLabel = ampName._incohQnLabel;
	}
	return *this;
}


void
waveName::setQnLabel(const isobarAmplitudePtr& amp)
{
	if (not amp) {
		printWarn << "null pointer for amplitude. cannot construct amplitude name." << endl;
		_cohQnLabel   = "";
		_incohQnLabel = "";
		return;
	} 
	const isobarDecayTopologyPtr& topo = amp->decayTopology();
	if (not topo) {
		printWarn << "null pointer for topology. cannot construct amplitude name." << endl;
		_cohQnLabel   = "";
		_incohQnLabel = "";
		return;
	} 
	if (not topo->checkTopology() or not topo->checkConsistency()) {
		printWarn << "decay topology has issues. cannot construct amplitude name." << endl;
		_cohQnLabel   = "";
		_incohQnLabel = "";
		return;
	}
	const particle& X = *(topo->XParticle());

	// X quantum numbers that are coherently summed for T and A
	ostringstream cohQn;
	cohQn << "I=" << spinQn(X.isospin());
	if (X.G() != 0)
		cohQn << ",G=" << parityQn(X.G());
	cohQn << ",J=" << spinQn(X.J())
	      << ",P=" << parityQn(X.P());
	if ((X.C() != 0))
		cohQn << ",C=" << parityQn(X.C());
	cohQn << ",M=" << spinQn(X.spinProj());

	// X decay chain
	const string chain = decayChain(*topo, *(topo->XIsobarDecayVertex()));
	cohQn << ",X" << chain << "";

	// set common coherent quantum numbers for T and A
	_cohQnLabel = cohQn.str();

	// set quantum numbers that are coherently summed for T and A
	if (amp->reflectivityBasis()) {
		ostringstream label;
		label << "refl=" << parityQn(X.reflectivity());
		_incohQnLabel = label.str();
	} else
		_incohQnLabel = "";

	//replace_all(_cohQnLabel,   "(", "_");
	//replace_all(_cohQnLabel,   ")", "_");
	//replace_all(_incohQnLabel, "(", "_");
	//replace_all(_incohQnLabel, ")", "_");
	return;
}


string
waveName::decayChain(const isobarDecayTopology& topo,
                     const isobarDecayVertex&   currentVertex)
{
	// recurse down decay chain
	ostringstream chain;
	// first daughter
	chain << "={" << currentVertex.daughter1()->name();
	if (not topo.isFsParticle(currentVertex.daughter1()))
		chain << decayChain(topo,
		  *static_pointer_cast<isobarDecayVertex>(topo.toVertex(currentVertex.daughter1())));
	// L, S
	chain << "[" << spinQn(currentVertex.L()) << "," << spinQn(currentVertex.S()) << "]";
	// second daughter
	chain << currentVertex.daughter2()->name();
	if (not topo.isFsParticle(currentVertex.daughter2()))
		chain << decayChain(topo,
			*static_pointer_cast<isobarDecayVertex>(topo.toVertex(currentVertex.daughter2())));
	chain << "}";
	return chain.str();
}
