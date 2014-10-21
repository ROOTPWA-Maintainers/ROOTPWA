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
//      production vertex virtual base class
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef PRODUCTIONVERTEX_H
#define PRODUCTIONVERTEX_H


#include <complex>

#include <boost/shared_ptr.hpp>

#include "Typedefs.hpp"
#include "interactionVertex.h"


class TClonesArray;


namespace rpwa {


	class productionVertex;
	typedef boost::shared_ptr<productionVertex> productionVertexPtr;


	class productionVertex : public interactionVertex {

	public:

		productionVertex();
		virtual ~productionVertex();

		// production specific accessors
		virtual const ParVector<LorentzVector>& referenceLzVecs() const = 0;  ///< returns Lorentz vectors for a number of events that defines z-axis for angular distributions
		virtual const particlePtr&              XParticle      () const = 0;  ///< returns X particle

		virtual ParVector<Complex> productionAmps() const;  ///< returns production amplitudes all events stored in particles

		virtual void setXFlavorQN() = 0;  ///< general interface to set flavor quantum numbers of X (baryon nmb., S, C, B) based on production mechanism

		virtual bool initKinematicsData(const TClonesArray& names                      ) = 0;  ///< general interface to initialize input data format
		virtual bool readKinematicsData(const std::vector<ParVector<Vector3> >& momenta) = 0;  ///< general interface to read multiple input data events

		virtual bool revertMomenta() = 0;  ///< general interface to reset momenta to the values of last read event block

		virtual std::string name() const { return "productionVertex"; }  ///< returns label used in graph visualization, reporting, and key file

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	private:

		static bool _debug;  ///< if set to true, debug messages are printed

	};


}  // namespace rpwa


#endif  // PRODUCTIONVERTEX_H
