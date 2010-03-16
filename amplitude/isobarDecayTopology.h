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
//      container class that holds all external information for
//      amplitude calculation of isobar decay
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef ISOBARDECAYTOPOLOGY_H
#define ISOBARDECAYTOPOLOGY_H


#include "isobarDecayVertex.h"
#include "decayTopology.h"


namespace rpwa {	

  class isobarDecayTopology : public decayTopology {
	
  public:
			
    isobarDecayTopology();
    isobarDecayTopology(const isobarDecayTopology&             topo);
    isobarDecayTopology(const std::vector<particle*>&          fsParticles,
			interactionVertex&                     productionVertex,
			const std::vector<interactionVertex*>& interactionVertices);
    isobarDecayTopology(const std::vector<particle*>&          fsParticles,
			interactionVertex&                     productionVertex,
			const std::vector<isobarDecayVertex*>& isobarDecayVertices);
    virtual ~isobarDecayTopology();

    isobarDecayTopology& operator = (const isobarDecayTopology& topo);
    
    isobarDecayTopology& constructDecay(const std::vector<particle*>&          fsParticles,
					interactionVertex&                     productionVertex,
					const std::vector<interactionVertex*>& interactionVertices);  ///< constructs the decay graph based on final state particles and vertices
    isobarDecayTopology& constructDecay(const std::vector<particle*>&          fsParticles,
					interactionVertex&                     productionVertex,
					const std::vector<isobarDecayVertex*>& isobarDecayVvertices);  ///< constructs the decay graph based on final state particles and vertices

    std::vector<isobarDecayVertex*>& isobarDecayVertices() { return _vertices; }  ///< returns all isobar decay vertices ordered by depth-first; first vertex is X-decay vertex
    isobarDecayVertex& xIsobarDecayVertex() { return *_vertices[0]; }  ///< returns X-decay vertex

    bool verifyTopology() const;  ///< returns whether decay has the correct topology

    bool checkConsistency() const; ///< checks conservation laws at vertices

    
    const TLorentzVector& updateIsobarLzVec();  ///< (re)calculates Lorentz-vectors of all isobars in the decay from final state particles and returns Lorentz-vector of X-system

    void clear();  ///< deletes all information

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


  private:
    
    std::vector<isobarDecayVertex*> _vertices;  ///< array of isobar-decay vertices excluding production vertex; ordered depth-first; this is a copy of the respective array in decayTopology

    static bool _debug;  ///< if set to true, debug messages are printed

  };


} // namespace rpwa


#endif  // ISOBARDECAYTOPOLOGY_H
