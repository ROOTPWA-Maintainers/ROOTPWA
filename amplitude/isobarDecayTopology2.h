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


#ifndef ISOBARDECAYTOPOLOGY2_H
#define ISOBARDECAYTOPOLOGY2_H


#include "isobarDecayVertex2.h"
#include "decayTopology2.h"


namespace rpwa {	

  class isobarDecayTopology2 : public decayTopology2 {
	
  public:
			
    isobarDecayTopology2();
    isobarDecayTopology2(const interactionVertexPtr&              productionVertex,
			 const std::vector<isobarDecayVertexPtr>& isobarDecayVertices,
			 const std::vector<particlePtr>&          fsParticles);
    isobarDecayTopology2(const interactionVertexPtr&              productionVertex,
			 const std::vector<interactionVertexPtr>& isobarDecayVertices,
			 const std::vector<particlePtr>&          fsParticles);
    isobarDecayTopology2(const isobarDecayTopology2&              topo);
    isobarDecayTopology2(const decayTopology2&                    topo);
    virtual ~isobarDecayTopology2();

    virtual isobarDecayTopology2& operator =(const isobarDecayTopology2& topo);
    virtual isobarDecayTopology2& operator =(const decayTopology2&       topo);
    virtual isobarDecayTopology2* clone(const bool cloneFsParticles      = false,
					const bool cloneProductionVertex = false) const;
    virtual void clear();  ///< deletes all information

    const std::vector<isobarDecayVertexPtr>& isobarDecayVertices() const { return _isobarVertices;    }  ///< returns all isobar decay vertices ordered by depth-first; first vertex is X-decay vertex
    const isobarDecayVertexPtr               xIsobarDecayVertex () const { return _isobarVertices[0]; }  ///< returns X-decay vertex

    bool checkTopology   () const;  ///< returns whether decay has the correct topology
    bool checkConsistency() const;  ///< checks conservation rules on all vertices

    isobarDecayTopology2 subDecay(const nodeDesc& startNd);  ///< returns sub-decay tree that starts at given vertex

    std::vector<isobarDecayTopology2> possibleDecays();  ///< constructs set of all possible decays given the final state particles and the constraints on I, J, L, and S
    
    const TLorentzVector& calcIsobarLzVec();  ///< (re)calculates Lorentz-vectors of all isobars in the decay from final state particles and returns Lorentz-vector of X-system

    virtual std::ostream& print(std::ostream& out) const;  ///< prints decay topology in human-readable form

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


  private:

    isobarDecayTopology2& constructDecay(const interactionVertexPtr&              productionVertex,
					 const std::vector<isobarDecayVertexPtr>& isobarDecayVertices,
					 const std::vector<particlePtr>&          fsParticles);  ///< constructs the decay graph based on final state particles and vertices
    isobarDecayTopology2& constructDecay(const interactionVertexPtr&              productionVertex,
					 const std::vector<interactionVertexPtr>& isobarDecayVertices,
					 const std::vector<particlePtr>&          fsParticles);  ///< constructs the decay graph based on final state particles and vertices

    void buildIsobarVertexArray();  ///< (re)builds array of isobar decay vertices

    decayTopologyGraphType joinDaughterGraphsX(const isobarDecayVertexPtr&   motherVertex,
					       const decayTopologyGraphType& daughterGraph1,
					       const decayTopologyGraphType& daughterGraph2,
					       nodeDesc&                     newTopNode);  ///< creates new graph with topNode as top node and graph1/2 as daughter graphs
    
    std::vector<isobarDecayVertexPtr> _isobarVertices;  ///< array of isobar-decay vertices excluding production vertex; ordered depth-first; this is a copy of the respective array in decayTopology

    static bool _debug;  ///< if set to true, debug messages are printed

  };
  

  inline
  std::ostream&
  operator <<(std::ostream&               out,
  	      const isobarDecayTopology2& topo)
  {
    return topo.print(out);
  }


} // namespace rpwa


#endif  // ISOBARDECAYTOPOLOGY2_H
