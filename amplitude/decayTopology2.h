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
//      amplitude calculation
//      internally the decay process is represented as a graph
//      the graph is constraint to contain exactly one production
//      vertex and at least one interaction veretx; in addtion for
//      each final state particle a corresponding final state vertex
//      will be created for internal use
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef DECAYTOPOLOGY2_H
#define DECAYTOPOLOGY2_H


// #include <vector>
// #include <map>

#include "particle.h"
#include "interactionVertex.h"
#include "fsVertex.h"
#include "decayGraph.hpp"


namespace rpwa {


  typedef decayGraph<interactionVertex, particle> decayGraphType;


  class decayTopology2 : public decayGraphType {
	
  public:
			
    decayTopology2();
    decayTopology2(const VPtr&                         productionVertex,
     		   const std::vector<VPtr>&            interactionVertices,
     		   const std::vector<rpwa::particle*>& fsParticles);
    virtual ~decayTopology2();



    // decayTopology2(const decayTopology2&                   topo);
  //   virtual decayTopology2& operator = (const decayTopology2& topo);
    
  //   virtual unsigned int nmbFsParticles() const { return _fsParticles.size(); }  ///< returne number of final state particles
  //   virtual unsigned int nmbVertices()    const { return _vertices.size();    }  ///< returns number of interaction vertices
  //   virtual std::vector<particle*>& fsParticles() { return _fsParticles; }  ///< returns final state particles
  //   virtual std::vector<interactionVertex*>& interactionVertices() { return _vertices; }  ///< returns all interaction vertices (excluding production vertex) ordered by depth-first; first vertex is X-decay vertex
  //   virtual interactionVertex& productionVertex() { return *_productionVertex; }  ///< returns production vertex
  //   virtual interactionVertex& xInteractionVertex() { return *_vertices[0]; }  ///< returns X-decay vertex

  //   virtual bool dataAreValid()   const { return false; }  ///< indicates whether data are complete and valid
  //   virtual bool verifyTopology() const { return true;  }  ///< returns whether decay has the correct topology
  //   virtual bool checkConsistency() const { return true; } ///< checks quantum number conservation rules on all vertices

  //   virtual std::ostream& print(std::ostream& out) const;  ///< prints decay topology in human-readable form
  //   virtual std::ostream& writeGraphViz(std::ostream& out) const;  ///< writes graph in GraphViz DOT format

    virtual void clear();  ///< deletes all information

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag

    
  //   //!!! bad quick hacks
  //   decayGraph subGraph(interactionVertex& startVertex);
  //   interactionVertex* vertex(particle& part);


  // protected:


  //   static decayGraph deepCopyGraph(const decayGraph& srcGraph,
  // 				    const bool        copyFsParticles = false);


  private:
    
    virtual decayTopology2& constructDecay(const VPtr&                         productionVertex,
					   const std::vector<VPtr>&            interactionVertices,
					   const std::vector<rpwa::particle*>& fsParticles);  ///< constructs the decay graph based on final state particles and vertices


    VPtr                                   _prodVertex;           ///< pointer to production vertex
    std::vector<VPtr>                      _intVertices;          ///< array of interaction vertices excluding production vertex; ordered depth-first
    std::map<rpwa::particle*, fsVertexPtr> _fsParticleVertexMap;  ///< maps final state particles to final state vertices
    
    static bool _debug;  ///< if set to true, debug messages are printed
    
  };
  

  // inline
  // std::ostream&
  // operator << (std::ostream&        out,
  // 	       const decayTopology2& topo) { return topo.print(out); }


} // namespace rpwa


#endif  // DECAYTOPOLOGY2_H
