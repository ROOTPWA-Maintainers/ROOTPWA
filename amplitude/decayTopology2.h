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
//      vertex and at least one interaction vertex; in addtion for
//      each final state particle a corresponding final state vertex
//      is created for internal use
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef DECAYTOPOLOGY2_H
#define DECAYTOPOLOGY2_H


#include <vector>

#include "particle.h"
#include "interactionVertex2.h"
#include "fsVertex.h"
#include "decayGraph.hpp"


namespace rpwa {


  typedef decayGraph<interactionVertex2, particle> decayTopologyGraphType;


  class decayTopology2 : public decayTopologyGraphType {
	
  public:
			
    decayTopology2();
    decayTopology2(const interactionVertexPtr&              productionVertex,
     		   const std::vector<interactionVertexPtr>& interactionVertices,
     		   const std::vector<particlePtr>&          fsParticles);
    decayTopology2(const decayTopology2&                    topo);
    decayTopology2(const decayTopologyGraphType&            graph);
    virtual ~decayTopology2();

    virtual decayTopology2& operator =(const decayTopology2&         topo);
    virtual decayTopology2& operator =(const decayTopologyGraphType& graph);
    virtual void clear();  ///< deletes all information
    
    virtual unsigned int nmbInteractionVertices() const { return _intVertices.size(); }  ///< returns number of interaction vertices
    virtual unsigned int nmbFsParticles        () const { return _fsParticles.size(); }  ///< returns number of final state particles

    virtual const std::vector<particlePtr>&          fsParticles        () const { return _fsParticles; }  ///< returns final state particles ordered depth-first
    virtual const std::vector<interactionVertexPtr>& interactionVertices() const { return _intVertices; }  ///< returns interaction vertices (excluding production vertex) ordered depth-first

    virtual const interactionVertexPtr& productionVertex  () const { return _prodVertex;     }  ///< returns production vertex
    virtual const interactionVertexPtr& xInteractionVertex() const { return _intVertices[0]; }  ///< returns X-decay vertex

    bool isProductionVertex (const interactionVertexPtr& vert) const { return (vert == _prodVertex); }
    bool isInteractionVertex(const interactionVertexPtr& vert) const;
    bool isFsVertex         (const interactionVertexPtr& vert) const;
    bool isFsParticle       (const particlePtr&          part) const;

    virtual bool checkTopology   () const;                  ///< returns whether decay has the correct topology
    virtual bool checkConsistency() const { return true; }  ///< checks conservation rules on all vertices

    virtual std::ostream& print(std::ostream& out) const;  ///< prints decay topology in human-readable form

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag

    
  //   //!!! bad quick hacks
  //   decayGraph subGraph(interactionVertex& startVertex);
  //   interactionVertex* vertex(particle& part);


  // protected:


  //   static decayGraph deepCopyGraph(const decayGraph& srcGraph,
  // 				    const bool        copyFsParticles = false);


  private:
    
    virtual decayTopology2& constructDecay(const interactionVertexPtr&              productionVertex,
					   const std::vector<interactionVertexPtr>& interactionVertices,
					   const std::vector<particlePtr>&          fsParticles);  ///< constructs the decay graph based on, production vertex, intermediate vertices, and final state particles


    interactionVertexPtr              _prodVertex;   ///< pointer to production vertex
    std::vector<interactionVertexPtr> _intVertices;  ///< array of interaction vertices excluding production vertex; ordered depth-first
    std::vector<particlePtr>          _fsParticles;  ///< array of final state particles; ordered depth-first
    
    static bool _debug;  ///< if set to true, debug messages are printed
    
  };
  

  inline
  std::ostream&
  operator <<(std::ostream&        out,
  	      const decayTopology2& topo)
  {
    return topo.print(out);
  }


}  // namespace rpwa


#endif  // DECAYTOPOLOGY2_H
