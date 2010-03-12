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
//      internally the decay process is represented as a graph using
//      the Boost Graph Library
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef DECAYTOPOLOGY_H
#define DECAYTOPOLOGY_H


#include <vector>
#include <map>

#include <boost/config.hpp>
#include <boost/utility.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "particle.h"
#include "interactionVertex.h"


namespace rpwa {	

  class decayTopology {
	
  public:
			
    decayTopology();
    decayTopology(const decayTopology&                   topo);
    decayTopology(const std::vector<particle*>&          fsParticles,
		  interactionVertex&                     productionVertex,
		  const std::vector<interactionVertex*>& vertices);
    virtual ~decayTopology();

    virtual decayTopology& operator = (const decayTopology& topo);
    
    virtual decayTopology& constructDecay(const std::vector<particle*>&          fsParticles,
					  interactionVertex&                     productionVertex,
					  const std::vector<interactionVertex*>& vertices);  ///< constructs the decay graph based on final state particles and vertices

    virtual unsigned int nmbFsParticles() const { return _fsParticles.size(); }  ///< returne number of final state particles
    virtual unsigned int nmbVertices()    const { return _vertices.size();    }  ///< returns number of interaction vertices
    virtual std::vector<particle*>& fsParticles() { return _fsParticles; }  ///< returns final state particles
    virtual std::vector<interactionVertex*>& interactionVertices() { return _vertices; }  ///< returns all interaction vertices (excluding production vertex) ordered by depth first; first vertex is X-decay vertex
    virtual interactionVertex& productionVertex() { return *_productionVertex; }  ///< returns production vertex
    virtual interactionVertex& xInteractionVertex() { return *_vertices[0]; }  ///< returns X-decay vertex

    virtual bool dataAreValid()   const { return false; }  ///< indicates whether data are complete and valid
    virtual bool verifyTopology() const { return true;  }  ///< returns whether decay has the correct topology

    virtual std::ostream& print(std::ostream& out) const;  ///< prints decay topology in human-readable form
    virtual std::ostream& writeGraphViz(std::ostream& out) const;  ///< writes graph in GraphViz DOT format

    virtual void clear();  ///< deletes all information

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag
    

  protected:
    
    // graph definition
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS,
				  boost::property<boost::vertex_name_t, interactionVertex*,
				    boost::property<boost::vertex_color_t, boost::default_color_type> >,
				  boost::property<boost::edge_name_t, particle*> > decayGraph;
    typedef boost::graph_traits<decayGraph> graphTraits;
    // node and edge property types
    typedef boost::property_map<decayGraph, boost::vertex_name_t>::type  nodePropType;
    typedef boost::property_map<decayGraph, boost::vertex_index_t>::type nodeIndexType;
    typedef boost::property_map<decayGraph, boost::edge_name_t>::type    edgePropType;
    typedef boost::property_map<decayGraph, boost::edge_index_t>::type   edgeIndexType;
    // node and edge descriptor types
    typedef graphTraits::vertex_descriptor nodeDesc;
    typedef graphTraits::edge_descriptor   edgeDesc;
    // iterators
    typedef graphTraits::vertex_iterator    nodeIterator;
    typedef graphTraits::edge_iterator      edgeIterator;
    typedef graphTraits::adjacency_iterator adjIterator;
    typedef graphTraits::out_edge_iterator  outEdgeIterator;
    typedef graphTraits::in_edge_iterator   inEdgeIterator;
    
    decayGraph   _graph;     ///< graph that represents particle decay
    nodePropType _nodeProp;  ///< node properties; pointers to interaction vertices
    edgePropType _edgeProp;  ///< edge properties; pointer to particles
    
    std::map<interactionVertex*, nodeDesc> _vertexNodeMap;     ///< maps vertex pointers to graph nodes
    std::map<particle*,          edgeDesc> _particleEdgeMap;   ///< maps particle pointers to graph edges
    typedef std::map<interactionVertex*, nodeDesc>::iterator vertexNodeMapIt;
    typedef std::map<particle*,          edgeDesc>::iterator particleEdgeMapIt;
    

  private:
    
    interactionVertex*                     _productionVertex;  ///< pointer to production vertex
    std::vector<interactionVertex*>        _vertices;          ///< array of interaction vertices excluding production vertex; ordered depth-first
    std::vector<particle*>                 _fsParticles;       ///< number of final state particles
    
    static bool _debug;  ///< if set to true, debug messages are printed
    
    // recorder for depth first visitor
    template<class nodeOutIterator>
    class dfsRecorder;
    
    // object generator function
    template<class nodeOutIterator>
    dfsRecorder<nodeOutIterator> makeDfsRecorder(nodeOutIterator nodeIt)
    { return dfsRecorder<nodeOutIterator>(nodeIt); }

  };
  

  inline
  std::ostream&
  operator << (std::ostream&        out,
	       const decayTopology& topo) { return topo.print(out); }


} // namespace rpwa


#endif  // DECAYTOPOLOGY_H
