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

#include <boost/config.hpp>
#include <boost/utility.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "particle.h"
#include "interactionVertex.h"


namespace rpwa {	

  class decayTopology {
	
  public:
			
    decayTopology();
    decayTopology(const decayTopology& topo);
    virtual ~decayTopology();

    decayTopology& operator = (const decayTopology& topo);
    
    decayTopology& constructDecay(const std::vector<particle*>&          fsParticles,
				  const std::vector<interactionVertex*>& vertices);  ///< constructs the dacy graph based on final state particles and vertices

    unsigned int nmbVertices()    { return _nmbVertices;    }  ///< returns number of interaction vertices
    unsigned int nmbFsParticles() { return _nmbFsParticles; }  ///< returne number of final state particles

    bool dataAreValid() const { return false; }  ///< indicates whether data are complete and valid

    std::ostream& print(std::ostream& out) const;  ///< prints decay topology in human-readable form

    void clear();  ///< deletes all information

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


  private:

    // graph definition
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
				  boost::property<boost::vertex_name_t, interactionVertex*>,
				  boost::property<boost::edge_name_t,   particle*> > decayGraph;
    // vertex and edge property types
    typedef boost::property_map<decayGraph, boost::vertex_name_t>::type vertexPropType;
    typedef boost::property_map<decayGraph, boost::edge_name_t>::type   edgePropType;
    // vertex and edge descriptor types
    typedef boost::graph_traits<decayGraph>::vertex_descriptor vertexDesc;
    typedef boost::graph_traits<decayGraph>::edge_descriptor   edgeDesc;
    // iterators
    typedef boost::graph_traits<decayGraph>::vertex_iterator    vertexIterator;
    typedef boost::graph_traits<decayGraph>::edge_iterator      edgeIterator;
    typedef boost::graph_traits<decayGraph>::adjacency_iterator adjIterator;
    
    decayGraph     _graph;           ///< graph that represents particle decay
    vertexPropType _vertexProp;      ///< vertex properties; pointers to interaction vertices
    edgePropType   _edgeProp;        ///< edge properties; pointer to particles
    unsigned int   _nmbVertices;     ///< number of interaction vertices
    unsigned int   _nmbFsParticles;  ///< number of final state particles

    static bool _debug;  ///< if set to true, debug messages are printed

  };


  inline
  std::ostream&
  operator << (std::ostream&        out,
	       const decayTopology& topo) { return topo.print(out); }


} // namespace rpwa


#endif  // DECAYTOPOLOGY_H
