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
//      general class that represents decay process in form of a graph using
//      the Boost Graph Library
//      * each graph node is associated with a pointer to a vertex of type V
//      * each graph edge is associated with a pointer to a particle of type P
//      * optionally aribtrary structs of type GBundleData,
//        NBundleData, and EBundleData can be assigned to the graph,
//        each node, and each edge, respectively
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef DECAYGRAPH_HPP
#define DECAYGRAPH_HPP


#include <vector>
#include <map>

#include <boost/config.hpp>
#include <boost/utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/copy.hpp>
#include <boost/shared_ptr.hpp>

#include "utilities.h"


// add custom vertex and edge properties to BGL graph
namespace boost {
  enum graph_bundle_t { graph_bundle };
  BOOST_INSTALL_PROPERTY(graph, bundle);
  enum vertex_VPointer_t { vertex_VPointer };
  BOOST_INSTALL_PROPERTY(vertex, VPointer);
  enum edge_PPointer_t { edge_PPointer };
  BOOST_INSTALL_PROPERTY(edge, PPointer);
}


namespace rpwa {	


  struct emptyGraphBundleData {
  };  ///< empty default structure for graph, node, and edge bundled properties


  template<class V,
	   class P,
	   class GBundleData = emptyGraphBundleData,
	   class NBundleData = emptyGraphBundleData,
	   class EBundleData = emptyGraphBundleData>
  class decayGraph {

  public:

    // vertex and particle pointers
    typedef typename boost::shared_ptr<V> VPtr;
    typedef typename boost::shared_ptr<P> PPtr;


  private:
    
    // node and edge properties
    typedef typename boost::default_color_type color_t;
    typedef typename boost::property<boost::graph_bundle_t, GBundleData,
		       boost::property<boost::graph_name_t, std::string> > graphProperties;
    typedef typename boost::property<boost::vertex_bundle_t,        NBundleData,
		       boost::property<boost::vertex_VPointer_t,    VPtr,
		         boost::property<boost::vertex_index_t,     std::size_t,
		           boost::property<boost::vertex_name_t,    std::string,
			     boost::property<boost::vertex_color_t, color_t> > > > > nodeProperties;
    typedef typename boost::property<boost::edge_bundle_t,        EBundleData,
		       boost::property<boost::edge_PPointer_t,    PPtr,
		         boost::property<boost::edge_index_t,     std::size_t,
			   boost::property<boost::edge_name_t,    std::string,
			     boost::property<boost::edge_color_t, color_t> > > > > edgeProperties;
    // graph definition
    typedef typename boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS,
					   nodeProperties, edgeProperties, graphProperties> graphType;
    typedef typename boost::subgraph<graphType> graph;
    typedef typename boost::graph_traits<graph> graphTraits;


  public:
			
    // node and edge descriptor types
    typedef typename graphTraits::vertex_descriptor nodeDesc;  ///< node descriptor type
    typedef typename graphTraits::edge_descriptor   edgeDesc;  ///< edge descriptor type
    // node and edge iterators
    typedef typename graphTraits::vertex_iterator    nodeIterator;     ///< node iterator type
    typedef typename graphTraits::adjacency_iterator adjIterator;      ///< node iterator type for adjacent nodes
    typedef typename graphTraits::edge_iterator      edgeIterator;     ///< edge iterator type
    typedef typename graphTraits::in_edge_iterator   inEdgeIterator;   ///< edge iterator type for edges going into a node
    typedef typename graphTraits::out_edge_iterator  outEdgeIterator;  ///< edge iterator type for edges coming out of a node
    // node and edge property types
    typedef typename boost::property_map<graph, boost::vertex_bundle_t  >::type nodeDataMapType;      ///< type of map [node descriptor] -> [node bundled property    ]
    typedef typename boost::property_map<graph, boost::vertex_VPointer_t>::type nodeVertexMapType;    ///< type of map [node descriptor] -> [vertex pointer property  ]
    typedef typename boost::property_map<graph, boost::vertex_index_t   >::type nodeIndexMapType;     ///< type of map [node descriptor] -> [node index property      ]
    typedef typename boost::property_map<graph, boost::vertex_name_t    >::type nodeNameMapType;      ///< type of map [node descriptor] -> [node name property       ]
    typedef typename boost::property_map<graph, boost::vertex_color_t   >::type nodeColorMapType;     ///< type of map [node descriptor] -> [node color property      ]
    typedef typename boost::property_map<graph, boost::edge_bundle_t    >::type edgeDataMapType;      ///< type of map [edge descriptor] -> [edge bundled property    ]
    typedef typename boost::property_map<graph, boost::edge_PPointer_t  >::type edgeParticleMapType;  ///< type of map [edge descriptor] -> [particle pointer property]
    typedef typename boost::property_map<graph, boost::edge_index_t     >::type edgeIndexMapType;     ///< type of map [edge descriptor] -> [edge index property      ]
    typedef typename boost::property_map<graph, boost::edge_name_t      >::type edgeNameMapType;      ///< type of map [edge descriptor] -> [edge name property       ]
    typedef typename boost::property_map<graph, boost::edge_color_t     >::type edgeColorMapType;     ///< type of map [edge descriptor] -> [edge color property      ]

    typedef typename boost::property_map<graph, boost::vertex_index_t>::const_type nodeIndexMapConstType;  ///< type of map [node descriptor] -> [node index property]
    typedef typename boost::property_map<graph, boost::vertex_name_t >::const_type nodeNameMapConstType;   ///< type of map [node descriptor] -> [node name property ]
    typedef typename boost::property_map<graph, boost::edge_index_t  >::const_type edgeIndexMapConstType;  ///< type of map [edge descriptor] -> [edge index property]
    typedef typename boost::property_map<graph, boost::edge_name_t   >::const_type edgeNameMapConstType;   ///< type of map [edge descriptor] -> [edge name property ]

    
  private:

    // reverse map types
    typedef typename std::map<VPtr, nodeDesc>            vertexNodeMapType;
    typedef typename vertexNodeMapType::iterator         vertexNodeMapIt;
    typedef typename vertexNodeMapType::const_iterator   vertexNodeMapConstIt;
    typedef typename std::map<PPtr, edgeDesc>            particleEdgeMapType;
    typedef typename particleEdgeMapType::iterator       particleEdgeMapIt;
    typedef typename particleEdgeMapType::const_iterator particleEdgeMapConstIt;


  public:

    decayGraph() { }
    decayGraph(const decayGraph& graph) { *this = graph; }
    virtual ~decayGraph() { }

    inline
    decayGraph&
    operator =(const decayGraph& graph)
    {
      if (this != &graph) {
	_graph           = graph._graph;
	_vertexNodeMap   = graph._vertexNodeMap;
	_particleEdgeMap = graph._particleEdgeMap;
      }
      return *this;
    }


  protected:

    void
    buildReverseMaps()  ///< (re)builds vertex-node and particle-edge maps
    {
      nodeIterator iNode, iNodeEnd;
      _vertexNodeMap.clear();
      for (boost::tie(iNode, iNodeEnd) = nodes(); iNode != iNodeEnd; ++iNode) {
	const VPtr& v = vertex(*iNode);
	if (!v) {
	  printErr << "encountered null vertex pointer for node " << *iNode
		   << " in graph '" << name() << "'. aborting." << std::endl;
	  throw;
	}
	_vertexNodeMap[v] = *iNode;
      }
      _particleEdgeMap.clear();
      edgeIterator iEdge, iEdgeEnd;
      for (boost::tie(iEdge, iEdgeEnd) = edges(); iEdge != iEdgeEnd; ++iEdge) {
	const PPtr& p = particle(*iEdge);
	if (!p) {
	  printErr << "encountered null particle pointer for edge " << *iEdge
		   << " in graph '" << name() << "'. aborting." << std::endl;
	  throw;
	}
	_particleEdgeMap[p] = *iEdge;
      }
    }

    virtual
    VPtr
    cloneNode(const nodeDesc& nd)
    {
      const VPtr v = vertex(nd);
      if (!v) {
	printErr << "encountered null vertex pointer for node " << nd
		 << " in graph '" << name() << "'. aborting." << std::endl;
	throw;
      }
      // clone vertex
      if (_debug)
	printInfo << "cloning " << *v << std::endl;
      VPtr newV(v->clone(false, false));
      // store new vertex in node property
      vertex(nd) = newV;
      return newV;
    }

    virtual
    PPtr
    cloneEdge(const edgeDesc& ed)
    {
      const PPtr p = particle(ed);
      if (!p) {
	printErr << "encountered null particle pointer for edge " << ed
		 << " while coling graph '" << name() << "'. aborting." << std::endl;
	throw;
      }
      // clone particle
      if (_debug)
	printInfo << "cloning " << *p << std::endl;
      PPtr newP(p->clone());
      // store new particle in node property
      particle(ed) = newP;
      // store new particle in array of outgoing particles in fromVertex
      const VPtr& fromV = fromVertex(ed);
      for (unsigned int i = 0; i < fromV->nmbOutParticles(); ++i)
	if (fromV->outParticles()[i] == p)
	  fromV->outParticles()[i] = newP;
      // store new particle in array of incoming particles in toVertex
      const VPtr& toV = toVertex(ed);
      for (unsigned int i = 0; i < toV->nmbInParticles(); ++i)
	if (toV->inParticles()[i] == p)
	  toV->inParticles()[i] = newP;
      return newP;
    }


  public:

    virtual
    decayGraph*
    clone(const bool cloneParticles = true) const  ///< creates deep copy of graph
    {
      // copy graph data structure
      decayGraph* graphClone = new decayGraph(*this);
      nodeIterator iNode, iNodeEnd;
      for (boost::tie(iNode, iNodeEnd) = graphClone->nodes(); iNode != iNodeEnd; ++iNode)
	graphClone->cloneNode(*iNode);
      if (cloneParticles) {
	edgeIterator iEdge, iEdgeEnd;
	for (boost::tie(iEdge, iEdgeEnd) = graphClone->edges(); iEdge != iEdgeEnd; ++iEdge)
	  graphClone->cloneEdge(*iEdge);
      }
      graphClone->buildReverseMaps();
      return graphClone;
    }
    
    virtual
    void
    clear()  ///< deletes all information
    {
      //_graph.clear();  // not yet implemented; argh!
      // workaround
      _graph = graph();
      _vertexNodeMap.clear();
      _particleEdgeMap.clear();
    }


    nodeDesc
    addVertex(const VPtr& newV)  ///< adds node to graph, creates edges according to particles in vertex, and stores vertex and particle pointers
    {
      if (!newV) {
	printErr << "null pointer to vertex. aborting." << std::endl;
	throw;
      }
      // create graph node for vertex and store pointer
      nodeDesc newNd = boost::add_vertex(_graph);
      vertex(newNd) = newV;
      _vertexNodeMap[newV] = newNd;
      if (_debug)
	printInfo << "adding " << *newV << " to graph '" << name() << "'"
		  << " as node [" << newNd << "]" << std::endl;
      // create edges for particles coming into vertex
      for (unsigned int iInPart = 0; iInPart < newV->nmbInParticles(); ++iInPart) {
	const PPtr& p = newV->inParticles()[iInPart];
	if (!p) {
	  printErr << "null pointer to incoming particle[" << iInPart << "] "
		   << "in vertex " << *newV << ". aborting."<< std::endl;
	  throw;
	}
	// match with outgoing particles of other vertices
	nodeIterator iNode, iNodeEnd;
	for (boost::tie(iNode, iNodeEnd) = nodes(); iNode != iNodeEnd; ++iNode) {
	  if (*iNode == newNd)
	    continue;
	  const VPtr& otherV = vertex(*iNode);
	  for (unsigned int iOutPart = 0; iOutPart < otherV->nmbOutParticles(); ++iOutPart)
	    if (p == otherV->outParticles()[iOutPart]) {
	      bool     inserted;
	      edgeDesc ed;
	      boost::tie(ed, inserted) = boost::add_edge(*iNode, newNd, _graph);
	      if (inserted) {
		particle(ed)        = p;
		_particleEdgeMap[p] = ed;
		if (_debug)
		  printInfo << "added edge for particle '" << p->name() << "' "
			    << "from node [" << *iNode << "] to node ["
			    << newNd << "] of graph '" << name() << "'" << std::endl;
	      } else {
		printErr << "could not add edge for particle " << *p << " "
			 << "to graph '" << name() << "'. aborting." << std::endl;
		throw;
	      }
	    }
	}
      }
      // create edges for particles coming out of vertex
      for (unsigned int iOutPart = 0; iOutPart < newV->nmbOutParticles(); ++iOutPart) {
	const PPtr& p = newV->outParticles()[iOutPart];
	if (!p) {
	  printErr << "null pointer to outgoing particle[" << iOutPart << "] "
		   << "in vertex " << *newV << ". aborting."<< std::endl;
	  throw;
	}
	// match with incoming particles of other vertices
	nodeIterator iNode, iNodeEnd;
	for (boost::tie(iNode, iNodeEnd) = nodes(); iNode != iNodeEnd; ++iNode) {
	  if (*iNode == newNd)
	    continue;
	  const VPtr& otherV = vertex(*iNode);
	  for (unsigned int iInPart = 0; iInPart < otherV->nmbInParticles(); ++iInPart)
	    if (p == otherV->inParticles()[iInPart]) {
	      bool     inserted;
	      edgeDesc ed;
	      boost::tie(ed, inserted) = boost::add_edge(newNd, *iNode, _graph);
	      if (inserted) {
		particle(ed)        = p;
		_particleEdgeMap[p] = ed;
		if (_debug)
		  printInfo << "added edge for particle '" << p->name() << "' "
			    << "from node [" << *iNode << "] to node ["
			    << newNd << "] of graph '" << name() << "'" << std::endl;
	      } else {
		printErr << "could not add edge for particle " << *p << " "
			 << "to graph '" << name() << "'. aborting." << std::endl;
		throw;
	      }
	    }
	}
      }
      return newNd;
    }


    // iterators for nodes and edges
    inline unsigned int nmbNodes() const { return boost::num_vertices(_graph); }  ///< returns number of nodes in graph
    inline unsigned int nmbEdges() const { return boost::num_edges   (_graph); }  ///< returns number of edges in graph
    
    inline std::pair<nodeIterator, nodeIterator> nodes() const { return boost::vertices(_graph); }  ///< returns begin and end iterators for list of nodes
    inline std::pair<edgeIterator, edgeIterator> edges() const { return boost::edges   (_graph); }  ///< returns begin and end iterators for list of edges

    inline std::pair<adjIterator, adjIterator> adjacentVertices(const nodeDesc& nd) const
    { return boost::adjacent_vertices(nd, _graph); }  ///< returns begin and end iterators for list of nodes adjacent to given node
    inline std::pair<adjIterator, adjIterator> adjacentVertices(const VPtr&     v ) const
    { return adjacentVertices(node(v));            }  ///< returns begin and end iterators for list of nodes adjacent to given vertex

    inline unsigned int nmbInEdges (const nodeDesc& nd) const { return boost::in_degree (nd, _graph); }  ///< returns number of edges going into given node
    inline unsigned int nmbInEdges (const VPtr&     v ) const { return nmbInEdges (node(v));          }  ///< returns number of edges going into given vertex
    inline unsigned int nmbOutEdges(const nodeDesc& nd) const { return boost::out_degree(nd, _graph); }  ///< returns number of edges coming out of given node
    inline unsigned int nmbOutEdges(const VPtr&     v ) const { return nmbOutEdges(node(v));          }  ///< returns number of edges coming out of given vertex

    inline std::pair<inEdgeIterator, inEdgeIterator>   incomingEdges(const nodeDesc& nd) const
    { return boost::in_edges (nd, _graph); }  ///< returns begin and end iterators for list of edges going into the given node
    inline std::pair<inEdgeIterator, inEdgeIterator>   incomingEdges(const VPtr&     v ) const
    { return incomingEdges(node(v));       }  ///< returns begin and end iterators for list of edges going into the given vertex
    inline std::pair<outEdgeIterator, outEdgeIterator> outgoingEdges(const nodeDesc& nd) const
    { return boost::out_edges(nd, _graph); }  ///< returns begin and end iterators for list of edges coming out of the given node
    inline std::pair<outEdgeIterator, outEdgeIterator> outgoingEdges(const VPtr&     v ) const
    { return outgoingEdges(node(v));       }  ///< returns begin and end iterators for list of edges coming out of the given vertex


    // node and edge mapping
    inline const VPtr& vertex     (const nodeDesc& nd) const { return boost::get(boost::vertex_VPointer, _graph)[nd]; }  ///< returns vertex associated to node
    inline const VPtr& operator [](const nodeDesc& nd) const { return vertex(nd);                                     }  ///< returns vertex associated to node
    inline const PPtr& particle   (const edgeDesc& ed) const { return boost::get(boost::edge_PPointer,  _graph)[ed];  }  ///< returns particle associated to edge
    inline const PPtr& operator [](const edgeDesc& ed) const { return node(ed);                                       }  ///< returns particle associated to edge

    inline VPtr& vertex     (const nodeDesc& nd) { return boost::get(boost::vertex_VPointer, _graph)[nd]; }  ///< returns vertex associated to node
    inline VPtr& operator [](const nodeDesc& nd) { return vertex(nd);                                     }  ///< returns vertex associated to node
    inline PPtr& particle   (const edgeDesc& ed) { return boost::get(boost::edge_PPointer,  _graph)[ed];  }  ///< returns particle associated to edge
    inline PPtr& operator [](const edgeDesc& ed) { return particle(ed);                                   }  ///< returns particle associated to edge
    
    inline bool isNode(const VPtr& v) const  ///< returns whether vertex is associated to a node in the graph
    { return (_vertexNodeMap.find  (v) != _vertexNodeMap.end());   }
    inline bool isEdge(const PPtr& p) const  ///< returns whether particle is associated to an edge in the graph
    { return (_particleEdgeMap.find(p) != _particleEdgeMap.end()); }

    inline
    nodeDesc
    node(const VPtr& v) const  ///< gives node descriptor associated to vertex
    {
      if (!v) {
	printErr << "null pointer for vertex. aborting." << std::endl;
	throw;
      }
      vertexNodeMapConstIt pos = _vertexNodeMap.find(v);
      if (pos == _vertexNodeMap.end()) {
      	printErr << "vertex " << *v << " is not a node in graph '" << name() << "'. "
      		 << "aborting." << std::endl;
      	throw;
      }
      return pos->second;
    }
    inline nodeDesc operator [](const VPtr& v) const { return node(v); }  ///< gives node descriptor associated to vertex
    
    inline
    edgeDesc
    edge(const PPtr& p) const  ///< gives edge descriptor associated to particle
    {
      if (!p) {
	printErr << "null pointer for particle. aborting." << std::endl;
	throw;
      }
      particleEdgeMapConstIt pos = _particleEdgeMap.find(p);
      if (pos == _particleEdgeMap.end()) {
      	printErr << "particle " << *p << " is not an edge in graph '" << name() << "'. "
      		 << "aborting." << std::endl;
      	throw;
      }
      return pos->second;
    }
    inline edgeDesc operator [](const PPtr& p) const { return edge(p); }  ///< gives node descriptor associated to vertex


    // edge nodes and vertices
    inline nodeDesc fromNode(const edgeDesc& ed) const { return boost::source(ed, _graph); }  ///< returns descriptor of node where the given edge is coming from
    inline nodeDesc fromNode(const PPtr&     p ) const { return fromNode(edge(p));         }  ///< returns descriptor of node where the given particle is coming from
    inline nodeDesc toNode  (const edgeDesc& ed) const { return boost::target(ed, _graph); }  ///< returns descriptor of node where the given edge is going to
    inline nodeDesc toNode  (const PPtr&     p ) const { return toNode  (edge(p));         }  ///< returns descriptor of node where the given particle is going to

    inline const VPtr& fromVertex(const edgeDesc& ed) const { return vertex(fromNode(ed)); }  ///< returns vertex where the given edge is coming from
    inline const VPtr& fromVertex(const PPtr&     p ) const { return vertex(fromNode(p));  }  ///< returns vertex where the given particle is coming from
    inline const VPtr& toVertex  (const edgeDesc& ed) const { return vertex(toNode(ed));   }  ///< returns vertex where the given edge is going to
    inline const VPtr& toVertex  (const PPtr&     p ) const { return vertex(toNode(p));    }  ///< returns vertex where the given particle is going to

    inline VPtr& fromVertex(const edgeDesc& ed) { return vertex(fromNode(ed)); }  ///< returns vertex where the given edge is coming from
    inline VPtr& fromVertex(const PPtr&     p ) { return vertex(fromNode(p));  }  ///< returns vertex where the given particle is coming from
    inline VPtr& toVertex  (const edgeDesc& ed) { return vertex(toNode(ed));   }  ///< returns vertex where the given edge is going to
    inline VPtr& toVertex  (const PPtr&     p ) { return vertex(toNode(p));    }  ///< returns vertex where the given particle is going to

    inline
    bool
    edgeExists(const nodeDesc& fromNd,
	       const nodeDesc& toNd) const  ///< returns whether there is an edge going from fromNd to toNd
    {
      std::pair<edgeDesc, bool> edFound = boost::edge(fromNd, toNd, _graph);
      return edFound.second;
    }
    inline bool particleExists(const VPtr& fromV,
     			       const VPtr& toV) const { return edgeExists(node(fromV), node(toV)); }  ///< returns whether there is a particle going from fromV to toV

    inline
    bool
    edgeConnects(const edgeDesc& ed,
		 const nodeDesc& fromNd,
		 const nodeDesc& toNd) const  ///< returns wether given edge goes from fromNd to toNd
    {
      std::pair<edgeDesc, bool> edFound = boost::edge(fromNd, toNd, _graph);
      //!!! for some reason the edge descriptor comparison returns
      //    false although descriptors are the same
      printWarn << "!!! DOES NOT WORK: ed = " << ed << ", fromNd = " << fromNd << ", "
		<< "toNd = " << toNd << "; " << edFound.first << ", " << edFound.second << "; "
		<< (edFound.first == ed) << std::endl;
      if ((edFound.second) && (edFound.first == ed))
	return true;
      return false;	
    }
    inline bool particleConnects(const PPtr& p,
    				 const VPtr& fromV,
    				 const VPtr& toV) const { return edgeConnects(edge(p), node(fromV), node(toV)); }  ///< returns wether given particle goes from fromV to toV

    inline
    edgeDesc
    edge(const nodeDesc& fromNd,
	 const nodeDesc& toNd) const  ///< returns descriptor of edge that goes from fromNd to toNd
    {
      std::pair<edgeDesc, bool> edFound = boost::edge(fromNd, toNd, _graph);
      if (edFound.second)
	return edFound.first;
      //!!! not sure whether this makes sense
      edgeDesc defaultEd;
      return defaultEd;
    }
    inline
    const PPtr
    particle(const VPtr& fromV,
    	     const VPtr& toV) const  ///< returns particle that goes from fromV to toV
    {
      std::pair<edgeDesc, bool> edFound = boost::edge(node(fromV), node(toV), _graph);
      if (edFound.second)
	return edge(edFound.first);
      return 0;
    }


    // property accessors
    nodeIndexMapType nodeIndexMap() { return boost::get(boost::vertex_index, _graph); }  ///< returns map [node descriptor] -> [node index property]
    nodeNameMapType  nodeNameMap () { return boost::get(boost::vertex_name,  _graph); }  ///< returns map [node descriptor] -> [node name property ]
    edgeIndexMapType edgeIndexMap() { return boost::get(boost::edge_index,   _graph); }  ///< returns map [edge descriptor] -> [edge index property]
    edgeNameMapType  edgeNameMap () { return boost::get(boost::edge_name,    _graph); }  ///< returns map [edge descriptor] -> [edge name property ]

    nodeIndexMapConstType nodeIndexMap() const { return boost::get(boost::vertex_index, _graph); }  ///< returns map [node descriptor] -> [node index property]
    nodeNameMapConstType  nodeNameMap () const { return boost::get(boost::vertex_name,  _graph); }  ///< returns map [node descriptor] -> [node name property ]
    edgeIndexMapConstType edgeIndexMap() const { return boost::get(boost::edge_index,   _graph); }  ///< returns map [edge descriptor] -> [edge index property]
    edgeNameMapConstType  edgeNameMap () const { return boost::get(boost::edge_name,    _graph); }  ///< returns map [edge descriptor] -> [edge name property ]

    inline const GBundleData& data ()                   const { return boost::get_property(_graph, boost::graph_bundle); }  ///< returns bundled graph property structure
    inline const NBundleData& data (const nodeDesc& nd) const { return boost::get(boost::vertex_bundle,   _graph)[nd];   }  ///< returns bundled property structure for given node
    inline const NBundleData& data (const VPtr&     v ) const { return data(node(v));                                    }  ///< returns bundled property structure for given vertex
    inline const EBundleData& data (const edgeDesc& ed) const { return boost::get(boost::edge_bundle,     _graph)[ed];   }  ///< returns bundled property structure for given edge
    inline const EBundleData& data (const PPtr&     p ) const { return date(edge(p));                                    }  ///< returns bundled property structure for given particle
    inline const std::size_t  index(const nodeDesc& nd) const { return boost::get(boost::vertex_index,    _graph)[nd];   }  ///< returns index property for given node
    inline const std::size_t  index(const VPtr&     v ) const { return index(node(v));                                   }  ///< returns index property for given vertex
    inline const std::size_t  index(const edgeDesc& ed) const { return boost::get(boost::edge_index,      _graph)[ed];   }  ///< returns index property for given edge
    inline const std::size_t  index(const PPtr&     p ) const { return index(edge(p));                                   }  ///< returns index property for given particle
    inline const std::string& name ()                   const { return boost::get_property(_graph, boost::graph_name);   }  ///< returns graph name property
    inline const std::string& name (const nodeDesc& nd) const { return boost::get(boost::vertex_name,     _graph)[nd];   }  ///< returns name property for given node
    inline const std::string& name (const VPtr&     v ) const { return name(node(v));                                    }  ///< returns name property for given vertex
    inline const std::string& name (const edgeDesc& ed) const { return boost::get(boost::edge_name,       _graph)[ed];   }  ///< returns name property for given edge
    inline const std::string& name (const PPtr&     p ) const { return name(edge(p));                                    }  ///< returns name property for given particle
    inline const color_t&     color(const nodeDesc& nd) const { return boost::get(boost::vertex_color,    _graph)[nd];   }  ///< returns color property for given node
    inline const color_t&     color(const VPtr&     v ) const { return color(node(v));                                   }  ///< returns color property for given vertex
    inline const color_t&     color(const edgeDesc& ed) const { return boost::get(boost::edge_color,      _graph)[ed];   }  ///< returns color property for given edge
    inline const color_t&     color(const PPtr&     p ) const { return color(edge(p));                                   }  ///< returns color property for given particle               
    
    inline GBundleData& data ()                   { return boost::get_property(_graph, boost::graph_bundle); }  ///< returns bundled graph property structure
    inline NBundleData& data (const nodeDesc& nd) { return boost::get(boost::vertex_bundle,   _graph)[nd];   }  ///< returns bundled property structure for given node
    inline NBundleData& data (const VPtr&     v ) { return data(node(v));                                    }  ///< returns bundled property structure for given vertex
    inline EBundleData& data (const edgeDesc& ed) { return boost::get(boost::edge_bundle,     _graph)[ed];   }  ///< returns bundled property structure for given edge
    inline EBundleData& data (const PPtr&     p ) { return date(edge(p));                                    }  ///< returns bundled property structure for given particle
    inline std::string& name ()                   { return boost::get_property(_graph, boost::graph_name);   }  ///< returns graph name property
    inline std::string& name (const nodeDesc& nd) { return boost::get(boost::vertex_name,     _graph)[nd];   }  ///< returns name property for given node
    inline std::string& name (const VPtr&     v ) { return name(node(v));                                    }  ///< returns name property for given vertex
    inline std::string& name (const edgeDesc& ed) { return boost::get(boost::edge_name,       _graph)[ed];   }  ///< returns name property for given edge
    inline std::string& name (const PPtr&     p ) { return name(edge(p));                                    }  ///< returns name property for given particle
    inline color_t&     color(const nodeDesc& nd) { return boost::get(boost::vertex_color,    _graph)[nd];   }  ///< returns color property for given node
    inline color_t&     color(const VPtr&     v ) { return color(node(v));                                   }  ///< returns color property for given vertex
    inline color_t&     color(const edgeDesc& ed) { return boost::get(boost::edge_color,      _graph)[ed];   }  ///< returns color property for given edge
    inline color_t&     color(const PPtr&     p ) { return color(edge(p));                                   }  ///< returns color property for given particle


    // node and edge sorting routines
    std::vector<nodeDesc>
    sortNodesDfs(const nodeDesc& startNd) const  ///< returns nodes that can be reached from given start node sorted depth-first
    {
      std::vector<nodeDesc> sortedNds;
      std::vector<color_t>  colors(nmbNodes());
      boost::depth_first_visit(_graph, startNd,
			       makeDfsRecorder(std::back_inserter(sortedNds)),
			       boost::make_iterator_property_map(colors.begin(),
								 nodeIndexMap(), colors[0]));
      if (_debug) {
	printInfo << "depth-first node order starting at node[" << startNd << "]: ";
	for (unsigned int i = 0; i < sortedNds.size(); ++i)
	  std::cout << sortedNds[i] << ((i < sortedNds.size() - 1) ? ", " : "");
	std::cout << std::endl;
      }
      return sortedNds;
    }

    inline
    std::vector<nodeDesc>
    sortNodesDfs(const VPtr& startV) const  ///< returns nodes that can be reached from given start vertex sorted depth-first
    {
      return sortNodesDfs(node(startV));
    }


    // subgraph routines
    decayGraph
    dfsSubGraph(const nodeDesc& startNd)  ///< constructs subgraph that contains all nodes that can be reached from given start node
    {
      // find all nodes connected to startNd via depth-first search
      std::vector<nodeDesc> subGraphNds = sortNodesDfs(startNd);
      if (_debug)
    	printInfo << "creating subgraph of graph '" << name() << "' starting at node " << startNd
    		  << " using nodes: " << std::flush;
      decayGraph subGraph = _graph.create_subgraph();
      for (unsigned int i = 0; i < subGraphNds.size(); ++i) {
    	subGraph.addMotherGraphNode(subGraphNds[i]);
    	if (_debug)
    	  std::cout << subGraphNds[i] << "    ";
      }
      if (_debug)
    	std::cout << std::endl;
      subGraph.buildReverseMaps();
      return subGraph;
    }

    inline
    decayGraph
    dfsSubGraph(const VPtr& startV) const  ///< constructs subgraph that contains all nodes that can be reached from given start vertex
    {
      return dfsSubGraph(node(startV));
    }


    void
    addGraph(const decayGraph& graph)  ///< copies all nodes and edges (including properties) in graph into this graph and creates additional edges for the possible connections between graph and this _graph
    {
      boost::copy_graph(graph._graph, _graph);
      buildReverseMaps();
      createEdgesFromParticles();
    }


    decayGraph
    joinDaughterGraphs(const VPtr&                    motherVertex,
		       const std::vector<decayGraph>& daughterGraphs)  ///< joins
    {
      decayGraph newGraph;
      newGraph.addVertex(motherVertex);
      for (unsigned int i = 0; i < daughterGraphs.size(); ++i)
	newGraph.addGraph(daughterGraphs[i]);
      return newGraph;
    }

    decayGraph
    joinDaughterGraphs(const VPtr&       motherVertex,
		       const decayGraph& daughterGraph1,
		       const decayGraph& daughterGraph2)
    {
      std::vector<decayGraph> daughterGraphs(2);
      daughterGraphs[0] = daughterGraph1;
      daughterGraphs[1] = daughterGraph2;
      return joinDaughterGraphs(motherVertex, daughterGraphs);
    }


    // output member functions
    template<typename NLabel>
    std::ostream&
    print(std::ostream& out,
	  NLabel        nodeLabelMap) const  ///< prints decay graph nodes in human-readable form; node labels are determined from nodeLabelMap
    {
      out << "graph '" << name() << "' has " << nmbNodes() << " nodes:" << std::endl;
      nodeIterator iNode, iNodeEnd;
      for (boost::tie(iNode, iNodeEnd) = nodes(); iNode != iNodeEnd; ++iNode) {
      	out << "    " << boost::get(nodeLabelMap, *iNode) << " --> ";
      	outEdgeIterator iEdge, iEdgeEnd;
      	for (boost::tie(iEdge, iEdgeEnd) = boost::out_edges(*iNode, _graph);
	     iEdge != iEdgeEnd; ++iEdge)
      	  out << boost::get(nodeLabelMap, target(*iEdge, _graph)) << " ";
      	out << std::endl;
      }
      return out;
    }

    template<typename NLabel, typename ELabel>
    std::ostream&
    print(std::ostream& out,
	  NLabel        nodeLabelMap,
	  ELabel        edgeLabelMap) const  ///< prints decay graph nodes and edges in human-readable form; node/edge labels are determined from node/edgeLabelMap
    {
      print(out, nodeLabelMap);
      out << "graph '" << name() << "' has " << nmbEdges() << " edges:" << std::endl;
      edgeIterator iEdge, iEdgeEnd;
      for (boost::tie(iEdge, iEdgeEnd) = edges(); iEdge != iEdgeEnd; ++iEdge)
      	out << "    " << boost::get(edgeLabelMap, *iEdge)
      	    << " = [" << boost::get(nodeLabelMap, fromNode(*iEdge))
      	    << ", "   << boost::get(nodeLabelMap, toNode  (*iEdge)) << "]" << std::endl;
      return out;
    }

    inline
    std::ostream&
    print(std::ostream& out) const  ///< prints decay graph in human-readable form using node and edge indices
    { return print(out, nodeIndexMap(), edgeIndexMap()); }
    
    template<typename NLabel, typename ELabel>
    std::ostream&
    writeGraphViz(std::ostream& out,
		  const NLabel& nodeLabelMap,
		  const ELabel& edgeLabelMap) const  ///< writes graph in GraphViz DOT format with node/edge labels defined by node/edgeLabelMap
    {
      boost::write_graphviz(out, _graph,
      			    boost::make_label_writer(nodeLabelMap),
      			    boost::make_label_writer(edgeLabelMap));
      return out;
    }

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


  protected:

    void
    createEdgesFromParticles()  ///< (re)builds edges according to particles in vertices
    {
      nodeIterator iFromNode, iFromNodeEnd;
      for (boost::tie(iFromNode, iFromNodeEnd) = nodes(); iFromNode != iFromNodeEnd; ++iFromNode) {
	nodeIterator iToNode, iToNodeEnd;
	for (boost::tie(iToNode, iToNodeEnd) = nodes(); iToNode != iToNodeEnd; ++iToNode) {
	  if (!edgeExists(*iFromNode, *iToNode))
	    // add edge, if outParticles in fromNode matches inParticles in toNode
	    for (unsigned int iOut = 0; iOut < vertex(*iFromNode)->nmbOutParticles(); ++iOut)
	      for (unsigned int iIn = 0; iIn < vertex(*iToNode)->nmbInParticles(); ++iIn) {
		const PPtr& outPart = vertex(*iFromNode)->outParticles()[iOut];
		const PPtr& inPart  = vertex(*iToNode  )->inParticles ()[iIn];
		if ((outPart == inPart) && (!isEdge(outPart))) {
		  bool     inserted;
		  edgeDesc ed;
		  boost::tie(ed, inserted) = addEdge(*iFromNode, *iToNode);
		  if (inserted) {
		    particle(ed)              = outPart;
		    _particleEdgeMap[outPart] = ed;
		    if (_debug)
		      printInfo << "added edge for particle '" << outPart->name() << "' "
				<< "from node [" << *iFromNode << "] to node ["
				<< *iToNode << "] of graph '" << name() << "'" << std::endl;
		  } else {
		    printErr << "could not add edge for particle " << *outPart << " "
			     << "to graph '" << name() << "'. aborting." << std::endl;
		    throw;
		  }
		}
	      }
	}
      }
    }


  private:

    decayGraph(const graph& g) { *this = g; }

    inline
    decayGraph&
    operator =(const graph& g)
    {
      if (&_graph != &g) {
	_graph = g;
	buildReverseMaps();
      }
      return *this;
    }

    inline
    nodeDesc
    addMotherGraphNode(const nodeDesc& nd)  ///< adds a node from mother graph to graph
    {
      return boost::add_vertex(nd, _graph);
    }

    inline
    std::pair<edgeDesc, bool>
    addEdge(const nodeDesc& fromNd,
	    const nodeDesc& toNd)  ///< adds an edge to graph
    {
      return boost::add_edge(fromNd, toNd, _graph);
    }


    // recorder for depth-first visitor
    template<class nodeOutIterator> class dfsRecorder : public boost::default_dfs_visitor {
    public:
      dfsRecorder(const nodeOutIterator& nodeIt)
	: _nodeIt(nodeIt) { }
      template<typename node, typename graph> inline void discover_vertex(node n, const graph &)
      { *_nodeIt++ = n; }
    private:
      nodeOutIterator _nodeIt;
    };
    
    // object generator function
    template<class nodeOutIterator>
    inline
    dfsRecorder<nodeOutIterator>
    makeDfsRecorder(const nodeOutIterator& nodeIt) const
    { return dfsRecorder<nodeOutIterator>(nodeIt); }
    

    graph _graph;  ///< graph that represents particle decay

    vertexNodeMapType   _vertexNodeMap;    ///< maps vertex pointers to graph nodes
    particleEdgeMapType _particleEdgeMap;  ///< maps particle pointers to graph edges

    static bool _debug;  ///< if set to true, debug messages are printed
    
  };
  

  template<class V,
   	   class P,
   	   class GBundleData,
   	   class NBundleData,
   	   class EBundleData>
  bool decayGraph<V, P, GBundleData, NBundleData, EBundleData>::_debug = false;


  template<class V,
   	   class P,
   	   class GBundleData,
   	   class NBundleData,
   	   class EBundleData>
  inline
  std::ostream&
  operator <<(std::ostream& out,
	      const decayGraph<V, P, GBundleData, NBundleData, EBundleData>& graph)
  {
    return graph.print(out);
  }


}  // namespace rpwa


#endif  // DECAYGRAPH_HPP
