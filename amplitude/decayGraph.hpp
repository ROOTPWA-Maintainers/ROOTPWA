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
  };


  template<class V,
	   class P,
	   class GBundleData = emptyGraphBundleData,
	   class NBundleData = emptyGraphBundleData,
	   class EBundleData = emptyGraphBundleData>
  class decayGraph {
	
  private:

    // node and edge properties
    typedef typename boost::default_color_type color_t;
    typedef typename boost::property<boost::graph_bundle_t, GBundleData,
		       boost::property<boost::graph_name_t, std::string> > graphProperties;
    typedef typename boost::property<boost::vertex_bundle_t,        NBundleData,
		       boost::property<boost::vertex_VPointer_t,    V*,
		         boost::property<boost::vertex_index_t,     std::size_t,
		           boost::property<boost::vertex_name_t,    std::string,
			     boost::property<boost::vertex_color_t, color_t> > > > > nodeProperties;
    typedef typename boost::property<boost::edge_bundle_t,        EBundleData,
		       boost::property<boost::edge_PPointer_t,    P*,
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
    typedef typename graphTraits::vertex_descriptor nodeDesc;
    typedef typename graphTraits::edge_descriptor   edgeDesc;
    // node and edge iterators
    typedef typename graphTraits::vertex_iterator    nodeIterator;
    typedef typename graphTraits::adjacency_iterator adjIterator;
    typedef typename graphTraits::edge_iterator      edgeIterator;
    typedef typename graphTraits::in_edge_iterator   inEdgeIterator;
    typedef typename graphTraits::out_edge_iterator  outEdgeIterator;
    // node and edge property types
    typedef typename boost::property_map<graph, boost::vertex_bundle_t  >::type nodeDataMapType;
    typedef typename boost::property_map<graph, boost::vertex_VPointer_t>::type nodeVertexMapType;
    typedef typename boost::property_map<graph, boost::vertex_index_t   >::type nodeIndexMapType;
    typedef typename boost::property_map<graph, boost::vertex_name_t    >::type nodeNameMapType;
    typedef typename boost::property_map<graph, boost::vertex_color_t   >::type nodeColorMapType;
    typedef typename boost::property_map<graph, boost::edge_bundle_t    >::type edgeDataMapType;
    typedef typename boost::property_map<graph, boost::edge_PPointer_t  >::type edgeParticleMapType;
    typedef typename boost::property_map<graph, boost::edge_index_t     >::type edgeIndexMapType;
    typedef typename boost::property_map<graph, boost::edge_name_t      >::type edgeNameMapType;
    typedef typename boost::property_map<graph, boost::edge_color_t     >::type edgeColorMapType;

    typedef typename boost::property_map<graph, boost::vertex_index_t>::const_type nodeIndexMapConstType;
    typedef typename boost::property_map<graph, boost::vertex_name_t >::const_type nodeNameMapConstType;
    typedef typename boost::property_map<graph, boost::edge_index_t  >::const_type edgeIndexMapConstType;
    typedef typename boost::property_map<graph, boost::edge_name_t   >::const_type edgeNameMapConstType;

    
  private:

    // reverse map types
    typedef typename std::map<V*, nodeDesc>              vertexNodeMapType;
    typedef typename vertexNodeMapType::iterator         vertexNodeMapIt;
    typedef typename vertexNodeMapType::const_iterator   vertexNodeMapConstIt;
    typedef typename std::map<P*, edgeDesc>              particleEdgeMapType;
    typedef typename particleEdgeMapType::iterator       particleEdgeMapIt;
    typedef typename particleEdgeMapType::const_iterator particleEdgeMapConstIt;


  public:

    decayGraph() { }
    decayGraph(const decayGraph& graph) { *this = graph; }

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

    inline
    decayGraph&
    clone() const
    {
      decayGraph clone = _graph;
      return clone;
    }

    virtual ~decayGraph() { }

    // construct graph by adding vertices
    inline
    nodeDesc
    addVertex(V& newV)
    {
      // create graph node for vertex and store pointer
      nodeDesc newNd = boost::add_vertex(_graph);
      vertexPointer(newNd)  = &newV;
      _vertexNodeMap[&newV] = newNd;
      if (_debug)
	printInfo << "adding " << newV << " to graph '" << name() << "'"
		  << " as node [" << newNd << "]" << std::endl;
      // create edges for particles coming into vertex
      for (unsigned int iInPart = 0; iInPart < newV.nmbInParticles(); ++iInPart) {
	P* p = newV.inParticles()[iInPart];
	if (!p) {
	  printErr << "null pointer to incoming particle[" << iInPart << "] "
		   << "in vertex " << newV << ". aborting."<< std::endl;
	  throw;
	}
	// match with outgoing particles of other vertices
	nodeIterator iNode, iNodeEnd;
	for (boost::tie(iNode, iNodeEnd) = boost::vertices(_graph); iNode != iNodeEnd; ++iNode) {
	  if (*iNode == newNd)
	    continue;
	  const V& otherV = vertex(*iNode);
	  for (unsigned int iOutPart = 0; iOutPart < otherV.nmbOutParticles(); ++iOutPart)
	    if (p == otherV.outParticles()[iOutPart]) {
	      bool     inserted;
	      edgeDesc ed;
	      boost::tie(ed, inserted) = boost::add_edge(*iNode, newNd, _graph);
	      if (inserted) {
		particlePointer(ed) = p;
		_particleEdgeMap[p] = ed;
		if (_debug)
		  printInfo << "added edge from node [" << *iNode << "] to node ["
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
      for (unsigned int iOutPart = 0; iOutPart < newV.nmbOutParticles(); ++iOutPart) {
	P* p = newV.outParticles()[iOutPart];
	if (!p) {
	  printErr << "null pointer to outgoing particle[" << iOutPart << "] "
		   << "in vertex " << newV << ". aborting."<< std::endl;
	  throw;
	}
	// match with incoming particles of other vertices
	nodeIterator iNode, iNodeEnd;
	for (boost::tie(iNode, iNodeEnd) = boost::vertices(_graph); iNode != iNodeEnd; ++iNode) {
	  if (*iNode == newNd)
	    continue;
	  const V& otherV = vertex(*iNode);
	  for (unsigned int iInPart = 0; iInPart < otherV.nmbInParticles(); ++iInPart)
	    if (p == otherV.inParticles()[iInPart]) {
	      bool     inserted;
	      edgeDesc ed;
	      boost::tie(ed, inserted) = boost::add_edge(newNd, *iNode, _graph);
	      if (inserted) {
		particlePointer(ed) = p;
		_particleEdgeMap[p] = ed;
		if (_debug)
		  printInfo << "added edge from node [" << newNd << "] to node ["
			    << *iNode << "] of graph '" << name() << "'" << std::endl;
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

    inline unsigned int nmbNodes() const { return boost::num_vertices(_graph); }
    inline unsigned int nmbEdges() const { return boost::num_edges   (_graph); }
    
    // iterators for nodes and edges
    inline std::pair<edgeIterator, edgeIterator> edges() { return boost::edges   (_graph); }
    inline std::pair<nodeIterator, nodeIterator> nodes() { return boost::vertices(_graph); }

    inline std::pair<adjIterator, adjIterator> adjacentVertices(const nodeDesc& nd)
    { return boost::adjacent_vertices(nd, _graph); }
    inline std::pair<adjIterator, adjIterator> adjacentVertices(const V&        v )
    { return adjacentVertices(node(v));            }

    inline std::pair<inEdgeIterator, inEdgeIterator>   incomingEdges(const nodeDesc& nd)
    { return boost::in_edges (nd, _graph); }
    inline std::pair<inEdgeIterator, inEdgeIterator>   incomingEdges(const V&        v )
    { return incomingEdges(node(v));       }
    inline std::pair<outEdgeIterator, outEdgeIterator> outgoingEdges(const nodeDesc& nd)
    { return boost::out_edges(nd, _graph); }
    inline std::pair<outEdgeIterator, outEdgeIterator> outgoingEdges(const V&        v )
    { return outgoingEdges(node(v));       }

    // node and edge mapping
    inline const V& vertex     (const nodeDesc& nd) const { return *boost::get(boost::vertex_VPointer, _graph)[nd]; }
    inline const V& operator [](const nodeDesc& nd) const { return vertex(nd);                                      }
    inline const P& particle   (const edgeDesc& ed) const { return *boost::get(boost::edge_PPointer,   _graph)[ed]; }
    inline const P& operator [](const edgeDesc& ed) const { return node(ed);                                        }

    inline V& vertex     (const nodeDesc& nd) { return *boost::get(boost::vertex_VPointer, _graph)[nd]; }
    inline V& operator [](const nodeDesc& nd) { return vertex(nd);                                      }
    inline P& particle   (const edgeDesc& ed) { return *boost::get(boost::edge_PPointer,   _graph)[ed]; }
    inline P& operator [](const edgeDesc& ed) { return particle(ed);                                    }
    
    inline bool isNode(const V& v) const
    {
      //!!! the const cast is ugly, but I don't have a better solution
      //    at the moment; find() expects a non-const pointer as
      //    argument, but only to compare it to the pointers stored in
      //    the map
      return (_vertexNodeMap.find  (const_cast<V*>(&v)) != _vertexNodeMap.end());
    }
    inline bool isEdge(const P& p) const
    {
      //!!! the const cast is ugly, but I don't have a better solution
      //    at the moment; find() expects a non-const pointer as
      //    argument, but only to compare it to the pointers stored in
      //    the map
      return (_particleEdgeMap.find(const_cast<P*>(&p)) != _particleEdgeMap.end());
    }

    inline
    const nodeDesc
    node(const V& v) const
    {
      //!!! the const cast is ugly, but I don't have a better solution
      //    at the moment; find() expects a non-const pointer as
      //    argument, but only to compare it to the pointers stored in
      //    the map
      vertexNodeMapConstIt pos = _vertexNodeMap.find(const_cast<V*>(&v));
      if (pos == _vertexNodeMap.end()) {
	printErr << "vertex " << v << " is not a node in graph '" << name() << "'. "
		 << "aborting." << std::endl;
	throw;
      }
      return pos->second;
    }
    inline const nodeDesc operator [](const V& v) const { return node(v); }
    
    inline
    const edgeDesc
    edge(const P& p) const
    {
      //!!! the const cast is ugly, but I don't have a better solution
      //    at the moment; find() expects a non-const pointer as
      //    argument, but only to compare it to the pointers stored in
      //    the map
      particleEdgeMapConstIt pos = _particleEdgeMap.find(const_cast<P*>(&p));
      if (pos == _particleEdgeMap.end()) {
	printErr << "particle " << p << " is not an edge in graph '" << name() << "'. "
		 << "aborting." << std::endl;
	throw;
      }
      return pos->second;
    }
    inline const edgeDesc operator [](const P& p) const { return edge(p); }

    // edge nodes and vertices
    inline const nodeDesc fromNode(const edgeDesc& ed) const { return boost::source(ed, _graph); }
    inline const nodeDesc fromNode(const P&        p ) const { return fromNode(edge(p));         }
    inline const nodeDesc toNode  (const edgeDesc& ed) const { return boost::target(ed, _graph); }
    inline const nodeDesc toNode  (const P&        p ) const { return toNode  (edge(p));         }

    inline const V& fromVertex(const P&        p ) const { return vertex(fromNode(p));  }
    inline const V& fromVertex(const edgeDesc& ed) const { return vertex(fromNode(ed)); }
    inline const V& toVertex  (const P&        p ) const { return vertex(toNode(p));    }
    inline const V& toVertex  (const edgeDesc& ed) const { return vertex(toNode(ed));   }

    inline V& fromVertex(const P&        p ) { return vertex(fromNode(p));  }
    inline V& fromVertex(const edgeDesc& ed) { return vertex(fromNode(ed)); }
    inline V& toVertex  (const P&        p ) { return vertex(toNode(p));    }
    inline V& toVertex  (const edgeDesc& ed) { return vertex(toNode(ed));   }

    inline unsigned int nmbInEdges (const nodeDesc& nd) const { return boost::in_degree (nd, _graph); }
    inline unsigned int nmbOutEdges(const nodeDesc& nd) const { return boost::out_degree(nd, _graph); }

    inline unsigned int nmbInParticles (const V& v) const { return nmbInEdges (node(v)); }
    inline unsigned int nmbOutParticles(const V& v) const { return nmbOutEdges(node(v)); }

    inline
    bool
    edgeExists(const nodeDesc& fromNd,
	       const nodeDesc& toNd) const
    {
      if ((nmbOutEdges(fromNd) < 1) || (nmbInEdges(toNd) < 1))
	return false;
      adjIterator i, iEnd;
      for (boost::tie(i, iEnd) = boost::adjacent_vertices(fromNd, _graph); i != iEnd; ++i)
	if (*i == toNd)
	  return true;
      return false;
    }
    inline bool particleExists(const V& fromV,
			       const V& toV) const { return edgeExists(node(fromV), node(toV)); }

    inline
    bool
    edgeConnects(const edgeDesc& ed,
		 const nodeDesc& fromNd,
		 const nodeDesc& toNd) const
    {
      std::pair<edgeDesc, bool> edFound;
      edFound = boost::edge(fromNd, toNd, _graph);
      //!!! for some reason the edge descriptor comparison returns
      //    false although descriptors are the same
      printWarn << "!!! DOES NOT WORK: ed = " << ed << ", fromNd = " << fromNd << ", "
		<< "toNd = " << toNd << "; " << edFound.first << ", " << edFound.second << "; "
		<< (edFound.first == ed) << std::endl;
      if ((edFound.second) && (edFound.first == ed))
	return true;
      return false;	
    }
    inline bool particleConnects(const P& p,
				 const V& fromV,
				 const V& toV) const { return edgeConnects(edge(p), node(fromV), node(toV)); }

    inline
    const edgeDesc
    edge(const nodeDesc& fromNd,
	 const nodeDesc& toNd) const
    {
      outEdgeIterator i, iEnd;
      for (boost::tie(i, iEnd) = boost::out_edges(fromNd, _graph); i != iEnd; ++i)
	if (boost::target(*i, _graph) == toNd)
	  return *i;
      //!!! not sure whether this makes sense
      edgeDesc defaultEd;
      return defaultEd;
    }

    inline
    const P*
    particle(const V& fromV,
	     const V& toV) const
    {
      outEdgeIterator i, iEnd;
      nodeDesc        toNd = node(toV);
      for (boost::tie(i, iEnd) = boost::out_edges(node(fromV), _graph); i != iEnd; ++i)
	if (boost::target(*i, _graph) == toNd)
	  return &particle(*i);
      return 0;
    }

    // property accessors
    nodeIndexMapType nodeIndexMap() { return boost::get(boost::vertex_index, _graph); }  ///< returns map that maps graph node descriptors to node indices
    nodeNameMapType  nodeNameMap () { return boost::get(boost::vertex_name,  _graph); }  ///< returns map that maps graph node descriptors to node name
    edgeIndexMapType edgeIndexMap() { return boost::get(boost::edge_index,   _graph); }  ///< returns map that maps graph edge descriptors to edge indices
    edgeNameMapType  edgeNameMap () { return boost::get(boost::edge_name,    _graph); }  ///< returns map that maps graph edge descriptors to edge name

    nodeIndexMapConstType nodeIndexMap() const { return boost::get(boost::vertex_index, _graph); }  ///< returns map that maps graph node descriptors to node indices
    nodeNameMapConstType  nodeNameMap () const { return boost::get(boost::vertex_name,  _graph); }  ///< returns map that maps graph node descriptors to node name
    edgeIndexMapConstType edgeIndexMap() const { return boost::get(boost::edge_index,   _graph); }  ///< returns map that maps graph edge descriptors to edge indices
    edgeNameMapConstType  edgeNameMap () const { return boost::get(boost::edge_name,    _graph); }  ///< returns map that maps graph edge descriptors to edge name

    inline const GBundleData& data ()                   const { return boost::get_property(_graph, boost::graph_bundle); }
    inline const NBundleData& data (const nodeDesc& nd) const { return boost::get(boost::vertex_bundle,   _graph)[nd];   }
    inline const NBundleData& data (const V&        v ) const { return data(node(v));                                    }
    inline const EBundleData& data (const edgeDesc& ed) const { return boost::get(boost::edge_bundle,     _graph)[ed];   }
    inline const EBundleData& data (const P&        p ) const { return date(edge(p));                                    }
    inline const std::string& name ()                   const { return boost::get_property(_graph, boost::graph_name);   }
    inline const std::string& name (const nodeDesc& nd) const { return boost::get(boost::vertex_name,     _graph)[nd];   }
    inline const std::string& name (const V&        v ) const { return name(node(v));                                    }
    inline const std::string& name (const edgeDesc& ed) const { return boost::get(boost::edge_name,       _graph)[ed];   }
    inline const std::string& name (const P&        p ) const { return name(edge(p));                                    }
    inline const std::size_t  index(const nodeDesc& nd) const { return boost::get(boost::vertex_index,    _graph)[nd];   }
    inline const std::size_t  index(const V&        v ) const { return index(node(v));                                   }
    inline const std::size_t  index(const edgeDesc& ed) const { return boost::get(boost::edge_index,      _graph)[ed];   }
    inline const std::size_t  index(const P&        p ) const { return index(edge(p));                                   }
    inline const color_t&     color(const nodeDesc& nd) const { return boost::get(boost::vertex_color,    _graph)[nd];   }
    inline const color_t&     color(const V&        v ) const { return color(node(v));                                   }
    inline const color_t&     color(const edgeDesc& ed) const { return boost::get(boost::edge_color,      _graph)[ed];   }
    inline const color_t&     color(const P&        p ) const { return color(edge(p));                                   }
    
    inline GBundleData& data ()                   { return boost::get_property(_graph, boost::graph_bundle); }
    inline NBundleData& data (const nodeDesc& nd) { return boost::get(boost::vertex_bundle,   _graph)[nd];   }
    inline NBundleData& data (const V&        v ) { return data(node(v));                                    }
    inline EBundleData& data (const edgeDesc& ed) { return boost::get(boost::edge_bundle,     _graph)[ed];   }
    inline EBundleData& data (const P&        p ) { return date(edge(p));                                    }
    inline std::string& name ()                   { return boost::get_property(_graph, boost::graph_name);   }
    inline std::string& name (const nodeDesc& nd) { return boost::get(boost::vertex_name,     _graph)[nd];   }
    inline std::string& name (const V&        v ) { return name(node(v));                                    }
    inline std::string& name (const edgeDesc& ed) { return boost::get(boost::edge_name,       _graph)[ed];   }
    inline std::string& name (const P&        p ) { return name(edge(p));                                    }
    inline color_t&     color(const nodeDesc& nd) { return boost::get(boost::vertex_color,    _graph)[nd];   }
    inline color_t&     color(const V&        v ) { return color(node(v));                                   }
    inline color_t&     color(const edgeDesc& ed) { return boost::get(boost::edge_color,      _graph)[ed];   }
    inline color_t&     color(const P&        p ) { return color(edge(p));                                   }

    // subgraph routines
    decayGraph
    dfsSubGraph(const nodeDesc& startNd) const
    {
      // find all nodes connected to startNd via depth-first search
      std::vector<nodeDesc>                  subGraphNds;
      std::vector<boost::default_color_type> colors(nmbNodes());
      boost::depth_first_visit(_graph, startNd,
    			       makeDfsRecorder(std::back_inserter(subGraphNds)),
    			       boost::make_iterator_property_map(colors.begin(),
    								 boost::get(boost::vertex_index, _graph), colors[0]));
      if (_debug)
    	printInfo << "creating subgraph of graph '" << name() << "' starting at node " << startNd
    		  << " using nodes: " << std::flush;
      decayGraph subGraph = _graph.create_subgraph();
      for (unsigned int i = 0; i < subGraphNds.size(); ++i) {
    	boost::add_vertex(subGraphNds[i], subGraph);
    	if (_debug)
    	  std::cout << subGraphNds[i] << "    ";
      }
      if (_debug)
    	std::cout << std::endl;
      return subGraph;
    }
    

    // virtual void clear();  ///< deletes all information

    template<typename NLabel>
    std::ostream&
    print(std::ostream& out,
	  NLabel        nodeLabelMap) const  ///< prints decay graph nodes in human-readable form; node labels are determined from nodeLabelMap
    {
      out << "graph '" << name() << "' has " << nmbNodes() << " nodes:" << std::endl;
      nodeIterator iNode, iNodeEnd;
      for (boost::tie(iNode, iNodeEnd) = boost::vertices(_graph); iNode != iNodeEnd; ++iNode) {
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
      for (boost::tie(iEdge, iEdgeEnd) = boost::edges(_graph); iEdge != iEdgeEnd; ++iEdge)
      	out << "    " << boost::get(edgeLabelMap, *iEdge)
      	    << " = [" << boost::get(nodeLabelMap, fromNode(*iEdge))
      	    << ", " << boost::get(nodeLabelMap, toNode  (*iEdge)) << "]" << std::endl;
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


  private:

    decayGraph(const graph& g) { *this = g; }

    inline
    decayGraph&
    operator =(const graph& g)
    {
      if (_graph != &g) {
	_graph = g;
	// rebuild reverse maps
	nodeIterator iNode, iNodeEnd;
	for (boost::tie(iNode, iNodeEnd) = boost::vertices(_graph); iNode != iNodeEnd; ++iNode) {
	  V* v = vertexPointer(*iNode);
	  if (!v) {
	    printErr << "encountered null vertex pointer for node " << *iNode
		     << " when copying graph '" << g.name() << "'. aborting." << std::endl;
	    throw;
	  }
	  _vertexNodeMap[v] = *iNode;
	}
	edgeIterator iEdge, iEdgeEnd;
	for (boost::tie(iEdge, iEdgeEnd) = boost::edges(_graph); iEdge != iEdgeEnd; ++iEdge) {
	  P* p = particle(*iNode);
	  if (!p) {
	    printErr << "encountered null particle pointer for edge " << *iEdge
		     << " when copying graph '" << g.name() << "'. aborting." << std::endl;
	    throw;
	  }
	  _particleEdgeMap[particle(*iEdge)] = *iEdge;
	}
      }
      return *this;
    }

    inline const V*& vertexPointer  (const nodeDesc& nd) const { return boost::get(boost::vertex_VPointer, _graph)[nd]; }
    inline const P*& particlePointer(const edgeDesc& ed) const { return boost::get(boost::edge_PPointer,   _graph)[ed]; }

    inline V*& vertexPointer  (const nodeDesc& nd) { return boost::get(boost::vertex_VPointer, _graph)[nd]; }
    inline P*& particlePointer(const edgeDesc& ed) { return boost::get(boost::edge_PPointer,   _graph)[ed]; }

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
    makeDfsRecorder(const nodeOutIterator& nodeIt) { return dfsRecorder<nodeOutIterator>(nodeIt); }
    
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
  { return graph.print(out); }


} // namespace rpwa


#endif  // DECAYGRAPH_HPP
