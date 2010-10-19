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
	BOOST_INSTALL_PROPERTY(graph,  bundle);
	enum vertex_VPtr_t  { vertex_VPtr  };
	BOOST_INSTALL_PROPERTY(vertex, VPtr);
	enum edge_PPtr_t    { edge_PPtr    };
	BOOST_INSTALL_PROPERTY(edge,   PPtr);
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
		                                   boost::GraphvizGraphProperty> graphProperties;
		typedef typename boost::property<boost::vertex_bundle_t,                        NBundleData,
		                                   boost::property<boost::vertex_VPtr_t,        VPtr,
		                                     boost::property<boost::vertex_index_t,     std::size_t,
		                                       boost::property<boost::vertex_name_t,    std::string,
		                                         boost::property<boost::vertex_color_t, color_t,
		                                           boost::GraphvizVertexProperty> > > > > nodeProperties;
		typedef typename boost::property<boost::edge_bundle_t,                      EBundleData,
		                                   boost::property<boost::edge_PPtr_t,      PPtr,
		                                     boost::property<boost::edge_name_t,    std::string,
		                                       boost::property<boost::edge_color_t, color_t,
		                                         boost::GraphvizEdgeProperty> > > > edgeProperties;
		// graph definition
		typedef typename boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS,
		                                       nodeProperties, edgeProperties, graphProperties> graphType;
		typedef typename boost::subgraph<graphType> graph;
		typedef typename boost::graph_traits<graph> graphTraits;


	public:
			
		// graph properties
		typedef typename boost::graph_property<graph, boost::graph_graph_attribute_t >::type graphAttributeType;      ///< type of graphViz graph attributes
		typedef typename boost::graph_property<graph, boost::graph_vertex_attribute_t>::type graphNodeAttributeType;  ///< type of graphViz node attributes
		typedef typename boost::graph_property<graph, boost::graph_edge_attribute_t  >::type graphEdgeAttributeType;  ///< type of graphViz edge attributes
		typedef typename boost::graph_property<graph, boost::graph_name_t            >::type graphNameType;           ///< type of graph name
		// node and edge descriptor types
		typedef typename graphTraits::vertex_descriptor nodeDesc;  ///< node descriptor type
		typedef typename graphTraits::edge_descriptor   edgeDesc;  ///< edge descriptor type
		// node and edge iterators
		typedef typename graphTraits::vertex_iterator    nodeIterator;     ///< node iterator type
		typedef typename graphTraits::adjacency_iterator adjIterator;      ///< node iterator type for adjacent nodes
		typedef typename graphTraits::edge_iterator      edgeIterator;     ///< edge iterator type
		typedef typename graphTraits::in_edge_iterator   inEdgeIterator;   ///< edge iterator type for edges going into a node
		typedef typename graphTraits::out_edge_iterator  outEdgeIterator;  ///< edge iterator type for edges coming out of a node
		// node property types
		typedef typename boost::property_map<graph, boost::vertex_bundle_t   >::type       nodeDataMapType;         ///< type of map [node descriptor] -> [node bundled property           ]
		typedef typename boost::property_map<graph, boost::vertex_VPtr_t     >::type       nodeVertexMapType;       ///< type of map [node descriptor] -> [vertex pointer property         ]
		typedef typename boost::property_map<graph, boost::vertex_VPtr_t     >::const_type nodeVertexMapConstType;  ///< type of map [node descriptor] -> [vertex pointer property         ]
		typedef typename boost::property_map<graph, boost::vertex_index_t    >::type       nodeIndexMapType;        ///< type of map [node descriptor] -> [node index property             ]
		typedef typename boost::property_map<graph, boost::vertex_index_t    >::const_type nodeIndexMapConstType;   ///< type of map [node descriptor] -> [node index property             ]
		typedef typename boost::property_map<graph, boost::vertex_name_t     >::type       nodeNameMapType;         ///< type of map [node descriptor] -> [node name property              ]
		typedef typename boost::property_map<graph, boost::vertex_name_t     >::const_type nodeNameMapConstType;    ///< type of map [node descriptor] -> [node name property              ]
		typedef typename boost::property_map<graph, boost::vertex_color_t    >::type       nodeColorMapType;        ///< type of map [node descriptor] -> [node color property             ]
		typedef typename boost::property_map<graph, boost::vertex_attribute_t>::type       nodeAttributeMapType;    ///< type of map [node descriptor] -> [node graphVis attribute property]
		// node edge property types
		typedef typename boost::property_map<graph, boost::edge_bundle_t   >::type       edgeDataMapType;           ///< type of map [edge descriptor] -> [edge bundled property           ]
		typedef typename boost::property_map<graph, boost::edge_PPtr_t     >::type       edgeParticleMapType;       ///< type of map [edge descriptor] -> [particle pointer property       ]
		typedef typename boost::property_map<graph, boost::edge_PPtr_t     >::const_type edgeParticleMapConstType;  ///< type of map [edge descriptor] -> [particle pointer property       ]
		typedef typename boost::property_map<graph, boost::edge_name_t     >::type       edgeNameMapType;           ///< type of map [edge descriptor] -> [edge name property              ]
		typedef typename boost::property_map<graph, boost::edge_name_t     >::const_type edgeNameMapConstType;      ///< type of map [edge descriptor] -> [edge name property              ]
		typedef typename boost::property_map<graph, boost::edge_color_t    >::type       edgeColorMapType;          ///< type of map [edge descriptor] -> [edge color property             ]
		typedef typename boost::property_map<graph, boost::edge_attribute_t>::type       edgeAttributeMapType;      ///< type of map [edge descriptor] -> [edge graphViz attribute property]
		typedef typename boost::property_map<graph, boost::edge_index_t    >::type       edgeIndexMapType;          ///< type of map [edge descriptor] -> [edge index property             ]
		typedef typename boost::property_map<graph, boost::edge_index_t    >::const_type edgeIndexMapConstType;     ///< type of map [edge descriptor] -> [edge index property             ]

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
			nodeIterator iNd, iNdEnd;
			_vertexNodeMap.clear();
			for (boost::tie(iNd, iNdEnd) = nodes(); iNd != iNdEnd; ++iNd) {
				const VPtr& vert = vertex(*iNd);
				if (not vert) {
					printErr << "encountered null vertex pointer for node " << *iNd
					         << " in graph '" << name() << "'. aborting." << std::endl;
					throw;
				}
				_vertexNodeMap[vert] = *iNd;
			}
			_particleEdgeMap.clear();
			edgeIterator iEd, iEdEnd;
			for (boost::tie(iEd, iEdEnd) = edges(); iEd != iEdEnd; ++iEd) {
				const PPtr& part = particle(*iEd);
				if (not part) {
					printErr << "encountered null particle pointer for edge " << *iEd
					         << " in graph '" << name() << "'. aborting." << std::endl;
					throw;
				}
				_particleEdgeMap[part] = *iEd;
			}
		}


		virtual
		VPtr
		cloneNode(const nodeDesc& nd,
		          const bool      cloneInParticles  = false,
		          const bool      cloneOutParticles = false)
		{
			const VPtr vert = vertex(nd);  // this must not be a reference
			if (not vert) {
				printErr << "encountered null vertex pointer for node " << nd
				         << " in graph '" << name() << "'. aborting." << std::endl;
				throw;
			}
			// clone vertex
			VPtr newVert(vert->clone(cloneInParticles, cloneOutParticles));
			// store new vertex in node property
			vertex(nd) = newVert;
			// update vertex-node map
			_vertexNodeMap.erase(vert);
			_vertexNodeMap[newVert] = nd;
			return newVert;
		}

		virtual
		PPtr
		cloneEdge(const edgeDesc& ed)
		{
			const PPtr part = particle(ed);  // this must not be a reference
			if (not part) {
				printErr << "encountered null particle pointer for edge " << ed
				         << " while cloning graph '" << name() << "'. aborting." << std::endl;
				throw;
			}
			// clone particle
			PPtr newPart(part->clone());
			// store new particle in node property
			particle(ed) = newPart;
			// replace particle in array of outgoing particles in fromVertex
			const VPtr& fromVert = fromVertex(ed);
			for (unsigned int i = 0; i < fromVert->nmbOutParticles(); ++i)
				if (fromVert->outParticles()[i] == part)
					fromVert->outParticles()[i] = newPart;
			// replace particle in array of incoming particles in toVertex
			const VPtr& toVert = toVertex(ed);
			for (unsigned int i = 0; i < toVert->nmbInParticles(); ++i)
				if (toVert->inParticles()[i] == part)
					toVert->inParticles()[i] = newPart;
			// update particle-edge map
			_particleEdgeMap.erase(part);
			_particleEdgeMap[newPart] = ed;
			return newPart;
		}


	protected:

		virtual
		decayGraph*
		doClone(const bool cloneParticles,
		        const bool) const  ///< helper function to use covariant return types with smart pointers; needed for public clone()
		{
			if (_debug)
				printInfo << "cloning graph '" << name() << "' "
				          << ((cloneParticles   ) ? "in" : "ex") << "cluding particles" << std::endl;
			// copy graph data structure
			decayGraph*  graphClone = new decayGraph(*this);
			nodeIterator iNd, iNdEnd;
			for (boost::tie(iNd, iNdEnd) = graphClone->nodes(); iNd != iNdEnd; ++iNd)
				graphClone->cloneNode(*iNd);
			if (cloneParticles) {
				edgeIterator iEd, iEdEnd;
				for (boost::tie(iEd, iEdEnd) = graphClone->edges(); iEd != iEdEnd; ++iEd)
					graphClone->cloneEdge(*iEd);
			}
			graphClone->buildReverseMaps();
			return graphClone;
		}


	public:

		inline
		boost::shared_ptr<decayGraph>
		clone(const bool cloneParticles = false) const  ///< creates deep copy of decay graph; must not be virtual
		{
			return boost::shared_ptr<decayGraph>(doClone(cloneParticles));
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
		addVertex(const VPtr& newVert)  ///< adds node to graph, creates edges according to particles in vertex, and stores vertex and particle pointers
		{
			if (not newVert) {
				printErr << "null pointer to vertex. aborting." << std::endl;
				throw;
			}
			// create graph node for vertex and store pointer
			nodeDesc newNd = boost::add_vertex(_graph);
			vertex        (newNd)   = newVert;
			_vertexNodeMap[newVert] = newNd;
			if (_debug)
				printInfo << "adding " << *newVert << " to graph '" << name() << "'"
				          << " as node [" << newNd << "]" << std::endl;
			// create edges for particles coming into vertex
			for (unsigned int iInPart = 0; iInPart < newVert->nmbInParticles(); ++iInPart) {
				const PPtr& part = newVert->inParticles()[iInPart];
				if (not part) {
					printErr << "null pointer to incoming particle[" << iInPart << "] "
					         << "in vertex " << *newVert << ". aborting."<< std::endl;
					throw;
				}
				// match with outgoing particles of other vertices
				nodeIterator iNd, iNdEnd;
				for (boost::tie(iNd, iNdEnd) = nodes(); iNd != iNdEnd; ++iNd) {
					if (*iNd == newNd)
						continue;
					const VPtr& otherVert = vertex(*iNd);
					for (unsigned int iOutPart = 0; iOutPart < otherVert->nmbOutParticles(); ++iOutPart)
						if (part == otherVert->outParticles()[iOutPart]) {
							bool     inserted;
							edgeDesc ed;
							boost::tie(ed, inserted) = boost::add_edge(*iNd, newNd, _graph);
							if (inserted) {
								particle        (ed)   = part;
								_particleEdgeMap[part] = ed;
								if (_debug)
									printInfo << "added edge for particle '" << part->name() << "' "
									          << "from node [" << *iNd << "] to node ["
									          << newNd << "] of graph '" << name() << "'" << std::endl;
							} else {
								printErr << "could not add edge for particle " << *part << " "
								         << "to graph '" << name() << "'. aborting." << std::endl;
								throw;
							}
						}
				}
			}
			// create edges for particles coming out of vertex
			for (unsigned int iOutPart = 0; iOutPart < newVert->nmbOutParticles(); ++iOutPart) {
				const PPtr& part = newVert->outParticles()[iOutPart];
				if (not part) {
					printErr << "null pointer to outgoing particle[" << iOutPart << "] "
					         << "in vertex " << *newVert << ". aborting."<< std::endl;
					throw;
				}
				// match with incoming particles of other vertices
				nodeIterator iNd, iNdEnd;
				for (boost::tie(iNd, iNdEnd) = nodes(); iNd != iNdEnd; ++iNd) {
					if (*iNd == newNd)
						continue;
					const VPtr& otherVert = vertex(*iNd);
					for (unsigned int iInPart = 0; iInPart < otherVert->nmbInParticles(); ++iInPart)
						if (part == otherVert->inParticles()[iInPart]) {
							bool     inserted;
							edgeDesc ed;
							boost::tie(ed, inserted) = boost::add_edge(newNd, *iNd, _graph);
							if (inserted) {
								particle        (ed)   = part;
								_particleEdgeMap[part] = ed;
								if (_debug)
									printInfo << "added edge for particle '" << part->name() << "' "
									          << "from node [" << *iNd << "] to node ["
									          << newNd << "] of graph '" << name() << "'" << std::endl;
							} else {
								printErr << "could not add edge for particle " << *part << " "
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
		inline const VPtr& vertex     (const nodeDesc& nd) const { return boost::get(boost::vertex_VPtr, _graph)[nd]; }  ///< returns vertex associated to node
		inline const VPtr& operator [](const nodeDesc& nd) const { return vertex(nd);                                 }  ///< returns vertex associated to node
		inline const PPtr& particle   (const edgeDesc& ed) const { return boost::get(boost::edge_PPtr,   _graph)[ed]; }  ///< returns particle associated to edge
		inline const PPtr& operator [](const edgeDesc& ed) const { return particle(ed);                               }  ///< returns particle associated to edge

		inline VPtr& vertex     (const nodeDesc& nd) { return boost::get(boost::vertex_VPtr, _graph)[nd]; }  ///< returns vertex associated to node
		inline VPtr& operator [](const nodeDesc& nd) { return vertex(nd);                                 }  ///< returns vertex associated to node
		inline PPtr& particle   (const edgeDesc& ed) { return boost::get(boost::edge_PPtr,   _graph)[ed]; }  ///< returns particle associated to edge
		inline PPtr& operator [](const edgeDesc& ed) { return particle(ed);                               }  ///< returns particle associated to edge
    
		inline bool isNode(const VPtr& v) const { return (_vertexNodeMap.find  (v) != _vertexNodeMap.end());   }  ///< returns whether vertex is associated to a node in the graph
    
		inline bool isEdge(const PPtr& p) const { return (_particleEdgeMap.find(p) != _particleEdgeMap.end()); }  ///< returns whether particle is associated to an edge in the graph
    

		inline
		nodeDesc
		node(const VPtr& v) const  ///< gives node descriptor associated to vertex
		{
			if (not v) {
				printErr << "null pointer for vertex. aborting." << std::endl;
				throw;
			}
			vertexNodeMapConstIt entry = _vertexNodeMap.find(v);
			if (entry == _vertexNodeMap.end()) {
				printErr << "vertex " << *v << " is not a node in graph '" << name() << "'. "
				         << "aborting." << std::endl;
				throw;
			}
			return entry->second;
		}
		inline nodeDesc operator [](const VPtr& v) const { return node(v); }  ///< gives node descriptor associated to vertex
    
		inline
		edgeDesc
		edge(const PPtr& p) const  ///< gives edge descriptor associated to particle
		{
			if (not p) {
				printErr << "null pointer for particle. aborting." << std::endl;
				throw;
			}
			particleEdgeMapConstIt entry = _particleEdgeMap.find(p);
			if (entry == _particleEdgeMap.end()) {
				printErr << "particle " << *p << " is not an edge in graph '" << name() << "'. "
				         << "aborting." << std::endl;
				throw;
			}
			//!!! for some reason the edge descriptor stored in the map does not work correctly,
			//    however, the one returned by boost::edge works fine
			//    I'm at a loss, because the behaviour is really strange:
			//    e.g. mysteriously ed != ed1
			//    or particle(ed) = something does _not_ change the particle pointer, whereas
			//    particle(ed1) = something works as expected
			//    at least the workaround is not too ugly
			const edgeDesc ed  = entry->second;
			const edgeDesc ed1 = edge(fromNode(ed), toNode(ed));
			return ed1;
		}
		inline edgeDesc operator [](const PPtr& p) const { return edge(p); }  ///< gives node descriptor associated to vertex


		// edge nodes and vertices
		inline nodeDesc fromNode(const edgeDesc& ed) const { return boost::source(ed, _graph); }  ///< returns descriptor of node where the given edge is coming from
		inline nodeDesc fromNode(const PPtr&     p ) const { return fromNode(edge(p));         }  ///< returns descriptor of node where the given particle is coming from
		inline nodeDesc toNode  (const edgeDesc& ed) const { return boost::target(ed, _graph); }  ///< returns descriptor of node where the given edge is going to
		inline nodeDesc toNode  (const PPtr&     p ) const { return toNode(edge(p));           }  ///< returns descriptor of node where the given particle is going to

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
		inline
		bool
		particleExists(const VPtr& fromV,
		               const VPtr& toV) const  ///< returns whether there is a particle going from fromV to toV
		{
			const nodeDesc fromNd = node(fromV);
			const nodeDesc toNd   = node(toV);
			return (edgeExists(fromNd, toNd) and particle(edge(fromNd, toNd)));
		}

		inline
		bool
		edgeConnects(const edgeDesc& ed,
		             const nodeDesc& fromNd,
		             const nodeDesc& toNd) const  ///< returns wether given edge goes from fromNd to toNd
		{
			std::pair<edgeDesc, bool> edFound = boost::edge(fromNd, toNd, _graph);
			if ((edFound.second) and (edFound.first == ed))
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
			if (not edFound.second)
				printWarn << "there is no edge connecting nodes[" << fromNd << "] and "
				          << "[" << toNd << "]. returning default edge." << std::endl;
			return edFound.first;
		}

		inline
		const PPtr&
		particle(const VPtr& fromV,
		         const VPtr& toV) const  ///< returns particle that goes from fromV to toV
		{
			std::pair<edgeDesc, bool> edFound = boost::edge(node(fromV), node(toV), _graph);
			if (edFound.second)
				return particle(edFound.first);
			return 0;
		}


		// property accessors
		nodeVertexMapType   nodeVertexMap  () { return boost::get(boost::vertex_VPtr,  _graph); }  ///< returns map [node descriptor] -> [node vertex pointer property  ]
		nodeIndexMapType    nodeIndexMap   () { return boost::get(boost::vertex_index, _graph); }  ///< returns map [node descriptor] -> [node index property           ]
		nodeNameMapType     nodeNameMap    () { return boost::get(boost::vertex_name,  _graph); }  ///< returns map [node descriptor] -> [node name property            ]
		edgeParticleMapType edgeParticleMap() { return boost::get(boost::edge_PPtr,    _graph); }  ///< returns map [edge descriptor] -> [edge particle pointer property]
		edgeIndexMapType    edgeIndexMap   () { return boost::get(boost::edge_index,   _graph); }  ///< returns map [edge descriptor] -> [edge index property           ]
		edgeNameMapType     edgeNameMap    () { return boost::get(boost::edge_name,    _graph); }  ///< returns map [edge descriptor] -> [edge name property            ]

		nodeVertexMapConstType   nodeVertexMap  () const { return boost::get(boost::vertex_VPtr,  _graph); }  ///< returns map [node descriptor] -> [node vertex pointer property  ]
		nodeIndexMapConstType    nodeIndexMap   () const { return boost::get(boost::vertex_index, _graph); }  ///< returns map [node descriptor] -> [node index property           ]
		nodeNameMapConstType     nodeNameMap    () const { return boost::get(boost::vertex_name,  _graph); }  ///< returns map [node descriptor] -> [node name property            ]
		edgeParticleMapConstType edgeParticleMap() const { return boost::get(boost::edge_PPtr,    _graph); }  ///< returns map [edge descriptor] -> [edge particle pointer property]
		edgeIndexMapConstType    edgeIndexMap   () const { return boost::get(boost::edge_index,   _graph); }  ///< returns map [edge descriptor] -> [edge index property           ]
		edgeNameMapConstType     edgeNameMap    () const { return boost::get(boost::edge_name,    _graph); }  ///< returns map [edge descriptor] -> [edge name property            ]

		inline const GBundleData& data ()                   const { return boost::get_property(_graph, boost::graph_bundle); }  ///< returns bundled graph property structure
		inline const NBundleData& data (const nodeDesc& nd) const { return boost::get(boost::vertex_bundle,   _graph)[nd];   }  ///< returns bundled property structure for given node
		inline const NBundleData& data (const VPtr&     v ) const { return data(node(v));                                    }  ///< returns bundled property structure for given vertex
		inline const EBundleData& data (const edgeDesc& ed) const { return boost::get(boost::edge_bundle,     _graph)[ed];   }  ///< returns bundled property structure for given edge
		inline const EBundleData& data (const PPtr&     p ) const { return date(edge(p));                                    }  ///< returns bundled property structure for given particle
		inline std::size_t        index(const nodeDesc& nd) const { return boost::get(boost::vertex_index,    _graph)[nd];   }  ///< returns index property for given node
		inline std::size_t        index(const VPtr&     v ) const { return index(node(v));                                   }  ///< returns index property for given vertex
		inline std::size_t        index(const edgeDesc& ed) const { return boost::get(boost::edge_index,      _graph)[ed];   }  ///< returns index property for given edge
		inline std::size_t        index(const PPtr&     p ) const { return index(edge(p));                                   }  ///< returns index property for given particle
		inline const std::string& name ()                   const { return boost::get_property(_graph, boost::graph_name);   }  ///< returns graph name property
		inline const std::string& name (const nodeDesc& nd) const { return boost::get(boost::vertex_name,     _graph)[nd];   }  ///< returns name property for given node
		inline const std::string& name (const VPtr&     v ) const { return name(node(v));                                    }  ///< returns name property for given vertex
		inline const std::string& name (const edgeDesc& ed) const { return boost::get(boost::edge_name,       _graph)[ed];   }  ///< returns name property for given edge
		inline const std::string& name (const PPtr&     p ) const { return name(edge(p));                                    }  ///< returns name property for given particle
		inline const color_t&     color(const nodeDesc& nd) const { return boost::get(boost::vertex_color,    _graph)[nd];   }  ///< returns color property for given node
		inline const color_t&     color(const VPtr&     v ) const { return color(node(v));                                   }  ///< returns color property for given vertex
		inline const color_t&     color(const edgeDesc& ed) const { return boost::get(boost::edge_color,      _graph)[ed];   }  ///< returns color property for given edge
		inline const color_t&     color(const PPtr&     p ) const { return color(edge(p));                                   }  ///< returns color property for given particle

		inline const boost::GraphvizAttrList& graphAttribute    () const { return boost::get_property(_graph, boost::graph_graph_attribute);   }  //< returns graphViz attribute property for graph
		inline const boost::GraphvizAttrList& graphNodeAttribute() const { return boost::get_property(_graph, boost::graph_vertex_attribute);  }  //< returns graphViz node attribute property for graph
		inline const boost::GraphvizAttrList& graphEdgeAttribute() const { return boost::get_property(_graph, boost::graph_edge_attribute);    }  //< returns graphViz edge attribute property for graph

		inline const boost::GraphvizAttrList& nodeAttribute(const nodeDesc& nd) const { return boost::get(boost::vertex_attribute, _graph)[nd]; }  //< returns graphViz attribute property for given node
		inline const boost::GraphvizAttrList& nodeAttribute(const VPtr&     v ) const { return nodeAttribute(node(v));                          }  ///< returns graphViz attribute property for given vertex
		inline const boost::GraphvizAttrList& edgeAttribute(const edgeDesc& ed) const { return boost::get(boost::edge_attribute,  _graph)[ed];  }  ///< returns graphViz attribute property for given edge
		inline const boost::GraphvizAttrList& edgeAttribute(const PPtr&     p ) const { return edgeAttribute(edge(p));                          }  ///< returns graphViz attribute property for given particle
    
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

		inline boost::GraphvizAttrList& graphAttribute    () { return boost::get_property(_graph, boost::graph_graph_attribute);   }  //< returns graphViz attribute property for graph
		inline boost::GraphvizAttrList& graphNodeAttribute() { return boost::get_property(_graph, boost::graph_vertex_attribute);  }  //< returns graphViz node attribute property for graph
		inline boost::GraphvizAttrList& graphEdgeAttribute() { return boost::get_property(_graph, boost::graph_edge_attribute);    }  //< returns graphViz edge attribute property for graph

		inline boost::GraphvizAttrList& nodeAttribute(const nodeDesc& nd) { return boost::get(boost::vertex_attribute, _graph)[nd]; }  //< returns graphViz attribute property for given node
		inline boost::GraphvizAttrList& nodeAttribute(const VPtr&     v ) { return nodeAttribute(node(v));                          }  ///< returns graphViz attribute property for given vertex
		inline boost::GraphvizAttrList& edgeAttribute(const edgeDesc& ed) { return boost::get(boost::edge_attribute,  _graph)[ed];  }  ///< returns graphViz attribute property for given edge
		inline boost::GraphvizAttrList& edgeAttribute(const PPtr&     p ) { return edgeAttribute(edge(p));                          }  ///< returns graphViz attribute property for given particle
    

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
				printInfo << "depth-first node order of graph '" << name() 
				          << "' starting at node[" << startNd << "]: ";
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


		nodeDesc
		topNode() const  ///< returns node which has no incoming edge and from which all other nodes can be reached
		{
			nodeIterator iNd, iNdEnd;
			for (boost::tie(iNd, iNdEnd) = nodes(); iNd != iNdEnd; ++iNd)
				if (nmbInEdges(*iNd) == 0) {
					std::vector<nodeDesc> sortedNds = sortNodesDfs(*iNd);
					if (sortedNds.size() == nmbNodes())
						return *iNd;
				}
			// could not find such node
			printWarn << "could not find top node in graph '" << name() << "'. "
			          << "returning default node." << std::endl;
			nodeDesc defaultNd = *iNdEnd;
			return defaultNd;
		}


		// subgraph routines
		decayGraph
		dfsSubGraph(const nodeDesc& startNd,
		            const bool      linkToMotherGraph = false)  ///< constructs subgraph that contains all nodes that can be reached from given start node
		{
			// find all nodes connected to startNd via depth-first search
			std::vector<nodeDesc> subGraphNds = sortNodesDfs(startNd);
			if (_debug)
				printInfo << "creating subgraph of graph '" << name() << "' starting at node " << startNd;
			decayGraph subDecayGraph;
			if (linkToMotherGraph) {
				//!!! with the current design of the decayGraph class
				//!!! graph-subgraph relations might be lost, because
				//!!! create_subgraph() creates a new graph and appends it to
				//!!! the list of children of this::_graph, however, the graph
				//!!! is then copied to subDecayGraph::_graph which breaks the
				//!!! parent-child connection; any modification to
				//!!! subDecayGraph::_graph will _not_ be reflected in the
				//!!! corresponding child graph of this::_graph
				//
				//!!! the only way around this would be to let the decayGraph
				//!!! class inherit from subgraph<graph> and overload
				//!!! subgraph<graph>::create_subgraph(); however none of the
				//!!! subgraph<graph> member functions are virtual...
				//
				//!!! another solution would be to make _graph a pointer to
				//!!! the underlying graph, but this would raise ownership and
				//!!! lifetime guarantee problems
				//
				//!!! also mind that every call to this function leaves a
				//!!! subgraph copy in the list of children in this::_graph;
				//!!! if not taken care of properly this might lead to
				//!!! memory-leak-like effects
				std::cout << " using nodes: " << std::flush;
				graph& subGraph = _graph.create_subgraph();
				for (unsigned int i = 0; i < subGraphNds.size(); ++i) {
					// add node from mother graph to subgraph
					boost::add_vertex(subGraphNds[i], subGraph);
					if (_debug)
						std::cout << subGraphNds[i] << "  ";
				}
				if (_debug)
					std::cout << std::endl;
				subDecayGraph = subGraph;
				subDecayGraph.buildReverseMaps();
			} else {
				if (_debug)
					std::cout << std::endl;
				// copy of part of mother graph; only vertex and particle
				// pointers (and none of the other properties) are copied
				for (unsigned int i = 0; i < subGraphNds.size(); ++i)
					// add node from mother graph to subgraph
					subDecayGraph.addVertex(vertex(subGraphNds[i]));
			}
			subDecayGraph.name() = "subgraph";
			return subDecayGraph;
		}

		inline
		decayGraph
		dfsSubGraph(const VPtr& startV)  ///< constructs subgraph that contains all nodes that can be reached from given start vertex
		{
			return dfsSubGraph(node(startV));
		}

		inline
		unsigned int
		nmbChildGraphs() const
		{
			return _graph.num_children();
		}


		void
		addGraph(const decayGraph& graph)  ///< copies all nodes and edges (including properties) in graph into this graph and creates additional edges for the possible connections between graph and this _graph
		{
			if (_debug)
				printInfo << "adding graph '" << graph.name() << "' "
				          << "to graph '" << name() << "'" << std::endl;
			boost::copy_graph(graph._graph, _graph);
			buildReverseMaps();
			createEdgesFromParticles();
		}


		// output member functions
		template<typename NLabel>
		std::ostream&
		print(std::ostream& out,
		      NLabel        nodeLabelMap) const  ///< prints decay graph nodes in human-readable form; node labels are determined from nodeLabelMap
		{
			out << "graph '" << name() << "' has " << nmbNodes() << " node(s):" << std::endl;
			nodeIterator iNd, iNdEnd;
			for (boost::tie(iNd, iNdEnd) = nodes(); iNd != iNdEnd; ++iNd) {
				out << "    " << boost::get(nodeLabelMap, *iNd) << " --> ";
				outEdgeIterator iEd, iEdEnd;
				for (boost::tie(iEd, iEdEnd) = boost::out_edges(*iNd, _graph);
				     iEd != iEdEnd; ++iEd)
					out << boost::get(nodeLabelMap, target(*iEd, _graph)) << " ";
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
			out << "graph '" << name() << "' has " << nmbEdges() << " edge(s):" << std::endl;
			edgeIterator iEd, iEdEnd;
			for (boost::tie(iEd, iEdEnd) = edges(); iEd != iEdEnd; ++iEd)
				out << "    " << boost::get(edgeLabelMap, *iEd)
				    << " = [" << boost::get(nodeLabelMap, fromNode(*iEd))
				    << ", "   << boost::get(nodeLabelMap, toNode  (*iEd)) << "]" << std::endl;
			return out;
		}

		virtual
		inline
		std::ostream&
		print(std::ostream& out) const  ///< prints decay graph in human-readable form using node and edge indices
		{ return print(out, nodeIndexMap(), edgeIndexMap()); }
    
		inline
		std::ostream&
		printPointers(std::ostream& out) const  ///< prints particle and vertex pointers stored in graph
		{ return print(out, nodeVertexMap(), edgeParticleMap()); }


		virtual
		inline
		std::ostream&
		writeGraphViz(std::ostream& out)  ///< writes graph in GraphViz DOT format
		{
			boost::write_graphviz(out, _graph);
			return out;
		}

		template<typename NLabel, typename ELabel>
		std::ostream&
		writeGraphViz(std::ostream& out,
		              const NLabel& nodeLabelMap,
		              const ELabel& edgeLabelMap)  ///< writes graph in GraphViz DOT format with node/edge labels defined by node/edgeLabelMap
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
					if (not edgeExists(*iFromNode, *iToNode))
						// add edge, if outParticles in fromNode matches inParticles in toNode
						for (unsigned int iOut = 0; iOut < vertex(*iFromNode)->nmbOutParticles(); ++iOut)
							for (unsigned int iIn = 0; iIn < vertex(*iToNode)->nmbInParticles(); ++iIn) {
								const PPtr& outPart = vertex(*iFromNode)->outParticles()[iOut];
								const PPtr& inPart  = vertex(*iToNode  )->inParticles ()[iIn ];
								if ((outPart == inPart) and (not isEdge(outPart))) {
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
