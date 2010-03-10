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
//      basic test program for vertex and decay topology
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <boost/config.hpp>
#include <boost/utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/connected_components.hpp>


#include "TVector3.h"

#include "utilities.h"
#include "particleDataTable.h"
#include "particle.h"
#include "diffractiveDissVertex.h"
#include "isobarDecayVertex.h"
#include "decayTopology.h"


using namespace std;
using namespace boost;
using namespace rpwa;


struct cycleDetector : public dfs_visitor<> {

  cycleDetector(bool& hasCycle) 
    : _hasCycle(hasCycle)
  { }

  template <class edge, class graph> void back_edge(edge, graph&) { _hasCycle = true; }

protected:

  bool& _hasCycle;

};


int
main(int argc, char** argv)
{

  if (1) {
    // switch on debug output
    //particleProperties::setDebug(true);
    //particleDataTable::setDebug(true);
    particle::setDebug(true);
    interactionVertex::setDebug(true);
  }

  particleDataTable& pdt = particleDataTable::instance();
  pdt.readFile();

  // test construction of vertices
  if (0) {
    TVector3 mom;
    mom = TVector3(1, 2, 3);
    particle beam("pi", -1, 0, mom);
    particle X("X", -1);
    X.setName("X");
    printInfo << "created particles: " << endl
	      << beam << endl
	      << X    << endl;
    diffractiveDissVertex vert1(beam, X);
    printInfo << "created vertex: " << endl
	      << vert1;
    diffractiveDissVertex vert2 = vert1;
    printInfo << "copied vertex: " << endl
	      << vert2;

    mom = TVector3(3, 4, 5);
    particle daughter1("pi", -1, 0, mom);
    mom = TVector3(4, 5, 6);
    particle daughter2("pi0", 0, -1, mom);
    isobarDecayVertex vert3(X, daughter1, daughter2, 1, 2);
    printInfo << "created vertex: " << endl
	      << vert3;
    isobarDecayVertex vert4 = vert3;
    printInfo << "copied vertex: " << endl
	      << vert4;
  }

  // BOOST Graph prototype
  if (1) {
    // define final state particles
    particle pi0("pi", -1);
    particle pi1("pi", +1);
    particle pi2("pi", -1);
    particle pi3("pi", +1);
    particle pi4("pi", -1);
    // define isobars
    particle sigma("sigma",     0);
    particle a1   ("a1(1269)", +1);
    particle f1   ("f1(1285)",  0);
    // define X-system
    //              q   I   G  2J  P   C  2M
    particle X("X", -1, 1, -1, 4, +1, +1, 2);
    // define vertices
    isobarDecayVertex vert0(X,     pi4, f1,    1, 2);
    isobarDecayVertex vert1(f1,    pi2, a1,    1, 2);
    isobarDecayVertex vert2(a1,    pi3, sigma, 1, 2);
    isobarDecayVertex vert3(sigma, pi0, pi1,   0, 0);

    // build graph
    vector<particle*> fsParticles;
    fsParticles.push_back(&pi0);
    fsParticles.push_back(&pi1);
    fsParticles.push_back(&pi2);
    fsParticles.push_back(&pi3);
    fsParticles.push_back(&pi4);
    const unsigned int nmbFsParticles = fsParticles.size();
    vector<interactionVertex*> decayVertices;
    decayVertices.push_back(&vert3);
    decayVertices.push_back(&vert1);
    decayVertices.push_back(&vert2);
    decayVertices.push_back(&vert0);
    const unsigned int nmbVertices = decayVertices.size();
    decayTopology topo(fsParticles, decayVertices);
    cout << endl;
    printInfo << "decay toplogy:" << endl
	      << topo << endl;

    // loop over all pairs of interaction vertices
    typedef adjacency_list<vecS, vecS, directedS,
                           property<vertex_name_t, interactionVertex*>,
                           property<edge_name_t,   particle*> > decayGraph;
    decayGraph g(nmbVertices + nmbFsParticles);
    typedef property_map<decayGraph, vertex_name_t>::type vertexPropType;
    typedef property_map<decayGraph, edge_name_t>::type   edgePropType;
    vertexPropType vertexProp = get(vertex_name, g);
    edgePropType   edgeProp   = get(edge_name,   g);
    typedef graph_traits<decayGraph>::vertex_descriptor vertexDesc;
    typedef graph_traits<decayGraph>::edge_descriptor   edgeDesc;
    for (unsigned int iFromVert = 0; iFromVert < nmbVertices; ++iFromVert)
      for (unsigned int iToVert = 0; iToVert < nmbVertices; ++iToVert) {
	interactionVertex* fromVertex = decayVertices[iFromVert];
	for (unsigned int iOutPart = 0; iOutPart < fromVertex->nmbOutParticles(); ++iOutPart) {
	  interactionVertex* toVertex = decayVertices[iToVert];
	  for (unsigned int iInPart = 0; iInPart < toVertex->nmbInParticles(); ++iInPart)
	    if (fromVertex->outParticles()[iOutPart] == toVertex->inParticles()[iInPart]) {
	      bool     inserted;
	      edgeDesc edge;
	      tie(edge, inserted) = add_edge(iFromVert, iToVert, g);
	      if (inserted) {
		edgeProp[edge]     = fromVertex->outParticles()[iOutPart];
		vertexDesc vertex  = source(edge, g);
		vertexProp[vertex] = fromVertex;
	      }
	    }
	}
      }
    // loop over all pairs of interaction vertices and final state particles
    for (unsigned int iVert = 0; iVert < nmbVertices; ++iVert) {
      interactionVertex* fromVertex = decayVertices[iVert];
      for (unsigned int iOutPart = 0; iOutPart < fromVertex->nmbOutParticles(); ++iOutPart)
	for (unsigned int iFsPart = 0; iFsPart < nmbFsParticles; ++iFsPart)
	  if (fromVertex->outParticles()[iOutPart] == fsParticles[iFsPart]) {
	    bool     inserted;
	    edgeDesc edge;
	    tie(edge, inserted) = add_edge(iVert, nmbVertices + iFsPart, g);
	    if (inserted) {
	      edgeProp[edge]     = fromVertex->outParticles()[iOutPart];
	      vertexDesc vertex  = source(edge, g);
	      vertexProp[vertex] = fromVertex;
	      vertex             = target(edge, g);
	      vertexProp[vertex] = 0;  // final state vertex
	    }
	  }
    }

    typedef graph_traits<decayGraph>::vertex_iterator    vertexIterator;
    typedef graph_traits<decayGraph>::adjacency_iterator adjIterator;
    typedef graph_traits<decayGraph>::edge_iterator      edgeIterator;
    {
      vertexIterator iVert, vertEnd;
      adjIterator    iAdj,  adjEnd;
      property_map<decayGraph, vertex_index_t>::type index_map = get(vertex_index, g);
      for (tie(iVert, vertEnd) = vertices(g); iVert != vertEnd; ++iVert) {
	const unsigned int i = get(index_map, *iVert);
	if (i < nmbVertices)
	  cout << "vertex[" << i << "] ";
	else
	  cout << "FS particle [" << i - nmbVertices << "] ";
	tie(iAdj, adjEnd) = adjacent_vertices(*iVert, g);
	if (iAdj == adjEnd)
	  cout << "has no adjacent vertex" << endl;
	else
	  cout << "---> ";
	for (; iAdj != adjEnd; ++iAdj) {
	  const unsigned int i = get(index_map, *iAdj);
	  if (i < nmbVertices)
	    cout << "vertex[" << i << "]";
	  else
	    cout << "FS particle [" << i - nmbVertices << "]";
	  if (next(iAdj) != adjEnd)
	    cout << " + ";
	  else
	    cout << endl;
	}
      }
    }

    {
      vertexIterator iVert, vertEnd;
      for (tie(iVert, vertEnd) = vertices(g); iVert != vertEnd; ++iVert) {
	interactionVertex* v = vertexProp[*iVert];
	cout << "vertex[" << *iVert << "] ";
	if (v) {
	  cout << "decay ";
	  for (unsigned int i = 0; i < v->nmbInParticles(); ++i) {
	    cout << v->inParticles()[i]->name();
	    if (i < v->nmbInParticles() - 1)
	      cout << " + ";
	  }
	  cout << " -> ";
	  for (unsigned int i = 0; i < v->nmbOutParticles(); ++i) {
	    cout << v->outParticles()[i]->name();
	    if (i < v->nmbOutParticles() - 1)
	      cout << " + ";
	  }
	  cout << endl;
	} else
	  cout << "final state" << endl;
      }
      edgeIterator iEdge, edgeEnd;
      for (tie(iEdge, edgeEnd) = edges(g); iEdge != edgeEnd; ++iEdge) {
	particle* p = edgeProp[*iEdge];
	cout << "particle[" << *iEdge << "] ";
	if (p) {
	  cout << p->name() << sign(p->charge()) << " ";
	  if (!vertexProp[target(*iEdge, g)])
	    cout << "final state particle";
	  cout << endl;
	} else
	  printErr << "zero pointer to particle" << endl;
      }
    }

    {
      typedef list<vertexDesc> vertexList;
      vertexList orderedVertices;
      topological_sort(g, front_inserter(orderedVertices));
      cout << "vertex order after topological sort: ";
      for (vertexList::iterator i = orderedVertices.begin(); i != orderedVertices.end(); ++i) {
	const unsigned int vIndex = *i;
	if (vIndex < nmbVertices)
	  cout << "vertex[" << vIndex << "]  ";
	else
	  cout << "FS particle [" << vIndex - nmbVertices << "]  ";
      }
      cout << endl;
    }

    {
      vector<int> components(num_vertices(g));
      const unsigned int nmbComponents = connected_components(g, &(components[0]));
      cout << "total number of components = " << nmbComponents << endl;
      for (unsigned int i = 0; i != components.size(); ++i) {
	if (i < nmbVertices)
	  cout << "vertex[" << i << "] ";
	else
	  cout << "FS particle [" << i - nmbVertices << "] ";
	cout << "is in component " << components[i] << endl;
      }
      cout << endl;
    }
    
    { // test for cycles
      bool          hasCycle = false;
      cycleDetector cycleDet(hasCycle);
      depth_first_search(g, visitor(cycleDet));
      cout << "the graph has a cycle " << hasCycle << endl;
      
      // add edges that introduce cycle
      add_edge(0, 3, g);

      hasCycle = false;
      depth_first_search(g, visitor(cycleDet));
      cout << "the graph has a cycle " << hasCycle << endl;
    }

  }

}




// #include <boost/config.hpp> // put this first to suppress some VC++ warnings

// #include <iostream>
// #include <iterator>
// #include <algorithm>
// #include <time.h>

// #include <boost/utility.hpp>
// #include <boost/graph/adjacency_list.hpp>
// #include <boost/graph/topological_sort.hpp>
// #include <boost/graph/depth_first_search.hpp>
// #include <boost/graph/dijkstra_shortest_paths.hpp>
// #include <boost/graph/visitors.hpp>

// using namespace std;
// using namespace boost;

// enum files_e { dax_h, yow_h, boz_h, zow_h, foo_cpp, 
//                foo_o, bar_cpp, bar_o, libfoobar_a,
//                zig_cpp, zig_o, zag_cpp, zag_o, 
//                  libzigzag_a, killerapp, N };
// const char* name[] = { "dax.h", "yow.h", "boz.h", "zow.h", "foo.cpp",
//                        "foo.o", "bar.cpp", "bar.o", "libfoobar.a",
//                        "zig.cpp", "zig.o", "zag.cpp", "zag.o",
//                        "libzigzag.a", "killerapp" };


// struct print_visitor : public bfs_visitor<> {
//   template <class Vertex, class Graph>
//   void discover_vertex(Vertex v, Graph&) {
//     cout << name[v] << " ";
//   }
// };




// int main(int,char*[])
// {

//   typedef pair<int,int> Edge;
//   Edge used_by[] = {
//     Edge(dax_h, foo_cpp), Edge(dax_h, bar_cpp), Edge(dax_h, yow_h),
//     Edge(yow_h, bar_cpp), Edge(yow_h, zag_cpp),
//     Edge(boz_h, bar_cpp), Edge(boz_h, zig_cpp), Edge(boz_h, zag_cpp),
//     Edge(zow_h, foo_cpp), 
//     Edge(foo_cpp, foo_o),
//     Edge(foo_o, libfoobar_a),
//     Edge(bar_cpp, bar_o),
//     Edge(bar_o, libfoobar_a),
//     Edge(libfoobar_a, libzigzag_a),
//     Edge(zig_cpp, zig_o),
//     Edge(zig_o, libzigzag_a),
//     Edge(zag_cpp, zag_o),
//     Edge(zag_o, libzigzag_a),
//     Edge(libzigzag_a, killerapp)
//   };
//   const std::size_t nedges = sizeof(used_by)/sizeof(Edge);

//   typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;
// #if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
//   // VC++ can't handle the iterator constructor
//   Graph g(N);
//   for (std::size_t j = 0; j < nedges; ++j) {
//     graph_traits<Graph>::edge_descriptor e; bool inserted;
//     tie(e, inserted) = add_edge(used_by[j].first, used_by[j].second, g);
//   }
// #else
//   Graph g(used_by, used_by + nedges, N);
// #endif
//   typedef graph_traits<Graph>::vertex_descriptor Vertex;

//   // Determine ordering for a full recompilation
//   // and the order with files that can be compiled in parallel
//   {
//     typedef list<Vertex> MakeOrder;
//     MakeOrder::iterator i;
//     MakeOrder make_order;


//     // Parallel compilation ordering
//     std::vector<int> time(N, 0);
//     for (i = make_order.begin(); i != make_order.end(); ++i) {    
//       // Walk through the in_edges an calculate the maximum time.
//       if (in_degree (*i, g) > 0) {
//         Graph::in_edge_iterator j, j_end;
//         int maxdist=0;
//         // Through the order from topological sort, we are sure that every 
//         // time we are using here is already initialized.
//         for (tie(j, j_end) = in_edges(*i, g); j != j_end; ++j)
//           maxdist=(std::max)(time[source(*j, g)], maxdist);
//         time[*i]=maxdist+1;
//       }
//     }

//     cout << "parallel make ordering, " << endl
//          << "vertices with same group number can be made in parallel" << endl;
//     {
//       graph_traits<Graph>::vertex_iterator i, iend;
//       for (tie(i,iend) = vertices(g); i != iend; ++i)
//         cout << "time_slot[" << name[*i] << "] = " << time[*i] << endl;
//     }

//   }
//   cout << endl;

//   // if I change yow.h what files need to be re-made?
//   {
//     cout << "A change to yow.h will cause what to be re-made?" << endl;
//     print_visitor vis;
//     breadth_first_search(g, vertex(yow_h, g), visitor(vis));
//     cout << endl;
//   }
//   cout << endl;

//   // are there any cycles in the graph?
//   {
//     bool has_cycle = false;
//     cycle_detector vis(has_cycle);
//     depth_first_search(g, visitor(vis));
//     cout << "The graph has a cycle? " << has_cycle << endl;
//   }
//   cout << endl;

//   // add a dependency going from bar.cpp to dax.h
//   {
//     cout << "adding edge bar_cpp -> dax_h" << endl;
//     add_edge(bar_cpp, dax_h, g);
//   }
//   cout << endl;

//   // are there any cycles in the graph?
//   {
//     bool has_cycle = false;
//     cycle_detector vis(has_cycle);
//     depth_first_search(g, visitor(vis));
//     cout << "The graph has a cycle now? " << has_cycle << endl;
//   }

//   return 0;
// }
