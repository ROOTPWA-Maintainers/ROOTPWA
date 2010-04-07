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


#include <fstream>
#include <sstream>

#include "TVector3.h"
#include "TSystem.h"

#include "svnVersion.h"
#include "utilities.h"
#include "particleDataTable.h"
#include "particle.h"
#include "diffractiveDissVertex.h"
#include "isobarDecayVertex.h"
#include "isobarDecayTopology.h"
#include "isobarDecayTopology2.h"
#include "decayGraph.hpp"
#include "diffractiveDissVertex2.h"
#include "isobarDecayVertex2.h"
#include "decayTopology2.h"


using namespace std;
using namespace rpwa;
using namespace boost;


struct graphBundleData {
  string foo;
};
struct nodeBundleData {
  double bar;
};
struct edgeBundleData {
  int blah;
};


typedef decayGraph<interactionVertex2, particle, graphBundleData,
		   nodeBundleData, edgeBundleData> graphType;


int
main(int argc, char** argv)
{
  printSvnVersion();

  // switch on debug output
  // particle::setDebug(true);
  // interactionVertex::setDebug(true);
  // diffractiveDissVertex::setDebug(true);
  // isobarDecayVertex::setDebug(true);
  // decayTopology::setDebug(true);
  // isobarDecayTopology::setDebug(true);
  graphType::setDebug(true);
  decayTopologyGraphType::setDebug(true);
  decayTopology2::setDebug(true);
  isobarDecayVertex2::setDebug(true);
  isobarDecayTopology2::setDebug(true);

  particleDataTable& pdt = particleDataTable::instance();
  pdt.readFile();


  if (1) {
    particlePtr pi0 = createParticle("pi-");
    particlePtr pi1 = createParticle("pi+");
    particlePtr pi2 = createParticle("pi-");
    particlePtr pi3 = createParticle("pi+");
    particlePtr pi4 = createParticle("pi-");
    // define isobars
    particlePtr sigma = createParticle("sigma");
    particlePtr a1    = createParticle("a1(1269)+");
    particlePtr f1    = createParticle("f1(1285)");
    // define X-system
    //                                   I   G  2J  P   C  2M
    particlePtr X = createParticle("X+", 2, +1, 4, +1, -1, 2);
    // define production vertex
    particlePtr beam = createParticle("pi-");
    diffractiveDissVertexPtr prodVert = createDiffractiveDissVertex(beam, X);
    // define vertices
    isobarDecayVertexPtr vert0 = createIsobarDecayVertex(X,     pi4, f1,    0, 3);
    isobarDecayVertexPtr vert1 = createIsobarDecayVertex(f1,    pi2, a1,    2, 2);
    isobarDecayVertexPtr vert2 = createIsobarDecayVertex(a1,    pi3, sigma, 2, 0);
    isobarDecayVertexPtr vert3 = createIsobarDecayVertex(sigma, pi0, pi1,   0, 0);

    // construct graph
    graphType g;
    g.name() = "test graph";
    g.addVertex(vert0);
    g.addVertex(vert1);
    g.addVertex(vert2);
    g.addVertex(vert3);
    cout << g;

    for (graphType::nodeIterator i = g.nodes().first; i != g.nodes().second; ++i) {
      const isobarDecayVertexPtr v = static_pointer_cast<isobarDecayVertex2>(g.vertex(*i));
      g.name(v) = v->mother()->name();
    }
    for (graphType::edgeIterator i = g.edges().first; i != g.edges().second; ++i) {
      const particlePtr& p = g.particle(*i);
      stringstream n;
      n << "edge [" << g.index(*i) << "] " << p->name();
      g.name(*i) = n.str();
    }
    g.print(cout, g.nodeNameMap(), g.edgeNameMap());
    {
      ofstream graphVizFile("decay.dot");
      g.writeGraphViz(graphVizFile, g.nodeNameMap(), g.edgeNameMap());
      gSystem->Exec("dot -Tps -o decay.ps decay.dot");
    }

    graphType g2  = g;
    g2.name()     = "graph copy";
    g2.data().foo = "foo";
    for (graphType::adjIterator i = g2.adjacentVertices(vert1).first;
    	 i != g2.adjacentVertices(vert1).second; ++i) {
      g2.data (*i).bar = *i + 0.5;
      g2.name (*i)    += " !bar!";
      g2.color(*i)     = white_color;
      cout << "vert1 adjacent vertex[" << *i << "]: " << *g2[*i] << endl;
    }
    cout << "nmbInParticles(vert1): " << g2.nmbInEdges(vert1) << endl;
    for (graphType::inEdgeIterator i = g2.incomingEdges(vert1).first;
    	 i != g2.incomingEdges(vert1).second; ++i) {
      g2.data(*i).blah = g2.index(*i) + 2;
      g2.name (*i)    += " !blah!";
      g2.color(*i)     = gray_color;
      cout << "vert1 in edge[" << *i << "]: " << *g2[*i] << endl;
    }
    cout << "nmbOutParticles(vert1): " << g2.nmbOutEdges(vert1) << endl;
    for (graphType::outEdgeIterator i = g2.outgoingEdges(vert1).first;
    	 i != g2.outgoingEdges(vert1).second; ++i) {
      g2.data(*i).blah = g2.index(*i) + 4;
      g2.name (*i)    += " !blah2!";
      g2.color(*i)     = black_color;
      cout << "vert1 out edge[" << *i << "]: " << *g2[*i] << endl;
    }
    isobarDecayVertexPtr dummyVert = createIsobarDecayVertex(X, pi4, f1, 0, 3);
    cout << "isNode(vert0) = "     << g2.isNode(vert0) << ", "
    	 << "isNode(dummyVert) = " << g2.isNode(dummyVert) << endl;
    particlePtr dummyPart = createParticle("X+", 1, +1, 4, +1, -1, 2);
    cout << "isEdge(f1) = "        << g2.isEdge(f1) << ", "
    	 << "isEdge(dummyPart) = " << g2.isEdge(dummyPart) << endl;
    cout << "fromVertex(f1): " << *g2.fromVertex(f1) << endl;
    cout << "toVertex(f1): "   << *g2.toVertex  (f1) << endl;
    cout << "particleExists(vert0, vert1) = " << g2.particleExists(vert0, vert1) << ", "
     	 << "particleExists(vert1, vert0) = " << g2.particleExists(vert1, vert0) << endl;

    //!!! does not work
    cout << "particleConnects(f1, vert0, vert1) = " << g2.particleConnects(f1, vert0, vert1) << endl;

    cout << "graph data = " << g2.data().foo << endl;
    for (graphType::nodeIterator i = g2.nodes().first; i != g2.nodes().second; ++i)
      cout << "node " << *i << ": data = " << g2.data(*i).bar << ", "
    	   << "name = '"  << g2.name(*i)  << "', "
    	   << "index = " << g2.index(*i) << ", "
    	   << "color = " << g2.color(*i) << ", "
    	   << endl;
    for (graphType::edgeIterator i = g2.edges().first; i != g2.edges().second; ++i)
      cout << "edge " << *i << ": data = " << g2.data(*i).blah << ", "
    	   << "name = '"  << g2.name(*i)  << "', "
    	   << "index = " << g2.index(*i) << ", "
    	   << "color = " << g2.color(*i) << ", "
    	   << endl;
    g2.print(cout, g2.nodeNameMap());
    g2.clear();
    printInfo << "after clear:" << endl << g2;

    decayTopologyGraphType g3;
    g3.name() = "test graph 3";
    g3.addVertex(vert0);
    g3.addVertex(vert1);
    g3.addVertex(vert2);
    g3.addVertex(vert3);
    g3.addVertex(prodVert);
    g3.addVertex(createFsVertex(pi0));
    g3.addVertex(createFsVertex(pi1));
    g3.addVertex(createFsVertex(pi2));
    g3.addVertex(createFsVertex(pi3));
    g3.addVertex(createFsVertex(pi4));
    cout << g3;
    decayTopology2 topo(g3);
    topo.name() = "topo graph copy";
    cout << topo;

    vector<isobarDecayVertexPtr> decayVertices;
    decayVertices.push_back(vert3);
    decayVertices.push_back(vert1);
    decayVertices.push_back(vert2);
    decayVertices.push_back(vert0);
    vector<particlePtr> fsParticles;
    fsParticles.push_back(pi0);
    fsParticles.push_back(pi1);
    fsParticles.push_back(pi2);
    fsParticles.push_back(pi3);
    fsParticles.push_back(pi4);
    isobarDecayTopology2 topo2(prodVert, decayVertices, fsParticles);
    topo2.name() = "topo";
    cout << topo2;
    const bool topologyOkay = topo2.checkTopology();
    cout << "topology okay = " << topologyOkay << endl;
    printInfo << endl << "performing consistency checks on isobar decay topology" << endl;
    const bool decayConsistent = topo2.checkConsistency();
    cout << "decay consistent = " << decayConsistent << endl;
    topo2.calcIsobarLzVec();

    {
      cout << endl << "testing cloning" << endl;
      isobarDecayTopology2 topo3 = *topo2.clone();
      isobarDecayVertexPtr v = topo3.isobarDecayVertices()[0];
      particlePtr          p = v->inParticles()[0];
      v->setL(2);
      v->setS(2);
      cout << *p << endl;
      p->setG(-1);
      p->setCharge(-1);
      p->setC(1);
      //topo3.checkConsistency();
      cout << topo3;
      decayTopologyGraphType::nodeIterator iNode, iNodeEnd;
      cout << "!!! topo2 vertex pointers: ";
      for (tie(iNode, iNodeEnd) = topo2.nodes(); iNode != iNodeEnd; ++iNode)
	cout << topo2.vertex(*iNode) << "    ";
      cout << endl;
      cout << "!!! topo3 vertex pointers: ";
      for (tie(iNode, iNodeEnd) = topo3.nodes(); iNode != iNodeEnd; ++iNode)
	cout << topo3.vertex(*iNode) << "    ";
      cout << endl;
      cout << "!!! topo3 interaction vertex pointers: ";
      for (unsigned int i = 0; i < topo3.nmbInteractionVertices(); ++i)
	cout << topo3.interactionVertices()[i] << "    ";
      cout << endl;
      cout << "!!! topo3 isobar decay vertex pointers: ";
      for (unsigned int i = 0; i < topo3.nmbInteractionVertices(); ++i)
	cout << topo3.isobarDecayVertices()[i] << "    ";
      cout << endl;
      decayTopologyGraphType::edgeIterator iEdge, iEdgeEnd;
      cout << "!!! topo2 particle pointers: ";
      for (tie(iEdge, iEdgeEnd) = topo2.edges(); iEdge != iEdgeEnd; ++iEdge)
	cout << topo2.particle(*iEdge) << "    ";
      cout << endl;
      cout << "!!! topo3 particle pointers: ";
      for (tie(iEdge, iEdgeEnd) = topo3.edges(); iEdge != iEdgeEnd; ++iEdge)
	cout << topo3.particle(*iEdge) << "    ";
      cout << endl;
      cout << "!!! topo3 FS particle pointers: ";
      for (unsigned int i = 0; i < topo3.nmbFsParticles(); ++i)
	cout << topo3.fsParticles()[i] << "    ";
      cout << endl;
    }

    {
      cout << endl << "testing subdecays and merge" << endl;
      isobarDecayTopology2 subDecay1 = topo2.subDecay(2);
      cout << "subdecay from node[2]: " << subDecay1;
      const isobarDecayTopology2& subDecay2 = topo2.subDecay(1);
      cout << "subdecay from node[1]: " << subDecay2;
      const isobarDecayTopology2& subDecay3 = topo2.subDecay(9);
      cout << "subdecay from node[9]: " << subDecay3;
      subDecay1.addDecay(subDecay3);
      cout << "subdecay from node[2] + subdecay from node[9]: " << subDecay1;
      isobarDecayTopology2 decay;
      decay.addVertex(vert0);
      subDecay1.addDecay(decay);
      cout << "subdecay from node[2] + subdecay from node[9] + vert0: " << subDecay1;
    }
    
    //topo2.possibleDecays();
  }

  // test construction of vertices
  if (0) {
    TVector3 mom;
    mom = TVector3(1, 2, 3);
    particle beam("pi-", mom);
    particle X("X-");
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
    particle daughter1("pi-", mom);
    mom = TVector3(4, 5, 6);
    particle daughter2("pi0", mom);
    isobarDecayVertex vert3(X, daughter1, daughter2, 1, 2);
    printInfo << "created vertex: " << endl
	      << vert3;
    isobarDecayVertex vert4 = vert3;
    printInfo << "copied vertex: " << endl
	      << vert4;
  }

  // test decay topology
  if (0) {
    // define final state particles
    particle pi0("pi-");
    particle pi1("pi+");
    particle pi2("pi-");
    particle pi3("pi+");
    particle pi4("pi-");
    // define isobars
    particle sigma("sigma");
    particle a1   ("a1(1269)+");
    particle f1   ("f1(1285)");
    // define X-system
    //               I   G  2J  P   C  2M
    particle X("X+", 1, +1, 4, +1, -1, 2);
    // define production vertex
    particle beam("pi-");
    diffractiveDissVertex prodVert(beam, X);
    // define vertices
    isobarDecayVertex vert0(X,     pi4, f1,    0, 3);
    isobarDecayVertex vert1(f1,    pi2, a1,    2, 2);
    isobarDecayVertex vert2(a1,    pi3, sigma, 2, 2);
    isobarDecayVertex vert3(sigma, pi0, pi1,   0, 0);

       
    // build graph
    vector<particle*> fsParticles;
    fsParticles.push_back(&pi0);
    fsParticles.push_back(&pi1);
    fsParticles.push_back(&pi2);
    fsParticles.push_back(&pi3);
    fsParticles.push_back(&pi4);
    vector<isobarDecayVertex*> decayVertices;
    decayVertices.push_back(&vert3);
    decayVertices.push_back(&vert1);
    decayVertices.push_back(&vert2);
    decayVertices.push_back(&vert0);
    isobarDecayTopology topo(fsParticles, prodVert, decayVertices);
    cout << endl;
    printInfo << "decay toplogy:" << endl
	      << topo << endl;
    
    const bool topologyOkay = topo.verifyTopology();
    cout << "topology okay = " << topologyOkay << endl << endl;

    printInfo << endl << "Performing consistency checks on Isobar-Vertices.." << endl;
    const bool consistency = topo.checkConsistency();
    cout << "consistency = " << consistency << endl << endl;


    for (unsigned int i = 0; i < topo.nmbFsParticles(); ++i)
      topo.fsParticles()[i]->setMomentum(TVector3(i + 1, 0, 0));
    printInfo << "updating Lorentz-vectors:" << endl
	      << topo.updateIsobarLzVec() << endl;
    
    ofstream graphVizFile("decay.dot");
    topo.writeGraphViz(graphVizFile);
    gSystem->Exec("dot -Tps -o decay.ps decay.dot");
  }
}
