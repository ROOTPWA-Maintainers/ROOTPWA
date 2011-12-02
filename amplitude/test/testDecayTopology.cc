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

#include "particleDataTable.h"
#include "../particle.h"
#include "decayGraph.hpp"
#include "decayTopology.h"
#include "diffractiveDissVertex.h"
#include "isobarDecayVertex.h"
#include "isobarDecayTopology.h"


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


typedef decayGraph<interactionVertex, particle, graphBundleData,
                   nodeBundleData, edgeBundleData> graphType;


int
main(int argc, char** argv)
{
	printSvnVersion();

	// switch on debug output
	particle::setDebug(true);
	graphType::setDebug(true);
	decayTopologyGraphType::setDebug(true);
	decayTopology::setDebug(true);
	interactionVertex::setDebug(true);
	diffractiveDissVertex::setDebug(true);
	isobarDecayVertex::setDebug(true);
	isobarDecayTopology::setDebug(true);

	particleDataTable& pdt = particleDataTable::instance();
	pdt.readFile();

	// test construction of vertices
	if (0) {
		TVector3 mom;
		mom = TVector3(1, 2, 3);
		particlePtr beam = createParticle("pi-");
		beam->setMomentum(mom);
		particlePtr target = createParticle("p+");
		particlePtr X      = createParticle("X-");
		printInfo << "created particles: " << endl
		          << *beam   << endl
		          << *target << endl
		          << *X      << endl;
		diffractiveDissVertexPtr vert1 = createDiffractiveDissVertex(beam, target, X);
		printInfo << "created vertex: " << endl
		          << *vert1 << endl;
		diffractiveDissVertexPtr vert2(vert1);
		printInfo << "copied vertex: " << endl
		          << *vert2 << endl;

		mom = TVector3(3, 4, 5);
		particlePtr daughter1 = createParticle("pi-");
		daughter1->setMomentum(mom);
		mom = TVector3(4, 5, 6);
		particlePtr daughter2 = createParticle("pi0");
		daughter2->setMomentum(mom);
		isobarDecayVertexPtr vert3 = createIsobarDecayVertex(X, daughter1, daughter2, 1, 2);
		printInfo << "created vertex: " << endl
		          << *vert3 << endl;
		isobarDecayVertexPtr vert4(vert3);
		printInfo << "copied vertex: " << endl
		          << *vert4 << endl;
	}

	if (1) {
		particlePtr pi0 = createParticle("pi-");
		particlePtr pi1 = createParticle("pi+");
		particlePtr pi2 = createParticle("pi-");
		particlePtr pi3 = createParticle("pi+");
		particlePtr pi4 = createParticle("pi-");
		// define isobars
		particlePtr sigma = createParticle("sigma");
		particlePtr a1    = createParticle("a1(1260)+");
		particlePtr f1    = createParticle("f1(1285)");
		// define X-system
		//                                   2I  G  2J  P   C  2M
		particlePtr X = createParticle("X+", 2, +1, 4, +1, -1, 2);
		// define production vertex
		particlePtr              beam     = createParticle("pi-");
		particlePtr              target   = createParticle("p+");
		diffractiveDissVertexPtr prodVert = createDiffractiveDissVertex(beam, target, X);
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
		printInfo << g;

		for (graphType::nodeIterator i = g.nodes().first; i != g.nodes().second; ++i) {
			const isobarDecayVertexPtr v = static_pointer_cast<isobarDecayVertex>(g.vertex(*i));
			g.name(v) = v->parent()->name();
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
			printInfo << "vert1 adjacent vertex[" << *i << "]: " << *g2[*i] << endl;
		}
		printInfo << "nmbInParticles(vert1): " << g2.nmbInEdges(vert1) << endl;
		for (graphType::inEdgeIterator i = g2.incomingEdges(vert1).first;
		     i != g2.incomingEdges(vert1).second; ++i) {
			g2.data(*i).blah = g2.index(*i) + 2;
			g2.name (*i)    += " !blah!";
			g2.color(*i)     = gray_color;
			printInfo << "vert1 in edge[" << *i << "]: " << *g2[*i] << endl;
		}
		printInfo << "nmbOutParticles(vert1): " << g2.nmbOutEdges(vert1) << endl;
		for (graphType::outEdgeIterator i = g2.outgoingEdges(vert1).first;
		     i != g2.outgoingEdges(vert1).second; ++i) {
			g2.data(*i).blah = g2.index(*i) + 4;
			g2.name (*i)    += " !blah2!";
			g2.color(*i)     = black_color;
			printInfo << "vert1 out edge[" << *i << "]: " << *g2[*i] << endl;
		}
		isobarDecayVertexPtr dummyVert = createIsobarDecayVertex(X, pi4, f1, 0, 3);
		printInfo << "isNode(vert0) = "     << g2.isNode(vert0) << ", "
		          << "isNode(dummyVert) = " << g2.isNode(dummyVert) << endl;
		particlePtr dummyPart = createParticle("X+", 1, +1, 4, +1, -1, 2);
		printInfo << "isEdge(f1) = "        << g2.isEdge(f1) << ", "
		          << "isEdge(dummyPart) = " << g2.isEdge(dummyPart) << endl;
		cout << "fromVertex(f1): " << *g2.fromVertex(f1) << endl;
		cout << "toVertex(f1): "   << *g2.toVertex  (f1) << endl;
		cout << "particleExists(vert0, vert1) = " << g2.particleExists(vert0, vert1) << ", "
		     << "particleExists(vert1, vert0) = " << g2.particleExists(vert1, vert0) << endl;
		cout << "particleConnects(f1, vert0, vert1) = " << g2.particleConnects(f1, vert0, vert1) << endl;
		cout << "particleConnects(a1, vert0, vert1) = " << g2.particleConnects(a1, vert0, vert1) << endl;
		printInfo << "graph data = " << g2.data().foo << endl;
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
		printInfo << g3;
		decayTopology topo(g3);
		topo.name() = "topo graph copy";
		printInfo << topo;

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
		isobarDecayTopology topo2(prodVert, decayVertices, fsParticles);
		topo2.name() = "topo";
		printInfo << topo2;
		const bool topologyOkay = topo2.checkTopology();
		printInfo << "topology okay = " << topologyOkay << endl;
		printInfo << "performing consistency checks on isobar decay topology" << endl;
		const bool decayConsistent = topo2.checkConsistency();
		printInfo << "decay consistent = " << decayConsistent << endl;
		topo2.calcIsobarLzVec();

		{
			cout << "-------------------------------------------------------------------------------"
			     << endl;
			printInfo << "testing cloning" << endl;
			isobarDecayTopologyPtr topo3 = topo2.clone();
			isobarDecayVertexPtr   v     = topo3->isobarDecayVertices()[0];
			particlePtr            p     = v->inParticles()[0];
			v->setL(2);
			v->setS(2);
			cout << *p << endl;
			p->setG(-1);
			p->setCharge(-1);
			p->setC(1);
			//topo3->checkTopology();
			//topo3->checkConsistency();
			cout << *topo3;
			decayTopologyGraphType::nodeIterator iNode, iNodeEnd;
			cout << "!!! topo2 vertex pointers: ";
			for (tie(iNode, iNodeEnd) = topo2.nodes(); iNode != iNodeEnd; ++iNode)
				cout << topo2.vertex(*iNode) << "    ";
			cout << endl;
			cout << "!!! topo3 vertex pointers: ";
			for (tie(iNode, iNodeEnd) = topo3->nodes(); iNode != iNodeEnd; ++iNode)
				cout << topo3->vertex(*iNode) << "    ";
			cout << endl;
			cout << "!!! topo3 interaction vertex pointers: ";
			for (unsigned int i = 0; i < topo3->nmbDecayVertices(); ++i)
				cout << topo3->decayVertices()[i] << "    ";
			cout << endl;
			cout << "!!! topo3 isobar decay vertex pointers: ";
			for (unsigned int i = 0; i < topo3->nmbDecayVertices(); ++i)
				cout << topo3->isobarDecayVertices()[i] << "    ";
			cout << endl;
			decayTopologyGraphType::edgeIterator iEdge, iEdgeEnd;
			cout << "!!! topo2 particle pointers: ";
			for (tie(iEdge, iEdgeEnd) = topo2.edges(); iEdge != iEdgeEnd; ++iEdge)
				cout << topo2.particle(*iEdge) << "    ";
			cout << endl;
			cout << "!!! topo3 particle pointers: ";
			for (tie(iEdge, iEdgeEnd) = topo3->edges(); iEdge != iEdgeEnd; ++iEdge)
				cout << topo3->particle(*iEdge) << "    ";
			cout << endl;
			cout << "!!! topo3 FS particle pointers: ";
			for (unsigned int i = 0; i < topo3->nmbFsParticles(); ++i)
				cout << topo3->fsParticles()[i] << "    ";
			cout << endl;
		}

		{
			cout << "-------------------------------------------------------------------------------"
			     << endl;
			printInfo << "testing subdecays and merge" << endl;
			isobarDecayTopology subDecay1 = topo2.subDecay(2);
			subDecay1.decayTopologyGraphType::printPointers(cout);
			subDecay1.decayTopology::print(cout);
			subDecay1.checkTopology();
			cout << "subdecay from node[2]: " << subDecay1;
			const isobarDecayTopology& subDecay2 = topo2.subDecay(1);
			cout << "subdecay from node[1]: " << subDecay2;
			const isobarDecayTopology& subDecay3 = topo2.subDecay(9);
			cout << "subdecay from node[9]: " << subDecay3;
			subDecay1.addDecay(subDecay3);
			cout << "subdecay from node[2] + subdecay from node[9]: " << subDecay1;
			isobarDecayTopology decay;
			//                                    2I  G  2J  P  C  2M
			particlePtr X2 = createParticle("X-", 2, -1, 4, +1, 0, 2);
			isobarDecayVertexPtr vertX= createIsobarDecayVertex(X2, pi4, f1, 2, 2);
			decay.addVertex(vertX);
			subDecay1.addDecay(decay);
			cout << "subdecay from node[2] + subdecay from node[9] + vert0: " << subDecay1;
			const isobarDecayTopology& subDecayX = topo2.subDecay(2);
			isobarDecayTopology newDecay = isobarDecayTopology::joinDaughterDecays(vertX, subDecayX, subDecay3);
			diffractiveDissVertexPtr prodVertX = createDiffractiveDissVertex(beam, target, X2);
			newDecay.setProductionVertex(prodVertX);
			cout << "joined graph: " << newDecay;
			newDecay.checkTopology();
			newDecay.checkConsistency();
		}
	}
}
