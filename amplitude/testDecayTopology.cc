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
#include "decayGraph.hpp"


using namespace std;
using namespace rpwa;


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
  decayGraph<interactionVertex, particle>::setDebug(true);

  particleDataTable& pdt = particleDataTable::instance();
  pdt.readFile();


  if (1) {
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

    // construct graph
    decayGraph<interactionVertex, particle> g;
    g.name() = "test graph";
    g.addVertex(vert0);
    g.addVertex(vert1);
    g.addVertex(vert2);
    g.addVertex(vert3);
    cout << g << endl;

    for (decayGraph<interactionVertex, particle>::nodeIterator i = g.nodes().first;
	 i != g.nodes().second; ++i) {
      const isobarDecayVertex& v = *static_cast<isobarDecayVertex*>(&g.vertex(*i));
      g.name(v) = v.mother().name();
    }
    for (decayGraph<interactionVertex, particle>::edgeIterator i = g.edges().first;
	 i != g.edges().second; ++i) {
      stringstream n;
      n << "edge " << g.index(*i);
      g.name(*i) = n.str();
    }
    g.print(cout, g.nodeNameMap(), g.edgeNameMap());
    ofstream graphVizFile("decay.dot");
    g.writeGraphViz(graphVizFile, g.nodeNameMap(), g.edgeNameMap());
    gSystem->Exec("dot -Tps -o decay.ps decay.dot");
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
