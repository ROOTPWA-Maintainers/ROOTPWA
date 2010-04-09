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
//      Sebastian Neubert    TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <fstream>

#include "TVector3.h"
#include "TSystem.h"

#include "utilities.h"
#include "particleDataTable.h"
#include "particle.h"
#include "diffractiveDissVertex.h"
#include "diffractiveDissVertex2.h"
#include "isobarDecayVertex.h"
#include "isobarDecayVertex2.h"
#include "isobarDecayTopology.h"
#include "isobarDecayTopology2.h"


using namespace std;
using namespace rpwa;


int
main(int argc, char** argv)
{
  // switch on debug output
  //particleProperties::setDebug(true);
  //particleDataTable::setDebug(true);
  particle::setDebug(true);
  //interactionVertex::setDebug(true);
  //diffractiveDissVertex::setDebug(true);
  //isobarDecayVertex::setDebug(true);
  //decayTopology::setDebug(true);
  isobarDecayTopology::setDebug(true);

  //decayTopologyGraphType::setDebug(true);
  //isobarDecayVertex2::setDebug(true);
  //decayTopology2::setDebug(true);
  isobarDecayTopology2::setDebug(true);

  particleDataTable& pdt = particleDataTable::instance();
  pdt.readFile();

  if (0) {
    particleProperties prop1("bla",2,-1,0,+1,+1);
    particleProperties prop2("blub",2,+1,2,-1,-1);
  
    string opt="IGJPC";
 
    cout << "Comparison  result: " << (prop1 == prop2) << endl;
    cout << "Comparison with opt="<<opt<<"  result: "
	 << (prop1 ==  pair<particleProperties,string>(prop2,opt)) << endl;

    vector<const particleProperties*> selection=pdt.entriesMatching(prop2,opt);
    cout << "Matching entries in pdt"<<" with prototype: "<< endl;
    cout << prop2 << endl;
    cout <<" with option "<<opt<< endl;

    for(unsigned int i=0;i<selection.size();++i){
      cout << *selection[i] << endl;
    }


    // define final state particles
    particle pi0("pi-");
    particle pi1("pi+");
    particle pi2("pi-");
 
    // define isobars
    particle sigma("sigma");
    //               I   G  2J  P   C  2M
    particle X("X-", 2, -1, 0, -1, +1, 0);
  
    particle beam("pi-");
    diffractiveDissVertex prodVert(beam, X);
  
    int lmax=8;

    isobarDecayVertex vert0(sigma, pi0,   pi1, 0, 0);  
    isobarDecayVertex vert1(X,     sigma, pi2, 0, 0);

    //build graph
    vector<particle*> fsParticles;
    fsParticles.push_back(&sigma);
    fsParticles.push_back(&pi2);
    //fsParticles.push_back(&pi2);
  
    vector<isobarDecayVertex*> decayVertices;
    decayVertices.push_back(&vert1);
    //decayVertices.push_back(&vert0);
    isobarDecayTopology topo(fsParticles, prodVert, decayVertices);
    cout << endl;
    printInfo << "decay toplogy:" << endl
	      << topo << endl;
  
    cout << "Number of vertices: " << topo.interactionVertices().size() << endl;
    cout << "Number of vertices: " << topo.isobarDecayVertices().size() << endl;


    fsParticles.clear();
    fsParticles.push_back(&pi0);
    fsParticles.push_back(&pi1);
    decayVertices.clear();
    decayVertices.push_back(&vert0);
    isobarDecayTopology topo2(fsParticles, prodVert, decayVertices);
    cout << endl;
    printInfo << "decay toplogy2:" << endl
	      << topo2 << endl;

    cerr << "Adding..." << endl;
    topo.addSubDecay(topo2);
  
    printInfo << "Joined decay toplogy:" << endl
	      << topo << endl;


    
    const bool topologyOkay = topo.verifyTopology();
    cout << "topology okay = " << topologyOkay << endl << endl;

    printInfo << endl << "Performing consistency checks on Isobar-Vertices.." << endl;
    const bool consistency = topo.checkConsistency();
    cout << "consistency = " << consistency << endl << endl;

   

    vector<isobarDecayVertex*> comb1;
    vert0.getListOfValidDecays(comb1,lmax);

    cout << "Created "<<comb1.size()<<" new decays"<<endl;

    cout << *(comb1[0]) << endl;

    vector<isobarDecayVertex*> dummy;dummy.push_back(NULL);

    vector<isobarDecayVertex*> comb2;
    vert1.getListOfValidDecays(comb1,dummy,comb2,lmax);

    cout << "Created "<<comb2.size()<<" new decays"<<endl;
  }

  if (0) {
    // define final state particles
    particle pi0("pi-");
    particle pi1("pi+");
    particle pi2("pi-");
    particle pi3("pi+");
    particle pi4("pi-");
    // define isobars
    particle i0("isobarA0");
    particle i1("isobarB-");
    particle i2("isobarC0");
    // define X-system
    particle X("X-");
    X.setMass(2.5);
    X.setWidth(0.3);
    // define production vertex
    particle beam("pi-");
    diffractiveDissVertex prodVert(beam, X);
    // define vertices
    isobarDecayVertex vert0(X,  i0, i1);
    isobarDecayVertex vert1(i0, pi0, pi1);
    isobarDecayVertex vert2(i1, pi2, i2);
    isobarDecayVertex vert3(i2, pi3, pi4);
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
    printInfo << "decay topology:" << endl
	      << topo << endl;

    vector<isobarDecayTopology*> decays = topo.possibleDecays();
    unsigned int consistentDecays=0;
    unsigned int inconsistentDecays=0;
    for (unsigned int i = 0; i < decays.size(); ++i) {
      cout << *decays[i];
      decays[i]->verifyTopology();
      isobarDecayVertex::setDebug(true);
      bool isConsistent = decays[i]->checkConsistency();
      isobarDecayVertex::setDebug(false);
      if (isConsistent){
	cout << "isobar decay topology is consistent" << endl;
	++consistentDecays;
      }
      else {
	cout << "isobar decay topology is NOT consistent" << endl;
	++inconsistentDecays;
      }
    }
    cout << "Got "<<inconsistentDecays<<" inconsistent"<<endl
	 << "and "<<consistentDecays<<" valid decays"<<endl
	 << "out of "<<decays.size()<<" constructed decays"<<endl;
  }

  if (1) {
    particlePtr pi0 = createParticle("pi-");
    particlePtr pi1 = createParticle("pi+");
    particlePtr pi2 = createParticle("pi-");
    particlePtr pi3 = createParticle("pi+");
    particlePtr pi4 = createParticle("pi-");
    // define isobars
    particlePtr i0 = createParticle("isobarA0");
    particlePtr i1 = createParticle("isobarB-");
    particlePtr i2 = createParticle("isobarC0");
    // define X-system
    particlePtr X = createParticle("X-");
    X->setMass(2.5);
    X->setWidth(0.3);
    // define production vertex
    particlePtr beam = createParticle("pi-");
    diffractiveDissVertexPtr prodVert = createDiffractiveDissVertex(beam, X);
    // define vertices
    isobarDecayVertexPtr vert0 = createIsobarDecayVertex(X,  i0,  i1);
    isobarDecayVertexPtr vert1 = createIsobarDecayVertex(i0, pi0, pi1);
    isobarDecayVertexPtr vert2 = createIsobarDecayVertex(i1, pi2, i2);
    isobarDecayVertexPtr vert3 = createIsobarDecayVertex(i2, pi3, pi4);
    vector<particlePtr> fsParticles;
    fsParticles.push_back(pi0);
    fsParticles.push_back(pi1);
    fsParticles.push_back(pi2);
    fsParticles.push_back(pi3);
    fsParticles.push_back(pi4);
    vector<isobarDecayVertexPtr> decayVertices;
    decayVertices.push_back(vert3);
    decayVertices.push_back(vert1);
    decayVertices.push_back(vert2);
    decayVertices.push_back(vert0);

    isobarDecayTopology2 topo(prodVert, decayVertices, fsParticles);
    cout << endl;
    printInfo << "decay topology:" << topo;
    vector<isobarDecayTopology2> decays             = topo.possibleDecays();
    unsigned int                 consistentDecays   = 0;
    unsigned int                 inconsistentDecays = 0;
    for (unsigned int i = 0; i < decays.size(); ++i) {
      cout << decays[i];
      // decays[i].printPointers(cout);
      // for (decayTopologyGraphType::nodeIterator j = decays[i].nodes().first;
      // 	   j != decays[i].nodes().second; ++j)
      //   decays[i].vertex(*j)->printPointers(cout);
      //isobarDecayVertex2::setDebug(true);
      bool isConsistent = decays[i].checkTopology() && decays[i].checkConsistency();
      //isobarDecayVertex2::setDebug(false);
      if (isConsistent){
      	cout << "isobar decay topology is consistent" << endl;
      	++consistentDecays;
      }
      else {
      	cout << "isobar decay topology is NOT consistent" << endl;
      	++inconsistentDecays;
      }
    }
    cout << "got " << inconsistentDecays << " inconsistent" << endl
	 << "and " << consistentDecays << " valid decays" << endl
	 << "out of " << decays.size() << " constructed decays" << endl;
  }
}
