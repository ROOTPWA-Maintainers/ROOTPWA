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
// $Rev:: 168                         $: revision of last commit
// $Author:: bgrube                   $: author of last commit
// $Date:: 2010-03-12 22:45:16 +0100 #$: date of last commit
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

#include <boost/config.hpp>
#include <boost/utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/connected_components.hpp>

#include "TVector3.h"
#include "TSystem.h"

#include "utilities.h"
#include "particleDataTable.h"
#include "particle.h"
#include "diffractiveDissVertex.h"
#include "isobarDecayVertex.h"
#include "isobarDecayTopology.h"


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
  // switch on debug output
  //particle::setDebug(true);
  //interactionVertex::setDebug(true);
  diffractiveDissVertex::setDebug(true);
  //isobarDecayVertex::setDebug(true);
  //decayTopology::setDebug(true);
  //isobarDecayTopology::setDebug(true);

  particleDataTable& pdt = particleDataTable::instance();
  pdt.readFile();

  // define final state particles
  particle pi0("pi-");
  particle pi1("pi+");
  particle pi2("pi-");
 
  // define isobars
  particle sigma("sigma");
  //              I   G  2J  P   C  2M
  particle X("X-", 2, -1, 0, -1, +1, 0);
  
  particle beam("pi-");
  diffractiveDissVertex prodVert(beam, X);
  
  int lmax=8;

  isobarDecayVertex vert0(sigma,     pi0, pi1, 0 , 0);  
  isobarDecayVertex vert1(X,         sigma, pi2, 0, 0);

// build graph
    vector<particle*> fsParticles;
    fsParticles.push_back(&pi0);
    fsParticles.push_back(&pi1);
    fsParticles.push_back(&pi2);
  
    vector<isobarDecayVertex*> decayVertices;
    decayVertices.push_back(&vert1);
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

   

    vector<isobarDecayVertex*> comb1;
    vert0.getListOfValidDecays(comb1,lmax);

    cout << "Created "<<comb1.size()<<" new decays"<<endl;

    cout << *(comb1[0]) << endl;

    vector<isobarDecayVertex*> dummy;dummy.push_back(NULL);

    vector<isobarDecayVertex*> comb2;
    vert1.getListOfValidDecays(comb1,dummy,comb2,lmax);

    cout << "Created "<<comb2.size()<<" new decays"<<endl;
    
    

}
