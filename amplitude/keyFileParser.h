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
//      singleton class that constructs decay topologies according to key files
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef KEYFILEPARSER_H
#define KEYFILEPARSER_H


//#include <iostream>
#include <string>
#include <vector>

#include "libconfig.h++"

#include "isobarDecayTopology.h"


namespace rpwa {


  class keyFileParser {

  public:

    static keyFileParser& instance() { return _instance; }  ///< get singleton instance
  
    isobarDecayTopologyPtr parse(const std::string& keyFileName);

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


  private:

    keyFileParser () { }
    ~keyFileParser() { }
    keyFileParser (const keyFileParser&);
    keyFileParser& operator =(const keyFileParser&);

    const libconfig::Setting* findGroup(const libconfig::Setting& parent,
					const std::string&        groupName,
					const bool                mustExist = true);  ///< finds field in keyfile and makes sure it is a group

    const libconfig::Setting* findList(const libconfig::Setting& parent,
				       const std::string&        listName,
				       const bool                mustExist = true);  ///< finds field in keyfile and makes sure it is a non-empty list

    bool constructParticle(const libconfig::Setting& particleKey,
			   particlePtr&              particle);  ///< creates particle using name in particle key
    
    bool constructDecayVertex(const libconfig::Setting& parentKey,
			      const particlePtr&        parentParticle);  ///< recursively traverses decay chain and creates decay vertices and final state particles

    bool constructProductionVertex(const libconfig::Setting& rootKey,
				   const particlePtr&        X);  ///< creates production vertex


    static keyFileParser                     _instance;       ///< singleton instance

    static interactionVertexPtr              _prodVert;       ///< production vertex
    static std::vector<isobarDecayVertexPtr> _decayVertices;  ///< isobar decay vertices
    static std::vector<particlePtr>          _fsParticles;    ///< final state particles

    static bool _debug;     ///< if set to true, debug messages are printed

  };


}  // namespace rpwa


#endif  // KEYFILEPARSER_H
