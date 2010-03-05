///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009
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
//      singleton class that manages all particle data
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "utilities.h"
#include "particleDataTable.h"

	
using namespace std;
using namespace rpwa;


// static member initialization
particleDataTable particleDataTable::_instance;

	






void
particleDataTable::initialize(const char* PDTfile)
{
  string name;
  double mass, width;
  int isospin, gparity, spin, parity,cparity;
  ifstream ifs(PDTfile);
	
  while (!(ifs >> name).eof()) {
    ifs >> mass;
    ifs >> width;
    ifs >> isospin;
    ifs >> gparity;
    ifs >> spin;
    ifs >> parity;
    ifs >> cparity;
    insert(particleData(name, mass, width, isospin, gparity, spin, parity, cparity));
  }
}


void
particleDataTable::insert(const particleData& p)
{
  head = new tableEntry(p, head);
}


void
particleDataTable::print() const
{
  tableEntry* te = this->head;
  while (te != NULL) {
    te->print();
    te = te->next();
  }
}


void
particleDataTable::dump() const
{
  tableEntry* te = this->head;
  while (te != NULL) {
    te->dump();
    te = te->next();
  }
}


int
particleDataTable::ListLen() const
{
  int len = 0;
  tableEntry* te = this->head;
  while (te != NULL) {
    len++;
    te = te->next();
  }
  return(len);
}
	

char**
particleDataTable::List() const
{
  int particleno = 0;
  int listlen = ListLen();
  particleData p;
  char** list;
  list = (char**) malloc(listlen * sizeof(char*));
  tableEntry *te = this->head;
  while (te != NULL) {
    p = te->Particle();
    list[particleno] = (char*) malloc(((p.Name().length()) + 1) * sizeof(char));
    strcpy(list[particleno], p.Name().c_str());
    //     p.Name().to_char_type(list[particleno]);
    particleno++;
    te = te->next();
  }
  return(list);
}
 

particleData
particleDataTable::get(const string& name) const
{
  tableEntry* te = this->head;
  while (te != NULL) {
    if ((te->Particle()).Name() == name)
      return te->Particle();
    te = te->next();
  }
  return particleData(); 
}


double
particleDataTable::mass(const string& name) const
{
  return this->get(name).Mass();
}


double
particleDataTable::width(const string& name) const
{
  return this->get(name).Width();
}
