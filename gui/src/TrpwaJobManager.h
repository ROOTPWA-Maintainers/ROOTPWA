/*
 * TrpwaJobManager.h
 *
 *  Created on: Oct 12, 2010
 *      Author: Promme
 *
 *      class to manages jobs on a cluster farm
 *
 *      (12.10.10)
 *      - Declaration of singleton class
 */

#include <set>
#include <string>

using namespace std;

#ifndef TRPWAJOBMANAGER_H_
#define TRPWAJOBMANAGER_H_

// these are the available farm types
// that might differ in job handling
static const string arrfarmtypes [] = {
		"local ",
		"cern ",
		"gridka ",
		"mainz ",
		"munich "
};

static const set<string> farmtypes(arrfarmtypes, arrfarmtypes+5);

class TrpwaJobManager{
public:
	// return the instance to this object
	static TrpwaJobManager* Instance();

	// return the available farm type
	// (local is always available)
	string GetFarmType(){return _available_farmtype;};

private:
	string _available_farmtype;
	static TrpwaJobManager* _pinstance;

	// farm type is checked on creation
	TrpwaJobManager();
	virtual ~TrpwaJobManager();

	// searches for signatures showing the
	// type of farm
	void CheckFarmType();
};

#endif /* TRPWAJOBMANAGER_H_ */
