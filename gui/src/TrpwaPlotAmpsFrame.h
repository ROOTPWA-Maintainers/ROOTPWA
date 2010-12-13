/*
 * TrpwaPlotAmpsFrame.h
 *
 *  Created on: Dec 13, 2010
 *      Author: Promme
 */
#include <TGFrame.h>
#include <vector>
#include <map>
#include <string>
#include "TGButton.h"
#include "TGComboBox.h"

using namespace std;

#ifndef TRPWAPLOTAMPSFRAME_H_
#define TRPWAPLOTAMPSFRAME_H_

class TrpwaPlotAmpsFrame : public TGTransientFrame {
private:

public:

	// provide a selection of waves to be modified by the user here
	TrpwaPlotAmpsFrame();

	// create all buttons etc.
	void Build();

	virtual ~TrpwaPlotAmpsFrame();

	// call the root script for class definition
	ClassDef(TrpwaPlotAmpsFrame,0);
};

#endif /* TRPWAPLOTAMPSFRAME_H_ */
