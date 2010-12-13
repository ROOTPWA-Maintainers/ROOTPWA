/*
 * TrpwaWaveSelectFrame.cc
 *
 *  Created on: Aug 26, 2010
 *      Author: Promme
 */

#include "TrpwaPlotAmpsFrame.h"
#include <sstream>
#include <iostream>
#include "TGButtonGroup.h"
#include "TGButton.h"
#include "TGTab.h"

using namespace std;

#ifndef TRPWAPLOTAMPSFRAME_CC_
#define TRPWAPLOTAMPSFRAME_CC_

TrpwaPlotAmpsFrame::TrpwaPlotAmpsFrame()
	: TGTransientFrame(gClient->GetRoot(),0,600,600, kHorizontalFrame){
	Build();
	fClient->WaitFor(this); // wait till the user closes this window
}

TrpwaPlotAmpsFrame::~TrpwaPlotAmpsFrame(){
	Cleanup();
}

#include "TrpwaCommonTools.h"

void TrpwaPlotAmpsFrame::Build(){


	SetWindowName("Inspect Results");
	Resize(GetDefaultSize());
	SetWMSize(fWidth, fHeight);
	CenterOnParent();
	MapSubwindows();
	MapWindow();
}


#endif /* TRPWAPLOTAMPSFRAME_CC_ */
