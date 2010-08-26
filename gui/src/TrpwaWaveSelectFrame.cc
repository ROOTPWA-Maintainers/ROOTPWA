/*
 * TrpwaWaveSelectFrame.cc
 *
 *  Created on: Aug 26, 2010
 *      Author: Promme
 */

#include "TrpwaWaveSelectFrame.h"

using namespace std;

#ifndef TRPWAWAVESELECTFRAME_CC_
#define TRPWAWAVESELECTFRAME_CC_

TrpwaWaveSelectFrame::TrpwaWaveSelectFrame(TWaveSelections& waveselections)
	: TGTransientFrame(gClient->GetRoot(),0,600,600, kTransientFrame){
	_waveselections = &waveselections; // give a reference to be kept after window distruction
	Build();
	fClient->WaitFor(this); // wait till the user closes this window
}

TrpwaWaveSelectFrame::~TrpwaWaveSelectFrame(){
	Cleanup();
}

void TrpwaWaveSelectFrame::Build(){
	SetWindowName("Wave selection");
	SetWMSize(fWidth, fHeight);
	CenterOnParent();
	MapWindow();
}


#endif /* TRPWAWAVESELECTFRAME_CC_ */
