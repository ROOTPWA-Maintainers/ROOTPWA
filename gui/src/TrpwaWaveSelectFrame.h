/*
 * TrpwaWaveSelectFrame.h
 *
 *  Created on: Aug 26, 2010
 *      Author: Promme
 */
#include <TGFrame.h>
#include <vector>
#include <string>

using namespace std;

#ifndef TRPWAWAVESELECTFRAME_H_
#define TRPWAWAVESELECTFRAME_H_

// struct to hold the wave selection per bin
struct TWaveSelection {
	vector<string> available_waves; // all amplitudes that are available to be fit
	vector<bool> selected_waves;    // Users selection for this bin
	int bin_low;                    // lower bound for this bin
	int bin_high;                   // upper bound for this bin
};

typedef vector<TWaveSelection> TWaveSelections;

class TrpwaWaveSelectFrame : public TGTransientFrame {
private:

	TWaveSelections* _waveselections; // set and selected waves

public:

	// provide a selection of waves to be modified by the user here
	TrpwaWaveSelectFrame(TWaveSelections& waveselections);

	// create all buttons etc.
	void Build();

	virtual ~TrpwaWaveSelectFrame();

	// call the root script for class definition
	ClassDef(TrpwaWaveSelectFrame,0);
};

#endif /* TRPWAWAVESELECTFRAME_H_ */
