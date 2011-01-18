/*
 * TrpwaFitOptionsFrame.h
 *
 *  Created on: Jan 18, 2011
 *      Author: Promme
 *
 *      Frame to set fit options before sending the fits
 *      to the job manager
 */

#include <TGFrame.h>
#include <vector>
#include <map>
#include <string>
#include "TGButton.h"
//#include "TGComboBox.h"
//#include "TRootEmbeddedCanvas.h"
//#include "TCanvas.h"
//#include "TGDoubleSlider.h"
#include "TGNumberEntry.h"

using namespace std;

#ifndef TRPWAFITOPTIONSFRAME_H_
#define TRPWAFITOPTIONSFRAME_H_

struct Tfit_options{
	unsigned int rank; // rank of fit
	unsigned int niterations; // number of iterations per bin
	bool use_normalization;	// use normalization integral
	int seed; // -1 : random seed else as given
};

class TrpwaFitOptionsFrame : public TGTransientFrame {
private:
	Tfit_options* _fit_options;
	TGNumberEntry* entry_seed; // user number for seed is set here
	TGNumberEntry* entry_nfits; // number of fits per bin set
	TGNumberEntry* entry_rank; // rank of fit
	TGCheckButton* button_use_normalization; // whether to use the normalization integral
	// create all buttons etc.
	void Build();
public:
	TrpwaFitOptionsFrame(Tfit_options& fit_options); // will be set by the user
	virtual ~TrpwaFitOptionsFrame();
	// call the root script for class definition
	ClassDef(TrpwaFitOptionsFrame,0);

	void SetOptions(); // read all options set by the user
	void CallSetOptions(long value); // calls only SetOptions for number entries
};

#endif /* TRPWAFITOPTIONSFRAME_H_ */
