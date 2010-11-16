/*
 * TrpwaMainFrame.h
 *
 *  Created on: Aug 24, 2010
 *      Author: Promme
 *
 *      Main Frame for rootpwa gui managing all windows for all steps of
 *      rootpwa
 *
 *      (24.08.2010)
 *      - Implementation of basic layout as singleton
 *      - since only one instance must exist
 *      (15.11.2010)
 *      - Implementation of tree filtering
 *      (16.11.2010)
 *      - Some cosmetics
 */
#include <TGFrame.h>
#include "TRootEmbeddedCanvas.h"
#include <TGProgressBar.h>
#include <vector>
#include "TrpwaSessionManager.h"
#include "TrpwaWaveSelectFrame.h"
#include <TGButton.h>

using namespace std;

#ifndef TRPWAMAINFRAME_H_
#define TRPWAMAINFRAME_H_

class TrpwaMainFrame : public TGMainFrame {
private:
	TRootEmbeddedCanvas *fEcanvas;

	TrpwaMainFrame(const TGWindow *p,UInt_t w,UInt_t h);

	// create all buttons etc.
	void Build();

	virtual ~TrpwaMainFrame();

	static TrpwaMainFrame* pInstance; // singleton reference

	TGVerticalFrame *steplist; // frame holding frames for each step
	TGGroupFrame *frame_session_options;

	vector <TGProgressBar*> progressbars; // references to the progress bars
	vector <TGTextButton*>  stepbuttons;  // references to the buttons calling the step procedures

	TrpwaSessionManager* current_session;

	TrpwaWaveSelectFrame* frame_wave_select; // window for wave selection

public:
	// return the reference to the only instance of this class
	static TrpwaMainFrame& Instance(){
		if (!pInstance){
			pInstance = new TrpwaMainFrame(gClient->GetRoot(),800,800);
		}
		return *pInstance;
	};

	void CheckStatus();

	// save the current session and start a new one
	void NewSession();

	// save the current session
	void SaveSession();

	// save the current session as (opening a dialog)
	void SaveSessionAs();

	// load a session
	void LoadSession();

	// check the status of the session
	void Update();

	void CloseWindow(); // override to call the application abortion

	void SelectWaves(); // calls the interface for wave selection

	void FitPartialWaves(); // send jobs for partial wave fitting

	void ShowFitResults(); // send a command to show the available fit results

	void GenKeys(); // send a command to generate the key files

	void CalcAmps(); // send jobs to calculate amplitudes

	void IntAmps(); // send jobs to integrate available amplitudes

	void FilterData(); // filter data given by root files

	void SetupWorkspace(); // Setup workspace as specified in the init file

	void Dummy(); // will be called for every function not implemented yet

	// call the root script for class definition
	ClassDef(TrpwaMainFrame,0);
};


#endif /* TRPWAMAINFRAME_H_ */
