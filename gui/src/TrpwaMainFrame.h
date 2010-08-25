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
 */
#include <TGFrame.h>
#include "TRootEmbeddedCanvas.h"
#include <TGProgressBar.h>
#include <vector>
#include "TrpwaSessionManager.h"

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

	TrpwaSessionManager* current_session;

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

	// load a session
	void LoadSession();

	// check the status of the session
	void Update();

	void CloseWindow(); // override to call the application abortion

	// call the root script for class definition
	ClassDef(TrpwaMainFrame,0);
};


#endif /* TRPWAMAINFRAME_H_ */
