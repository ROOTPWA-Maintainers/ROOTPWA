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
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
// File and Version Information:
// $Rev:: 862                         $: revision of last commit
// $Author:: schmeing                 $: author of last commit
// $Date:: 2012-07-06 13:54:31 +0200 #$: date of last commit
//
// Description:
//      Header file for the GuiPwaMain class that provides
//		the main form for the pwa gui
//
//
// Author List:
//      Stephan Schmeing          TUM            (original author)
//
//
//-------------------------------------------------------------------------

#ifndef GuiPwaMain_H
#define GuiPwaMain_H

#include <QObject>
#include <QModelIndex>
#include <QtGui/QMainWindow>

#include <TH1F.h>
#include <TH2F.h>
#include <TQtWidget.h>
#include <TGraph.h>

#include <ui_GuiPwaMain.h>

#include "RootPwaDataObject.h"
#include "GuiWaveTreeModel.h"

namespace rpwa{

	class GuiPwaMain: public QMainWindow, private Ui::GuiPwaMain{
	Q_OBJECT
	private:
		// Variables
		RootPwaDataObject _RootDataObject;
		GuiWaveTreeModel _WaveTreeModel;
		QItemSelectionModel _WaveSelection;
		TH1F * _ShownHistogram;

		static bool _Debug; ///< if set to true, debug messages are printed

		// Functions
		void  ClearWaves(); ///< Clears the canvas from the currently loaded waves

	private slots:
		// Automatically connected slots
		void on_actionOpenTree_triggered(); ///< Loads the selected rootfile and tree into _RootDataObject and updates _WaveTreeModel when "Open Tree" is selected from the menu
		void on_actionParse_CompassPWA_txts_triggered(); ///< Parses the selected CompassPWA txt files into a rootfile and displays its tree like if it was opened with on_actionOpenTree_triggered

		// Manually connected slots
		void mon__WaveSelection_currentChanged( const QModelIndex & Current, const QModelIndex & Previous );

	protected:
		// Variables

		// Functions

	public:
		// Constructors + Destructors
		GuiPwaMain(QMainWindow *parent = 0); ///< Initializes the Gui

		// Get && Set

		// Functions
		void changeEvent(QEvent *event); ///< Event handling
	};

} // namespace rpwa

#endif /* GuiPwaMain_H */
