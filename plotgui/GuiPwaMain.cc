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
//      Code file for the GuiPwaMain class that provides
//		the main form for the pwa gui
//
//
// Author List:
//      Stephan Schmeing          TUM            (original author)
//
//
//-------------------------------------------------------------------------

#include <list>

#include <QMessageBox>

#include "GuiOpenMultipleFilesDialog.h"

#include "GuiPwaMain.h"

using namespace std;
using namespace rpwa;

bool GuiPwaMain::_Debug = false;

///< Clears the canvas from the currently loaded waves
void GuiPwaMain::ClearWaves(){
	_RootDataObject.Clear();
	_WaveTreeModel.Clear();
	_PlotWidget->GetCanvas()->Clear();
	delete _ShownHistogram;
	_ShownHistogram = 0;
	_PlotWidget->Refresh();
}

// Loads the selected rootfile and tree into _RootDataObject and updates _WaveTreeModel when "Open Tree" is selected from the menu
void GuiPwaMain::on_actionOpenTree_triggered(){
	QString FileName = QFileDialog::getOpenFileName(this, tr("Select Root Tree"), "", tr("Root Files (*.root)"));
	if( !FileName.isEmpty() ){
		ClearWaves();

		if( _RootDataObject.LoadFile( (FileName.toLatin1()).constData() ) ){
			list<string> TreeList;
			int NumberOfTrees = _RootDataObject.TreesInFile( TreeList );
			if( 0 == NumberOfTrees ){
				QMessageBox::warning(this, tr("Error"), "File does not contain any tree. Please load another file." );
			}
			else if( 1 == NumberOfTrees ){
				if( _RootDataObject.SelectTree( TreeList.front() ) ){
					list<string> BranchList;
					int NumberOfBranches = _RootDataObject.BranchesInTree( BranchList );
					if( 0 == NumberOfBranches ){
						QMessageBox::warning(this, tr("Error"), "Tree does not contain branches of class fitResult." );
					}
					else if( 1 == NumberOfBranches ){
						if( _RootDataObject.SelectBranch( BranchList.front() ) ){
							_RootDataObject.MapTreeByMassWithHighestLikelihood();
							_WaveTreeModel.AddTree(_RootDataObject);
						}
					}
					else{
						QString NumberStr;
						QMessageBox::information(this, tr("Test"), "Tree contains "+NumberStr.setNum(NumberOfBranches)+" branches. Branch selection is not implemented yet." );
					}
				}
			}
			else{
				QString NumberStr;
				QMessageBox::information(this, tr("Test"), "File contains "+NumberStr.setNum(NumberOfTrees)+" trees. Tree selection is not implemented yet." );
			}
		}
	}
}

// Parses the selected CompassPWA txt files into a rootfile and displays its tree like if it was opened with on_actionOpenTree_triggered
void GuiPwaMain::on_actionParse_CompassPWA_txts_triggered(){
	GuiOpenMultipleFilesDialog* CompassFilesDialog = new GuiOpenMultipleFilesDialog( this, tr("Select CompassPWA text files (Fitresults and Integrals)"), "", tr("Text Files(*.txt)") );

	if( CompassFilesDialog->exec() ){
		QString RootFileName = QFileDialog::getSaveFileName(this, tr("Save Parsed Fitresults As"), "", tr("Root Files (*.root)"));
		if( !RootFileName.isEmpty() ){
			// Parsing the CompassPwa files

			// Loading parsed file
		}
	}

	delete CompassFilesDialog;
}

void GuiPwaMain::mon__WaveSelection_currentChanged( const QModelIndex& Current, const QModelIndex & Previous ){
	_PlotWidget->GetCanvas()->Clear();
	delete _ShownHistogram;
	_ShownHistogram = 0;

	if( static_cast<GuiWaveTreeModelItem *>( Current.internalPointer() )->Type() == GuiWaveTreeModelItem::Wave ){
		_ShownHistogram = _RootDataObject.IntensityHist( Current.data().toString().toLatin1().constData(), "MassBin [Gev]" );

		if( _ShownHistogram ){
			_ShownHistogram->Draw("e1");
		}
	}

	_PlotWidget->Refresh();
}
/*void GuiPwaMain::on__WaveTreeView_clicked( const QModelIndex& Index ){
	_PlotWidget->GetCanvas()->Clear();
	delete _ShownHistogram;
	_ShownHistogram = 0;

	if( static_cast<GuiWaveTreeModelItem *>( Index.internalPointer() )->Type() == GuiWaveTreeModelItem::Wave ){
		_ShownHistogram = _RootDataObject.IntensityHist( Index.data().toString().toLatin1().constData(), "MassBin [Gev]" );

		_ShownHistogram->Draw("E1");
	}

	_PlotWidget->Refresh();
}*/

// Initializes the Gui
GuiPwaMain::GuiPwaMain(QMainWindow *parent):
		QMainWindow(parent),
		_WaveTreeModel(this),
		_WaveSelection(&_WaveTreeModel),
		_ShownHistogram(0){
//	RootPwaDataObject::SetDebug(true);
	setupUi( this );

	 //WaveTreeView defined in ui_GuiPwaMain.h
	_WaveTreeView->header()->hide();
	_WaveTreeView->setModel(&_WaveTreeModel);
	_WaveTreeView->setSelectionModel(&_WaveSelection);

	// Make the the embedded TCanvas to be the current ROOT TCanvas
	_PlotWidget->Refresh();
	_PlotWidget->GetCanvas()->cd();

	// Set up connections
	QObject::connect( &_WaveSelection, SIGNAL( currentChanged( const QModelIndex &, const QModelIndex & ) ), this,  SLOT( mon__WaveSelection_currentChanged( const QModelIndex &, const QModelIndex & ) ) );
}

// Event handling
void GuiPwaMain::changeEvent(QEvent *event){
	QWidget::changeEvent(event);
	switch( event->type() ){
	case QEvent::LanguageChange:
		retranslateUi( this );
		break;
	default:
		break;
    }
 }
