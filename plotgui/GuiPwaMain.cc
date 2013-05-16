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

#include <vector>
#include <list>

#include <QMessageBox>

#include <TLegend.h>

#include "GuiOpenMultipleFilesDialog.h"

#include "GuiPwaMain.h"

using namespace std;
using namespace rpwa;

bool GuiPwaMain::_Debug = false;

// Clears the canvas and deletes the histogram
void GuiPwaMain::ClearHistogram(){
	_PlotWidget->Clear();

	// Deletes all Histograms in list
	while( !_ShownHistograms.empty() ){
		delete _ShownHistograms.back();
		_ShownHistograms.pop_back();
	}
}

// Clears the canvas from the currently loaded waves
void GuiPwaMain::ClearWaves(){
	_RootDataObject.Clear();
	_WaveTreeModel.Clear();
	_PlotWidget->GetCanvas()->Clear();
	ClearHistogram();
	_PlotWidget->Refresh();
}

// Creates and draws the histogram depending on the selected waves and the plotting mode
void GuiPwaMain::DrawHistogram(){
	QModelIndexList selectedWaves = _WaveSelection.selectedIndexes();

	if( selectedWaves.isEmpty() ){
		ClearHistogram();
		_PlotWidget->Refresh();
	}
	else{
		switch( _PlotMode->currentIndex() ){
		case 0: // Single Intensity
		{
			ClearHistogram();

			_ShownHistograms.push_back( _RootDataObject.IntensityHist( selectedWaves[0].data().toString().toLatin1().constData(), "MassBin [Gev]" ) );
			if( _ShownHistograms.back() ){
				_ShownHistograms.back()->Draw("e1");
			}

			for(int i=1; i < selectedWaves.size(); ++i){
				_ShownHistograms.push_back( _RootDataObject.IntensityHist( selectedWaves[i].data().toString().toLatin1().constData(), "MassBin [Gev]" ) );

				if( _ShownHistograms.back() ){
					if( 9 > i ){
						_ShownHistograms.back()->SetMarkerColor(i);
						_ShownHistograms.back()->SetLineColor(i);
					}
					else{ // Color 10 is white/background
						_ShownHistograms.back()->SetMarkerColor(i+1);
						_ShownHistograms.back()->SetLineColor(i+1);
					}

					_ShownHistograms.back()->Draw("e1same");
				}
			}

			if( 1 < selectedWaves.size() ){
				_PlotWidget->GetCanvas()->BuildLegend();
			}

			_PlotWidget->Refresh();

			break;
		}
		case 1: // Coherent Sum
		{
			ClearHistogram();

			list<string> WaveNames;

			for(int i=0; i < selectedWaves.size(); ++i){
				if( _Debug ){
					printDebug << "Added Wave " << selectedWaves[i].data().toString().toLatin1().constData() << '\n';
				}
				WaveNames.push_back( selectedWaves[i].data().toString().toLatin1().constData() );
			}

			_ShownHistograms.push_back( _RootDataObject.CoherentSumHist( WaveNames, "MassBin [Gev]" ) );
			if( _ShownHistograms.back() ){
				_ShownHistograms.back()->Draw("e1");
			}

			_PlotWidget->Refresh();
			break;
		}
		case 2: // Phase Shift
		{
			ClearHistogram();

			// Sorts the list of selected waves to make correlation between the waves in the histogram and position of histogram easier
			qSort( selectedWaves.begin(), selectedWaves.end() );

			// Divides the Pad into subpads for each phase shift
			_PlotWidget->GetCanvas()->Divide(selectedWaves.size(),selectedWaves.size());

			// Loops over all pairs of waves and creates a phase shift histogram for each pair
			int k;
			for(int i=0; i < selectedWaves.size(); ++i){
				for(k=0; k < selectedWaves.size(); ++k){
					if( _Debug ){
						printDebug << "Processing Pad (" << i << ';' << k << ")\n";
					}

					_ShownHistograms.push_back( _RootDataObject.PhaseShiftHist( selectedWaves[i].data().toString().toLatin1().constData(), selectedWaves[k].data().toString().toLatin1().constData(), "MassBin [Gev]" ) );

					if( _ShownHistograms.back() ){
						if( _Debug ){
							printDebug << "Drawing Pad (" << i << ';' << k << ")\n";
						}

						_PlotWidget->cd( i*selectedWaves.size() + k + 1 );

						_ShownHistograms.back()->Draw("e1");
					}
				}
			}

			_PlotWidget->Refresh();

			break;
		}
		default:
			ClearHistogram();
			_PlotWidget->Refresh();
			break;
		}
	}
}

// Takes care of preparing _RootFataObject for plotting
void GuiPwaMain::ProcessOpenedFile(){
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
				QMessageBox::information(this, tr("Not Implemented"), "Tree contains "+NumberStr.setNum(NumberOfBranches)+" branches. Branch selection is not implemented yet." );
			}
		}
	}
	else{
		QString NumberStr;
		QMessageBox::information(this, tr("Not Implemented"), "File contains "+NumberStr.setNum(NumberOfTrees)+" trees. Tree selection is not implemented yet." );
	}
}

// Loads the selected rootfile and tree into _RootDataObject and updates _WaveTreeModel when "Open Tree" is selected from the menu
void GuiPwaMain::on_actionOpenTree_triggered(){
	QString FileName = QFileDialog::getOpenFileName(this, tr("Select Root Tree"), "", tr("Root Files (*.root)"));
	if( !FileName.isEmpty() ){
		ClearWaves();

		if( _RootDataObject.LoadFile( FileName.toLatin1().constData() ) ){
			ProcessOpenedFile();
		}
	}
}

// Parses the selected CompassPWA txt files into a rootfile and displays its tree like if it was opened with on_actionOpenTree_triggered
void GuiPwaMain::on_actionParse_CompassPWA_txts_triggered(){
	GuiOpenMultipleFilesDialog* CompassFilesDialog = new GuiOpenMultipleFilesDialog( this, tr("Select CompassPWA text files (Fitresults and Integrals)"), "", tr("Text Files(*.txt)") );

	if( CompassFilesDialog->exec() ){
		QString RootFileName = QFileDialog::getSaveFileName(this, tr("Save Parsed Fitresults As"), "", tr("Root Files (*.root)"));
		if( !RootFileName.isEmpty() ){
			ClearWaves();

			vector<string> DataFiles;

			// Parsing the CompassPwa files
			if( _RootDataObject.ParseFiles( "/nfs/hicran/project/compass/analysis/sschmeing/rootpwa/amplitude/particleDataTable.txt", RootFileName.toLatin1().constData(), CompassFilesDialog->GetFiles(DataFiles) ) ){
				ProcessOpenedFile();
			}

			if( _Debug ){
				printDebug << "Received FileList of size: " << DataFiles.size() << '\n';
			}
		}
	}

	delete CompassFilesDialog;
}

void GuiPwaMain::on__PlotMode_currentIndexChanged(){
	DrawHistogram();
}

void GuiPwaMain::mon__WaveSelection_selectionChanged( const QItemSelection & selected, const QItemSelection & deselected ){
	const QModelIndexList& selectedItems( selected.indexes() );
	int allWaves = true;
	QItemSelection unselect;

	for( int i=0; i < selectedItems.size(); i++ ){
		if( static_cast<GuiWaveTreeModelItem *>( selectedItems[i].internalPointer() )->Type() != GuiWaveTreeModelItem::Wave ){
			allWaves = false;
			unselect.select(selectedItems[i],selectedItems[i]);
		}
	}

	if( allWaves ){
		DrawHistogram();
	}
	else{
		_WaveSelection.select(unselect, QItemSelectionModel::Deselect);

	}
}

// Initializes the Gui
GuiPwaMain::GuiPwaMain(QMainWindow *parent):
		QMainWindow(parent),
		_WaveTreeModel(this),
		_WaveSelection(&_WaveTreeModel){
//	RootPwaDataObject::SetDebug(true);
//	GuiPwaMain::SetDebug(true);
//	GuiFileListTreeModel::SetDebug(true);
//	GuiStringTreeModelItem::SetDebug(true);
	setupUi( this );

	 //WaveTreeView defined in ui_GuiPwaMain.h
	_WaveTreeView->header()->hide();
	_WaveTreeView->setModel(&_WaveTreeModel);
	_WaveTreeView->setSelectionMode(QAbstractItemView::ExtendedSelection);
	_WaveTreeView->setSelectionModel(&_WaveSelection);

	// Make the the embedded TCanvas to be the current ROOT TCanvas
	_PlotWidget->Refresh();
	_PlotWidget->cd();
	_PlotWidget->GetCanvas()->SetFillColor(10);

	// Set up connections
	QObject::connect( &_WaveSelection, SIGNAL( selectionChanged( const QItemSelection &, const QItemSelection & ) ), this,  SLOT( mon__WaveSelection_selectionChanged( const QItemSelection &, const QItemSelection & ) ) );
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
