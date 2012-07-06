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
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
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
#include <QFileDialog>

#include "GuiPwaMain.h"

using namespace std;
using namespace rpwa;

bool GuiPwaMain::_Debug = false;

// Loads the selected rootfile and tree into _RootDataObject and updates _WaveTreeModel when "Open Tree" is selected from the menu
void GuiPwaMain::on_actionOpenTree_triggered(){
	QString fileName = QFileDialog::getOpenFileName(this, tr("Select Root Tree"), "", tr("Root Files (*.root)"));
	if( !fileName.isEmpty() ){
		_RootDataObject.Clear();
		_WaveTreeModel.Clear();
		_PlotWidget->GetCanvas()->Clear();
		delete _ShownHistogram;
		_ShownHistogram = 0;
		_PlotWidget->Refresh();

		if( _RootDataObject.LoadFile( (fileName.toLatin1()).constData() ) ){
			list<string> TreeList;
			int NumberOfTrees = _RootDataObject.TreesInFile( TreeList );
			if( 0 == NumberOfTrees ){
				QMessageBox::warning(this, tr("Error"), "File does not contain any tree. Please load another file." );
			}
			else if( 1 == NumberOfTrees ){
				_RootDataObject.SelectTree( TreeList.front() );
				_RootDataObject.MapTreeByMassWithHighestLikelihood();
				_WaveTreeModel.AddTree(_RootDataObject);
			}
			else{
				QString NumberStr;
				QMessageBox::information(this, tr("Test"), "File contains "+NumberStr.setNum(NumberOfTrees)+" trees. Tree selection is not implemented yet." );
			}
		}
	}
}


void GuiPwaMain::mon__WaveSelection_currentChanged( const QModelIndex& Current, const QModelIndex & Previous ){
	_PlotWidget->GetCanvas()->Clear();
	delete _ShownHistogram;
	_ShownHistogram = 0;

	if( static_cast<GuiWaveTreeModelItem *>( Current.internalPointer() )->Type() == GuiWaveTreeModelItem::Wave ){
		_ShownHistogram = _RootDataObject.IntensityHist( Current.data().toString().toLatin1().constData(), "MassBin [Gev]" );

		_ShownHistogram->Draw("E1");
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
