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
//      Code file for the GuiOpenMultipleFilesDialog class that provides
//		a dialog form where several files can be selected one after
//		another in different directories
//
//
// Author List:
//      Stephan Schmeing          TUM            (original author)
//
//
//-------------------------------------------------------------------------

#include <QMessageBox>

#include "GuiOpenMultipleFilesDialog.h"

using namespace std;
using namespace rpwa;

bool GuiOpenMultipleFilesDialog::GuiSelectedFilesWidgetAdd::_Debug = false;
bool GuiOpenMultipleFilesDialog::_Debug = false;

// Sets up the Widget for the selection view, which is added to the file dialog
GuiOpenMultipleFilesDialog::GuiSelectedFilesWidgetAdd::GuiSelectedFilesWidgetAdd( QWidget * parent ):
		QWidget( parent ),
		_SelectedFilesModel(this),
		_SelectedFilesSelection(&_SelectedFilesModel){
	setupUi( this );

	 //WaveTreeView defined in ui_GuiSelectedFilesWidgetAdd.h
	_SelectedFilesTreeView->header()->hide();
	_SelectedFilesTreeView->setModel(&_SelectedFilesModel);
	_SelectedFilesTreeView->setSelectionModel(&_SelectedFilesSelection);
	_SelectedFilesTreeView->setSelectionMode(QAbstractItemView::MultiSelection);
}

// Returns a list of all selected files
QModelIndexList GuiOpenMultipleFilesDialog::GuiSelectedFilesWidgetAdd::SelectedFiles() const{
	return _SelectedFilesSelection.selectedIndexes();
}

// Adds a file string to the model
void GuiOpenMultipleFilesDialog::GuiSelectedFilesWidgetAdd::AddFileString(const QString& Data){
	_SelectedFilesModel.AddFileString( Data );
}

// Adds a complete folder to the model
void GuiOpenMultipleFilesDialog::GuiSelectedFilesWidgetAdd::AddFolderString(const QDir& Folder){
	_SelectedFilesModel.AddFolderString( Folder );
}

// Deletes an item and all parents that are left without a child (except the _RootItem of course)
void GuiOpenMultipleFilesDialog::GuiSelectedFilesWidgetAdd::DeleteItem( const QModelIndex& Item ){
	_SelectedFilesModel.DeleteItem( Item );
}

// Deletes all selected items and all parents that are left without a child (except the _RootItem of course)
void GuiOpenMultipleFilesDialog::GuiSelectedFilesWidgetAdd::DeleteSelectedItems(){
	QModelIndexList DeletionList = _SelectedFilesSelection.selectedIndexes();

	for( int i = 0; i < DeletionList.size(); i++){
		DeleteItem( DeletionList[i] );
	}
}

void GuiOpenMultipleFilesDialog::mon__BtAdd_clicked(){
	QString SelectionType;

	QModelIndexList SelectionIndexList = _FileDialogListView->selectionModel()->selectedIndexes();
	for( int i=0; i < SelectionIndexList.size(); ++i ){
		if( 0 == SelectionIndexList.at(i).column() ){
			SelectionType = ( _FileDialogListView->model()->index( SelectionIndexList.at(i).row(), 2, SelectionIndexList.at(i).parent() ).data() ).toString();
			if( "txt File" == SelectionType ){
				_SelectedFilesWidget.AddFileString( directory().absoluteFilePath( ( SelectionIndexList.at(i).data() ).toString() ) );
			}
			else if ( "Folder" == SelectionType ){
				_SelectedFilesWidget.AddFolderString( QDir( directory().absoluteFilePath( ( SelectionIndexList.at(i).data() ).toString() ) ) );
			}
			else{
				QMessageBox::warning(this, tr("Error"), "Selection is not of correct file type." );
			}
		}
	}
}

void GuiOpenMultipleFilesDialog::mon__BtAddAll_clicked(){
	_SelectedFilesWidget.AddFolderString( directory() );
}

void GuiOpenMultipleFilesDialog::mon__BtRemove_clicked(){
	_SelectedFilesWidget.DeleteSelectedItems();
}

// Initializes the file dialog
GuiOpenMultipleFilesDialog::GuiOpenMultipleFilesDialog( QWidget * parent, const QString & caption, const QString & directory, const QString & filter ):
		QFileDialog( parent, caption, directory, filter ){
	// Sets standard type for file dialog
	this->setFileMode(QFileDialog::ExistingFiles);

	// Adds the GuiSelectedFilesWidgetAdd widget to the file dialog
	QGridLayout *Layout = static_cast<QGridLayout *>( this->layout() );
	Layout->addWidget( &_SelectedFilesWidget, 0, 3, 3, 3 );

	// Enlarges the dialog to compensate the additional widget
	this->resize( this->width() + _SelectedFilesWidget.width(), this->height() );

	// Get the QListView from the FileDialog
	if( !(_FileDialogListView = this->findChild<QListView*>("listView")) ){
		QMessageBox::warning(this, tr("Error"), "Could not find list view of the file dialog. Buttons deactivated." );
		_SelectedFilesWidget._BtAdd->setEnabled(false);
		_SelectedFilesWidget._BtAddAll->setEnabled(false);
		_SelectedFilesWidget._BtRemove->setEnabled(false);
	}

	// Connects the additional buttons to functions
	QObject::connect( _SelectedFilesWidget._BtAdd, SIGNAL( clicked( bool ) ), this, SLOT( mon__BtAdd_clicked() ) );
	QObject::connect( _SelectedFilesWidget._BtAddAll, SIGNAL( clicked( bool ) ), this, SLOT( mon__BtAddAll_clicked() ) );
	QObject::connect( _SelectedFilesWidget._BtRemove, SIGNAL( clicked( bool ) ), this, SLOT( mon__BtRemove_clicked() ) );
}
