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
//      Code file for the GuiFileListTreeModel class that provides a model
//		for QTreeview for a file list view
//
//
// Author List:
//      Stephan Schmeing          TUM            (original author)
//
//
//-------------------------------------------------------------------------

#include <string.h>

#include <QStringList>

#include "reportingUtils.hpp"

#include "GuiFileListTreeModel.h"

using namespace std;
using namespace rpwa;

bool GuiFileListTreeModel::_Debug = false;

// Returns the RootItem for the base class
GuiStringTreeModelItem *GuiFileListTreeModel::RootItem() const{
	return _RootItem;
}

// Returns the RootItem for the base class
GuiStringTreeModelItem *GuiFileListTreeModel::RootItem(){
	return _RootItem;
}

// Adds a child of type Folder to Parent
GuiFileListTreeModelItem *GuiFileListTreeModel::AddFolder(const QString& Data, GuiFileListTreeModelItem *Parent){
	if( Parent ){
		GuiFileListTreeModelItem *Folder = new GuiFileListTreeModelItem(Data, GuiFileListTreeModelItem::Folder, Parent);
		if( Parent->AppendChild(Folder) ){
			if( _Debug ){
				printDebug << "Folder '"<<Data.toLatin1().constData()<<"' added to '"<<Parent->Data().toLatin1().constData()<<"'\n";
			}
			return Folder;
		}
		else{
			delete Folder;
			if( _Debug ){
				printDebug << "GuiFileListTreeModelItem::AddFolder AppendChild failed adding '"<<Data.toLatin1().constData()<<"' to '"<<Parent->Data().toLatin1().constData()<<"'\n";
			}
			return 0;
		}
	}
	else{
		return 0;
	}
}

// Adds a child of type File to Parent
GuiFileListTreeModelItem *GuiFileListTreeModel::AddFile(const QString& Data, GuiFileListTreeModelItem *Parent){
	if( Parent ){
		GuiFileListTreeModelItem *File = new GuiFileListTreeModelItem(Data, GuiFileListTreeModelItem::File, Parent);
		if( Parent->AppendChild(File) ){
			if( _Debug ){
				printDebug << "File '"<<Data.toLatin1().constData()<<"' added to '"<<Parent->Data().toLatin1().constData()<<"'\n";
			}
			return File;
		}
		else{
			delete File;
			if( _Debug ){
				printDebug << "GuiFileListTreeModelItem::AddFile AppendChild failed adding '"<<Data.toLatin1().constData()<<"' to '"<<Parent->Data().toLatin1().constData()<<"'\n";
			}
			return 0;
		}
	}
	else{
		return 0;
	}
}

// Calls GuiStringTreeModel constructor and creates _RootItem
GuiFileListTreeModel::GuiFileListTreeModel( QObject *Parent ):
		GuiStringTreeModel(Parent){
	_RootItem = new GuiFileListTreeModelItem("Master of the Universe", GuiFileListTreeModelItem::Master);
}

// Deletes _RootItem
GuiFileListTreeModel::~GuiFileListTreeModel(){
	delete _RootItem;
}

// Adds a file string to the model
void GuiFileListTreeModel::AddFileString(const QString& Data){
	QStringList DataList = Data.split('/');
	GuiFileListTreeModelItem *ParentItem = _RootItem;
	int ChildIndex = 0;

	// If the first element is emtpy no drive specification is included and therefore a / is added to the second element and the first is deleted
	if( DataList.size() > 1){
		if( DataList[0].isEmpty() ){
			DataList[1].push_front( QChar('/') );
			DataList.removeAt(0);
		}
	}

	// Finding the last folder that is already included in the model
	int i = 0;
	while( -1 != ( ChildIndex = ParentItem->Find(DataList[i]) ) ){
		ParentItem = ParentItem->Child( ChildIndex );
		++i;
	}

	// Adding all Folders that are not already bin added to the model
	for(; i < DataList.size() - 1; ++i ){
		ParentItem = AddFolder(DataList[i], ParentItem);
	}

	// Adding the file
	if( i < DataList.size() ){
		// File has not already bin added to the model
		AddFile(DataList[i], ParentItem);
	}

	if( _Debug ){
		printDebug << "File String added\n";
	}
	// Sends signal so other items can react to it(for example the QTreeView showing the data)
	emit GuiStringTreeModel::layoutChanged();
}

// Adds a complete folder to the model
void GuiFileListTreeModel::AddFolderString(const QDir& Folder){
	QStringList DataList = Folder.absolutePath().split('/');
	GuiFileListTreeModelItem *ParentItem = _RootItem;
	int ChildIndex = 0;

	// If the first element is emtpy no drive specification is included and therefore a / is added to the second element and the first is deleted
	if( DataList.size() > 1){
		if( DataList[0].isEmpty() ){
			DataList[1].push_front( QChar('/') );
			DataList.removeAt(0);
		}
	}

	// Finding the last folder that is already included in the model
	int i = 0;

	while(  ( i < DataList.size() ) && ( -1 != ( ChildIndex = ParentItem->Find(DataList[i]) ) ) ){
		ParentItem = ParentItem->Child( ChildIndex );
		++i;
	}

	// Adding all Folders that are not already bin added to the model
	for(; i < DataList.size(); ++i ){
		ParentItem = AddFolder(DataList[i], ParentItem);
	}

	// Finding all text files in the directory
	QStringList Filters;
	Filters << "*.txt";
	QStringList FileList = Folder.entryList(Filters);

	// Deleting all files from the list that are already included in the model
	i = 0;

	while( i < FileList.size() ){
		if( -1 == ParentItem->Find(FileList[i]) ){
			// File is not already included
			++i;
		}
		else{
			// File is already included
			FileList.removeAt(i);
		}
	}

	// Adding the files from the list to the model
	for(i = 0; i < FileList.size(); ++i ){
		AddFile(FileList[i], ParentItem);
	}

	if( _Debug ){
		printDebug << "File String added\n";
	}
	// Sends signal so other items can react to it(for example the QTreeView showing the data)
	emit GuiStringTreeModel::layoutChanged();
}

///< Deletes an Item and all parents that are left without a child (except the _RootItem of course)
void GuiFileListTreeModel::DeleteItem( const QModelIndex& Item ){
	GuiFileListTreeModelItem *ParentItem = static_cast<GuiFileListTreeModelItem *>( Item.parent().internalPointer() );
	int ItemIndex = Item.row();

	//Remove Item
	ParentItem->DeleteChild( ItemIndex );

	//Remove all parents that are left without a child (except the _RootItem of course)
	while( ( !ParentItem->NChildren() ) && ParentItem->Type() != GuiFileListTreeModelItem::Master ){
		ItemIndex = ParentItem->IndexAtParent();
		ParentItem = ParentItem->Parent();
		ParentItem->DeleteChild( ItemIndex );
	}
}
