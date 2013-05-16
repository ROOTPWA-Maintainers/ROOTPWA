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

#include <string>
#include <set>

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

///< Adds all folders in the given path in Data if they do not already exists and returns the GuiFileListTreeModelItem of the last folder in the string
GuiFileListTreeModelItem *GuiFileListTreeModel::AddCompletePath(const QString& Data){
	QStringList DataList = Data.split('/');
	GuiFileListTreeModelItem *ParentItem = _RootItem;
	int ChildIndex = 0;

	// If the first element is empty no drive specification is included and therefore a / is added to the second element and the first is deleted
	if( DataList.size() > 1){
		if( DataList[0].isEmpty() ){
			DataList[1].push_front( QChar('/') );
			DataList.removeAt(0);
		}
	}

	// Finding the last folder that is already included in the model
	int i = 0;
	int n;
	while( ( i < DataList.size() ) && ( -1 != ( ChildIndex = ParentItem->FindSub(DataList[i]) ) ) ){
		ParentItem = ParentItem->Child( ChildIndex );

		// FollowSub takes care of incrementing i to the starting value of the next Item
		if( 0 <= ( n = ParentItem->FollowSub(DataList, i) ) ){
			// Not the complete path from the new ParentItem is included in DataList and therefore ParentItem has to be splitted
			ParentItem->Split(n);

			break;
		}
	}

	// Adding all Folders that are not already bin added to the model
	if( i < DataList.size() ){
		QStringList NewDataList = DataList.mid(i);
		ParentItem = AddFolder(NewDataList.join("/"), ParentItem);
	}

	return ParentItem;
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
	int index = Data.lastIndexOf('/');
	QString Folder = Data.left(index);
	QString File = Data.mid(index+1); // index + 1 to be after the '/'

	GuiFileListTreeModelItem *ParentItem = AddCompletePath(Folder);

	// Adding the file
	if( 0 > ParentItem->Find(File) ){
		// File has not already been added to the model
		AddFile(File, ParentItem);
	}

	if( _Debug ){
		printDebug << "File String added\n";
	}

	// Sends signal so other items can react to it(for example the QTreeView showing the data)
	emit GuiStringTreeModel::layoutChanged();
}

// Adds a complete folder to the model
void GuiFileListTreeModel::AddFolderString(const QDir& Folder){
	GuiFileListTreeModelItem *ParentItem = AddCompletePath( Folder.absolutePath() );

	// Finding all text files in the directory
	QStringList Filters;
	Filters << "*.txt";
	QStringList FileList = Folder.entryList(Filters);

	// Deleting all files from the list that are already included in the model
	int i = 0;

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

///< Deletes all items on DeletionList and all parents that are left without a child (except the _RootItem of course)
void GuiFileListTreeModel::DeleteItems( QModelIndexList DeletionList ){
	if(_Debug){
		printDebug << Print(cout);
	}

	;

	if( _Debug ){
		printDebug << "DeletionList:\n";
		for( int i=0; i < DeletionList.size(); ++i ){
			printDebug << "Type-" << static_cast<GuiFileListTreeModelItem *>( DeletionList[i].internalPointer() )->Type() << "-item: "<<static_cast<GuiFileListTreeModelItem *>( DeletionList[i].internalPointer() )->Data().toLatin1().constData() <<'\n';
		}
		printDebug << endl;
	}

	// Remove all elements from DeletionList for which one of its parents is also in the list (First creating a set of all children of all elements in the list and then remove all elements from the list that are in the set)
	set<GuiStringTreeModelItem const*> ChildrenSet;
	for( QList<QModelIndex>::const_iterator it = DeletionList.constBegin(); it != DeletionList.constEnd(); ++it ){
		static_cast<GuiFileListTreeModelItem *>( it->internalPointer() )->AddChildrenToSet( ChildrenSet );
		if( _Debug ) printDebug << "Loop1\n";
	}
	if( _Debug ){
		printDebug << "ChildrenSet: (Size " << ChildrenSet.size() << ")\n";
		for( set<GuiStringTreeModelItem const*>::const_iterator sit = ChildrenSet.begin(); sit != ChildrenSet.end(); ++sit ){
			printDebug << (*sit) << '\n';
		}
		printDebug << endl;
	}
	for( QList<QModelIndex>::iterator it = DeletionList.begin(); it != DeletionList.end(); ){
		if( ChildrenSet.count( static_cast<GuiFileListTreeModelItem *>(it->internalPointer()) ) ){
			// Object is a child of another object in the list and has to be removed
			if( _Debug ){
				printDebug << "Erase Item: " << it->internalPointer() << '\n';
			}
			DeletionList.erase(it);
		}
		else{
			// Object is not a child
			if( _Debug ){
				printDebug << "Leave Item: " << it->internalPointer() << '\n';
			}
			++it;
		}
	}

	if( _Debug ){
		printDebug << "Reduced DeletionList:\n";
		for( int i=0; i < DeletionList.size(); ++i ){
			printDebug << "Type-" << static_cast<GuiFileListTreeModelItem *>( DeletionList[i].internalPointer() )->Type() << "-item: "<<static_cast<GuiFileListTreeModelItem *>( DeletionList[i].internalPointer() )->Data().toLatin1().constData()<<'\n';
		}
		printDebug << endl;
	}

	QModelIndex ParentItemModelIndex;
	GuiFileListTreeModelItem *ParentItem;
	while( !DeletionList.isEmpty() ){
		ParentItemModelIndex = DeletionList.last().parent();
		ParentItem = static_cast<GuiFileListTreeModelItem *>( DeletionList.last().internalPointer() )->Parent(); // ParentItemModelIndex.parent() does not work since the master object has no ModelIndex and cannot be accessed this way
		int ItemIndex = DeletionList.last().row();

		//Remove Item
		if( _Debug ){
			printDebug << "Remove Item: " << static_cast<GuiFileListTreeModelItem *>( DeletionList.last().internalPointer() )->Data().toLatin1().constData() << endl;
			printDebug << "Parent type: " << ParentItem->Type() << endl;
			printDebug << "Parent data: " << ParentItem->Data().toLatin1().constData() << endl;
			printDebug << "Item index: " << ItemIndex << endl;
		}

		emit beginRemoveRows( ParentItemModelIndex, ItemIndex, ItemIndex );
		ParentItem->DeleteChild( ItemIndex );
		emit endRemoveRows();
		DeletionList.removeLast();

		if( _Debug ){
			printDebug << "Item removed!" << endl;
		}

		//Remove all parents that are left without a child (except the _RootItem of course) and remove them from DeletionList in case they are on it
		while( ( !ParentItem->NChildren() ) && ParentItem->Type() != GuiFileListTreeModelItem::Master ){
			if( _Debug ){
				printDebug << "Remove ParentItem: " << ParentItem->Data().toLatin1().constData() << endl;
			}

			ItemIndex = ParentItemModelIndex.row();
			ParentItem = ParentItem->Parent();
			ParentItemModelIndex = ParentItemModelIndex.parent();

			if( _Debug ){
				printDebug << "Parent type: " << ParentItem->Type() << endl;
				printDebug << "Parent data: " << ParentItem->Data().toLatin1().constData() << endl;
				printDebug << "Item index: " << ItemIndex << endl;
			}

			emit beginRemoveRows( ParentItemModelIndex, ItemIndex, ItemIndex );
			ParentItem->DeleteChild( ItemIndex );
			emit endRemoveRows();

			if( _Debug ){
				printDebug << "ParentItem removed!" << endl;
			}
		}
	}
}

// Appends all files from list to the end of FileList and returns it
vector<string>& GuiFileListTreeModel::GetFiles( vector<string>& FileList ) const{
	return _RootItem->GetFiles( FileList );
}

// Prints the whole Tree
ostream& GuiFileListTreeModel::Print(ostream& out) const{
	QModelIndexList TreeList = ListOfTree();
	out << "Tree:\n";
	for(int i=0; i<TreeList.size(); ++i){
		out << '(' << static_cast<GuiFileListTreeModelItem *>( TreeList[i].internalPointer() )->Type() << ") " << static_cast<GuiFileListTreeModelItem *>( TreeList[i].internalPointer() )->Data().toLatin1().constData() << '\n';
	}

	return out;
}
