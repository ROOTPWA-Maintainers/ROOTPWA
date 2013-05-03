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
//      Code file for the GuiWaveTreeModel class that provides a model
//		for QTreeview for the wave view in GuiPwaMain
//
//
// Author List:
//      Stephan Schmeing          TUM            (original author)
//
//
//-------------------------------------------------------------------------

#include <string.h>

#include "reportingUtils.hpp"

#include "GuiWaveTreeModel.h"

using namespace std;
using namespace rpwa;

bool GuiWaveTreeModel::_Debug = false;

// Adds a child of type File to _RootItem
GuiWaveTreeModelItem *GuiWaveTreeModel::AddFile(QString Data){
	GuiWaveTreeModelItem *File = new GuiWaveTreeModelItem(Data, GuiWaveTreeModelItem::File, _RootItem);
	if( _RootItem->AppendChild(File) ){
		if( _Debug ){
			printDebug << "File '"<<Data.toLatin1().constData()<<"' added\n";
		}
		return File;
	}
	else{
		delete File;
		if( _Debug ){
			printDebug << "GuiWaveTreeModel::AddFile AppendChild failed\n";
		}
		return 0;
	}
}

// Adds a child of type Tree to Parent
GuiWaveTreeModelItem *GuiWaveTreeModel::AddTree(QString Data, GuiWaveTreeModelItem *Parent){
	if( Parent ){
		GuiWaveTreeModelItem *Tree = new GuiWaveTreeModelItem(Data, GuiWaveTreeModelItem::Tree, Parent);
		if( Parent->AppendChild(Tree) ){
			if( _Debug ){
				printDebug << "Tree '"<<Data.toLatin1().constData()<<"' added to '"<<Parent->Data().toLatin1().constData()<<"'\n";
			}
			return Tree;
		}
		else{
			delete Tree;
			if( _Debug ){
				printDebug << "GuiWaveTreeModel::AddTree AppendChild failed adding '"<<Data.toLatin1().constData()<<"' to '"<<Parent->Data().toLatin1().constData()<<"'\n";
			}
			return 0;
		}
	}
	else{
		return 0;
	}
}

// Adds a child of type Wave to Parent
GuiWaveTreeModelItem *GuiWaveTreeModel::AddWave(QString Data, GuiWaveTreeModelItem *Parent){
	if( Parent ){
		GuiWaveTreeModelItem *Wave = new GuiWaveTreeModelItem(Data, GuiWaveTreeModelItem::Wave, Parent);
		if( Parent->AppendChild(Wave) ){
			if( _Debug ){
				printDebug << "Wave '"<<Data.toLatin1().constData()<<"' added to '"<<Parent->Data().toLatin1().constData()<<"'\n";
			}
			return Wave;
		}
		else{
			delete Wave;
			if( _Debug ){
				printDebug << "GuiWaveTreeModel::AddWave AppendChild failed\n";
			}
			return 0;
		}
	}
	else{
		return 0;
	}
}

// Calls QAbstractItemModel constructor and creates _RootItem
GuiWaveTreeModel::GuiWaveTreeModel( QObject *Parent ):
	QAbstractItemModel(Parent){
	_RootItem = new GuiWaveTreeModelItem("Master of the Universe", GuiWaveTreeModelItem::Master);
}

// Deletes _RootItem
GuiWaveTreeModel::~GuiWaveTreeModel(){
	delete _RootItem;
}


// Deletes all items of the tree (except _RootItem, which should never be deleted, while object exists)
void GuiWaveTreeModel::Clear(){
	emit beginRemoveRows( QModelIndex(), 0, _RootItem->NChildren() );

	_RootItem->DeleteChildren();

	emit endRemoveRows();
}

// Adds a root tree selected in DataObject to the model
void GuiWaveTreeModel::AddTree(const RootPwaDataObject& DataObject){
	GuiWaveTreeModelItem *File = AddFile( DataObject.DataFileName().c_str() );
	GuiWaveTreeModelItem *Tree = AddTree( DataObject.DataTreeName().c_str(), File);

	const vector<string> *WaveList = DataObject.WavesInTree();
	if( WaveList ){
		for( unsigned int i=0; i < WaveList->size(); ++i ){
			AddWave( WaveList->at(i).c_str(), Tree );
		}
	}

	if( _Debug ){
		printDebug << "Tree added\n";
	}
	// Sends signal so other items can react to it(for example the QTreeView showing the data)
	emit QAbstractItemModel::layoutChanged();
}

// Returns the index of the item in the model specified by the given row, column and parent index.
QModelIndex GuiWaveTreeModel::index(int row, int column, const QModelIndex &parent) const{
	if ( !hasIndex(row, column, parent) ){
		if( _Debug ){
			QModelIndex debug_mindex;
			printDebug << "GuiWaveTreeModel::index("<<row<<';'<<column<<';'<<parent.data().toString().toLatin1().constData()<<"):["<<debug_mindex.row()<<';'<<debug_mindex.column()<<';'<<debug_mindex.parent().data().toString().toLatin1().constData()<<"]\n";
		}
		return QModelIndex();
	}

	GuiWaveTreeModelItem *ParentItem;

	if ( !parent.isValid() ){
		ParentItem = _RootItem;
	}
	else{
    	ParentItem = static_cast<GuiWaveTreeModelItem*>( parent.internalPointer() );
	}
	GuiWaveTreeModelItem *ChildItem = ParentItem->Child(row);
	if (ChildItem){
		if( _Debug ){
			QModelIndex debug_mindex = createIndex(row, column, ChildItem);
			printDebug << "GuiWaveTreeModel::index("<<row<<';'<<column<<';'<<parent.data().toString().toLatin1().constData()<<"):["<<debug_mindex.row()<<';'<<debug_mindex.column()<<';'<<debug_mindex.parent().data().toString().toLatin1().constData()<<"]\n";
		}
		return createIndex(row, column, ChildItem);
	}
	else{
		if( _Debug ){
			QModelIndex debug_mindex;
			printDebug << "GuiWaveTreeModel::index("<<row<<';'<<column<<';'<<parent.data().toString().toLatin1().constData()<<"):["<<debug_mindex.row()<<';'<<debug_mindex.column()<<';'<<debug_mindex.parent().data().toString().toLatin1().constData()<<"]\n";
		}
		return QModelIndex();
	}
}

// Returns the parent of the model item with the given index. If the item has no parent, an invalid QModelIndex is returned.
QModelIndex GuiWaveTreeModel::parent(const QModelIndex &index) const{
	if ( !index.isValid() ){
		if( _Debug ){
			QModelIndex debug_mindex;
			printDebug << "GuiWaveTreeModel::parent(["<<index.row()<<';'<<index.column()<<"]):["<<debug_mindex.row()<<';'<<debug_mindex.column()<<"]\n";
		}
		return QModelIndex();
	}

	GuiWaveTreeModelItem *ChildItem = static_cast<GuiWaveTreeModelItem*>( index.internalPointer() );
	GuiWaveTreeModelItem *ParentItem = ChildItem->Parent();

	if (ParentItem == _RootItem){
		if( _Debug ){
			QModelIndex debug_mindex;
			printDebug << "GuiWaveTreeModel::parent(["<<index.row()<<';'<<index.column()<<"]):["<<debug_mindex.row()<<';'<<debug_mindex.column()<<"]\n";
		}
		return QModelIndex();
	}

	if( _Debug ){
		QModelIndex debug_mindex = createIndex(ParentItem->IndexAtParent(), 0, ParentItem);
		printDebug << "GuiWaveTreeModel::parent(["<<index.row()<<';'<<index.column()<<"]):["<<debug_mindex.row()<<';'<<debug_mindex.column()<<"]\n";
	}
	return createIndex(ParentItem->IndexAtParent(), 0, ParentItem);
}

// Returns the number of rows under the given parent. When the parent is valid it means that rowCount is returning the number of children of parent.
int GuiWaveTreeModel::rowCount(const QModelIndex &parent) const{
	if ( parent.column() > 0 ){
		if( _Debug ){
			printDebug << "GuiWaveTreeModel::rowCount("<<parent.data().toString().toLatin1().constData()<<"):"<<0<<'\n';
		}
		return 0;
	}
	else{
		GuiWaveTreeModelItem *ParentItem;

		if ( !parent.isValid() ){
			ParentItem = _RootItem;
		}
		else{
			ParentItem = static_cast<GuiWaveTreeModelItem*>( parent.internalPointer() );
		}

		if( _Debug ){
			printDebug << "GuiWaveTreeModel::rowCount("<<parent.data().toString().toLatin1().constData()<<"):"<<ParentItem->NChildren()<<'\n';
		}
		return ParentItem->NChildren();
    }
}

// Returns the number of columns for the children of the given parent. (Which is always 1)
int GuiWaveTreeModel::columnCount(const QModelIndex &parent) const{
	// The GuiWaveTreeModelItems only have one string as data and therefore always one row
	if( _Debug ){
		printDebug << "GuiWaveTreeModel::columnCount("<<parent.data().toString().toLatin1().constData()<<"):"<<1<<'\n';
	}

	return 1;
}

// Returns the data stored under the given role for the item referred to by the index.
QVariant GuiWaveTreeModel::data(const QModelIndex &index, int role) const{
	if ( !index.isValid() ){
		if( _Debug ){
			printDebug << "GuiWaveTreeModel::data(["<<index.row()<<';'<<index.column()<<';'<<index.parent().data().toString().toLatin1().constData()<<"];"<<role<<"):"<<QString().toLatin1().constData()<<'\n';
		}
		return QVariant();
	}

	if (role != Qt::DisplayRole){
		if( _Debug ){
			printDebug << "GuiWaveTreeModel::data(["<<index.row()<<';'<<index.column()<<';'<<index.parent().data().toString().toLatin1().constData()<<"];"<<role<<"):"<<QString().toLatin1().constData()<<'\n';
		}
		return QVariant();
	}

	GuiWaveTreeModelItem *Item = static_cast<GuiWaveTreeModelItem*>( index.internalPointer() );

	if( _Debug ){
		printDebug << "GuiWaveTreeModel::data(["<<index.row()<<';'<<index.column()<<';'<<index.parent().data().toString().toLatin1().constData()<<"];"<<role<<"):"<<Item->Data().toLatin1().constData()<<'\n';
	}
	return Item->Data();
}

// Returns the item flags for the given index.
Qt::ItemFlags GuiWaveTreeModel::flags(const QModelIndex &index) const{
	if ( !index.isValid() ){
		if( _Debug ){
			printDebug << "GuiWaveTreeModel::flags(["<<index.row()<<';'<<index.column()<<';'<<index.parent().data().toString().toLatin1().constData()<<"]):"<<0<<'\n';
		}
		return 0;
	}

	if( _Debug ){
		printDebug << "GuiWaveTreeModel::flags("<<index.row()<<';'<<index.column()<<';'<<index.parent().data().toString().toLatin1().constData()<<"]):"<<(Qt::ItemIsEnabled | Qt::ItemIsSelectable)<<'\n';
	}
	return Qt::ItemIsEnabled | Qt::ItemIsSelectable;
}

// Returns the data for the given role and section in the header with the specified orientation. (Just an empty QVariant, since no headers are implemented)
QVariant GuiWaveTreeModel::headerData ( int section, Qt::Orientation orientation, int role ) const{
	return QVariant();
}
