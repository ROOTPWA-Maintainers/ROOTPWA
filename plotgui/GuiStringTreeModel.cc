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
//      Code file for the GuiStringTreeModel class that provides a model
//		for QTreeview for a string view in GuiPwaMain
//
//
// Author List:
//      Stephan Schmeing          TUM            (original author)
//
//
//-------------------------------------------------------------------------

#include <string.h>

#include "reportingUtils.hpp"

#include "GuiStringTreeModel.h"

using namespace std;
using namespace rpwa;

bool GuiStringTreeModel::_Debug = false;

///< Appends the ModelIndex of all children (and subchildren and subsubchildren ...) of Parent to List
void GuiStringTreeModel::AppendChildren( QModelIndexList& List, const QModelIndex& Parent) const{
	for(int i=0; i < rowCount(Parent); ++i){
		List.append( index( i, 0, Parent ) );
		AppendChildren( List, List.last() );
	}
}


///< Returns a list of all entries in the tree
QModelIndexList GuiStringTreeModel::ListOfTree() const{
	QModelIndexList ItemList;

	AppendChildren( ItemList, QModelIndex() ); // RootItem has no ModelIndex

	return ItemList;
}

// Calls QAbstractItemModel constructor
GuiStringTreeModel::GuiStringTreeModel( QObject *Parent ):
		QAbstractItemModel(Parent){
}

// Deletes all items of the tree (except _RootItem, which should never be deleted, while object exists)
void GuiStringTreeModel::Clear(){
	emit beginRemoveRows( QModelIndex(), 0, RootItem()->NChildren() );

	RootItem()->DeleteChildren();

	emit endRemoveRows();
}

// Returns the index of the item in the model specified by the given row, column and parent index.
QModelIndex GuiStringTreeModel::index(int row, int column, const QModelIndex &parent) const{
	if ( !hasIndex(row, column, parent) ){
		if( _Debug ){
			QModelIndex debug_mindex;
			printDebug << "GuiStringTreeModel::index("<<row<<';'<<column<<';'<<parent.data().toString().toLatin1().constData()<<"):["<<debug_mindex.row()<<';'<<debug_mindex.column()<<';'<<debug_mindex.parent().data().toString().toLatin1().constData()<<"]\n";
		}
		return QModelIndex();
	}

	GuiStringTreeModelItem *ParentItem;

	if ( !parent.isValid() ){
		ParentItem = RootItem();
	}
	else{
    	ParentItem = static_cast<GuiStringTreeModelItem*>( parent.internalPointer() );
	}
	GuiStringTreeModelItem *ChildItem = ParentItem->Child(row);
	if (ChildItem){
		if( _Debug ){
			QModelIndex debug_mindex = createIndex(row, column, ChildItem);
			printDebug << "GuiStringTreeModel::index("<<row<<';'<<column<<';'<<parent.data().toString().toLatin1().constData()<<"):["<<debug_mindex.row()<<';'<<debug_mindex.column()<<';'<<debug_mindex.parent().data().toString().toLatin1().constData()<<"]\n";
		}
		return createIndex(row, column, ChildItem);
	}
	else{
		if( _Debug ){
			QModelIndex debug_mindex;
			printDebug << "GuiStringTreeModel::index("<<row<<';'<<column<<';'<<parent.data().toString().toLatin1().constData()<<"):["<<debug_mindex.row()<<';'<<debug_mindex.column()<<';'<<debug_mindex.parent().data().toString().toLatin1().constData()<<"]\n";
		}
		return QModelIndex();
	}
}

// Returns the parent of the model item with the given index. If the item has no parent, an invalid QModelIndex is returned.
QModelIndex GuiStringTreeModel::parent(const QModelIndex &index) const{
	if ( !index.isValid() ){
		if( _Debug ){
			QModelIndex debug_mindex;
			printDebug << "GuiStringTreeModel::parent(["<<index.row()<<';'<<index.column()<<"]):["<<debug_mindex.row()<<';'<<debug_mindex.column()<<"]\n";
		}
		return QModelIndex();
	}

	GuiStringTreeModelItem *ChildItem = static_cast<GuiStringTreeModelItem*>( index.internalPointer() );
	GuiStringTreeModelItem *ParentItem = ChildItem->Parent();

	if (ParentItem == RootItem()){
		if( _Debug ){
			QModelIndex debug_mindex;
			printDebug << "GuiStringTreeModel::parent(["<<index.row()<<';'<<index.column()<<"]):["<<debug_mindex.row()<<';'<<debug_mindex.column()<<"]\n";
		}
		return QModelIndex();
	}

	if( _Debug ){
		QModelIndex debug_mindex = createIndex(ParentItem->IndexAtParent(), 0, ParentItem);
		printDebug << "GuiStringTreeModel::parent(["<<index.row()<<';'<<index.column()<<"]):["<<debug_mindex.row()<<';'<<debug_mindex.column()<<"]\n";
	}
	return createIndex(ParentItem->IndexAtParent(), 0, ParentItem);
}

// Returns the number of rows under the given parent. When the parent is valid it means that rowCount is returning the number of children of parent.
int GuiStringTreeModel::rowCount(const QModelIndex &parent) const{
	if ( parent.column() > 0 ){
		if( _Debug ){
			printDebug << "GuiStringTreeModel::rowCount("<<parent.data().toString().toLatin1().constData()<<"):"<<0<<'\n';
		}
		return 0;
	}
	else{
		GuiStringTreeModelItem *ParentItem;

		if ( !parent.isValid() ){
			ParentItem = RootItem();
		}
		else{
			ParentItem = static_cast<GuiStringTreeModelItem*>( parent.internalPointer() );
		}

		if( _Debug ){
			printDebug << "GuiStringTreeModel::rowCount("<<parent.data().toString().toLatin1().constData()<<"):"<<ParentItem->NChildren()<<'\n';
		}
		return ParentItem->NChildren();
    }
}

// Returns the number of columns for the children of the given parent. (Which is always 1)
int GuiStringTreeModel::columnCount(const QModelIndex &parent) const{
	// The GuiWaveTreeModelItems only have one string as data and therefore always one row
	if( _Debug ){
		printDebug << "GuiStringTreeModel::columnCount("<<parent.data().toString().toLatin1().constData()<<"):"<<1<<'\n';
	}

	return 1;
}

// Returns the data stored under the given role for the item referred to by the index.
QVariant GuiStringTreeModel::data(const QModelIndex &index, int role) const{
	if ( !index.isValid() ){
		if( _Debug ){
			printDebug << "GuiStringTreeModel::data(["<<index.row()<<';'<<index.column()<<';'<<index.parent().data().toString().toLatin1().constData()<<"];"<<role<<"):"<<QString().toLatin1().constData()<<'\n';
		}
		return QVariant();
	}

	if ( (role != Qt::DisplayRole) && (role != Qt::ToolTipRole) ){
		if( _Debug ){
			printDebug << "GuiStringTreeModel::data(["<<index.row()<<';'<<index.column()<<';'<<index.parent().data().toString().toLatin1().constData()<<"];"<<role<<"):"<<QString().toLatin1().constData()<<'\n';
		}
		return QVariant();
	}

	GuiStringTreeModelItem *Item = static_cast<GuiStringTreeModelItem*>( index.internalPointer() );

	if( _Debug ){
		printDebug << "GuiStringTreeModel::data(["<<index.row()<<';'<<index.column()<<';'<<index.parent().data().toString().toLatin1().constData()<<"];"<<role<<"):"<<Item->Data().toLatin1().constData()<<'\n';
	}
	return Item->Data();
}

// Returns the item flags for the given index.
Qt::ItemFlags GuiStringTreeModel::flags(const QModelIndex &index) const{
	if ( !index.isValid() ){
		if( _Debug ){
			printDebug << "GuiStringTreeModel::flags(["<<index.row()<<';'<<index.column()<<';'<<index.parent().data().toString().toLatin1().constData()<<"]):"<<0<<'\n';
		}
		return 0;
	}

	if( _Debug ){
		printDebug << "GuiStringTreeModel::flags("<<index.row()<<';'<<index.column()<<';'<<index.parent().data().toString().toLatin1().constData()<<"]):"<<(Qt::ItemIsEnabled | Qt::ItemIsSelectable)<<'\n';
	}
	return Qt::ItemIsEnabled | Qt::ItemIsSelectable;
}

// Returns the data for the given role and section in the header with the specified orientation. (Just an empty QVariant, since no headers are implemented)
QVariant GuiStringTreeModel::headerData ( int section, Qt::Orientation orientation, int role ) const{
	return QVariant();
}
