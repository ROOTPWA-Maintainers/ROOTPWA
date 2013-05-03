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
//      Header file for the GuiStringTreeModel class that provides a model
//		for QTreeview for a string view
//
//
// Author List:
//      Stephan Schmeing          TUM            (original author)
//
//
//-------------------------------------------------------------------------

#ifndef GuiStringTreeModel_H
#define GuiStringTreeModel_H

#include <QAbstractItemModel>
#include <QModelIndex>
#include <QVariant>

#include "GuiStringTreeModelItem.h"

namespace rpwa{

	class GuiStringTreeModel: public QAbstractItemModel{
	Q_OBJECT
	private:
		// Variables

		static bool _Debug; ///< if set to true, debug messages are printed

		// Functions
		virtual GuiStringTreeModelItem *RootItem() const = 0;
		virtual GuiStringTreeModelItem *RootItem() = 0;

	public:
		// Constructors + Destructors
		GuiStringTreeModel( QObject *Parent ); ///< Calls QAbstractItemModel constructor and creates _RootItem

		// Get && Set

		// Functions
		void Clear(); ///< Deletes all items of the tree (except _RootItem, which should never be deleted, while object exists)

		// Overloaded QAbstractItemModel functions
		QModelIndex index(int row, int column, const QModelIndex &parent = QModelIndex()) const; ///< Returns the index of the item in the model specified by the given row, column and parent index.
		QModelIndex parent(const QModelIndex &index) const; ///< Returns the parent of the model item with the given index. If the item has no parent, an invalid QModelIndex is returned.
		int rowCount(const QModelIndex &parent = QModelIndex()) const; ///< Returns the number of rows under the given parent. When the parent is valid it means that rowCount is returning the number of children of parent.
		int columnCount(const QModelIndex &parent = QModelIndex()) const; ///< Returns the number of columns for the children of the given parent. (Which is always 1)
		QVariant data(const QModelIndex &index, int role) const; ///< Returns the data stored under the given role for the item referred to by the index.
		Qt::ItemFlags flags(const QModelIndex &index) const; ///< Returns the item flags for the given index.
		QVariant headerData ( int section, Qt::Orientation orientation, int role = Qt::DisplayRole ) const; ///< Returns the data for the given role and section in the header with the specified orientation. (Just an empty QVariant, since no headers are implemented)
	};

} // namespace rpwa

#endif /* GuiStringTreeModel_H */
