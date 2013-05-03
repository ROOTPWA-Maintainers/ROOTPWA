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
//      Header file for the GuiFileListTreeModel class that provides a model
//		for QTreeview for a file list view
//
//
// Author List:
//      Stephan Schmeing          TUM            (original author)
//
//
//-------------------------------------------------------------------------

#ifndef GuiFileListTreeModel_H
#define GuiFileListTreeModel_H

#include <QModelIndex>
#include <QVariant>
#include <QString>
#include <QDir>

#include "GuiStringTreeModel.h"
#include "GuiFileListTreeModelItem.h"

namespace rpwa{

	class GuiFileListTreeModel: public GuiStringTreeModel{
	private:
		// Variables
		GuiFileListTreeModelItem *_RootItem;

		static bool _Debug; ///< if set to true, debug messages are printed

		// Functions
		GuiStringTreeModelItem *RootItem() const; ///< Returns the RootItem for the base class
		GuiStringTreeModelItem *RootItem(); ///< Returns the RootItem for the base class
		GuiFileListTreeModelItem *AddFolder(const QString& Data, GuiFileListTreeModelItem *Parent); ///< Adds a child of type Folder to Parent
		GuiFileListTreeModelItem *AddFile(const QString& Data, GuiFileListTreeModelItem *Parent); ///< Adds a child of type File to Parent

	public:
		// Constructors + Destructors
		GuiFileListTreeModel( QObject *Parent ); ///< Calls GuiStringTreeModel constructor and creates _RootItem
		~GuiFileListTreeModel(); ///< Deletes _RootItem

		// Get && Set

		// Functions
		void AddFileString(const QString& Data); ///< Adds a file string to the model
		void AddFolderString(const QDir& Folder); ///< Adds a complete folder to the model
		void DeleteItem( const QModelIndex& Item ); ///< Deletes an Item and all parents that are left without a child (except the _RootItem of course)
	};

} // namespace rpwa

#endif /* GuiFileListTreeModel_H */
