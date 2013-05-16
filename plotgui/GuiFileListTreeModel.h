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

#include <vector>
#include <string>

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
		GuiFileListTreeModelItem *AddCompletePath(const QString& Data); ///< Adds all folders in the given path in Data if they do not already exists and returns the GuiFileListTreeModelItem of the last folder in the string

	public:
		// Constructors + Destructors
		GuiFileListTreeModel( QObject *Parent ); ///< Calls GuiStringTreeModel constructor and creates _RootItem
		~GuiFileListTreeModel(); ///< Deletes _RootItem

		// Get && Set
		static bool Debug() { return _Debug; } ///< returns debug flag
		static void SetDebug(const bool Debug = true) { _Debug = Debug; } ///< sets debug flag

		// Functions
		void AddFileString(const QString& Data); ///< Adds a file string to the model
		void AddFolderString(const QDir& Folder); ///< Adds a complete folder to the model
		void DeleteItems( QModelIndexList DeletionList ); ///< Deletes all items on DeletionList and all parents that are left without a child (except the _RootItem of course)
		std::vector<std::string>& GetFiles( std::vector<std::string>& FileList ) const; ///< Appends all files from list to the end of FileList and returns it

		std::ostream& Print(std::ostream& out) const;
	};

} // namespace rpwa

#endif /* GuiFileListTreeModel_H */
