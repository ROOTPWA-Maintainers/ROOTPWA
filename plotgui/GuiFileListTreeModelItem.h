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
//      Header file for the GuiFileListTreeModelItem class that provides an
//		item for the model for QTreeview for a file list view
//
//
// Author List:
//      Stephan Schmeing          TUM            (original author)
//
//
//-------------------------------------------------------------------------

#ifndef GuiFileListTreeModelItem_H
#define GuiFileListTreeModelItem_H

#include <vector>
#include <string>

#include <QList>
#include <QVariant>
#include <QString>
#include <QStringList>

#include "GuiStringTreeModelItem.h"

namespace rpwa{

	class GuiFileListTreeModelItem : public GuiStringTreeModelItem{
	public:
		// Enums
		enum E_Types{
			Master,
			Folder,
			File,
		};

	private:
		// Variables
		E_Types _Type;
		QStringList _SubFolderList; ///< List of folders that are included in this folder string (in case it is a folder, otherwise it's empty)

		static bool _Debug; ///< if set to true, debug messages are printed

		// Functions

	public:
		// Constructors + Destructors
		GuiFileListTreeModelItem(const QString& Data, E_Types Type, GuiFileListTreeModelItem *Parent = 0); ///< Fills Data, Type and Parent into the member variables

		// Get && Set
		GuiFileListTreeModelItem *Child(int Row); ///< Returns child with index Index
		GuiFileListTreeModelItem *Parent(); ///< Returns the parent
		E_Types Type() const; ///< Returns the type of this item
		const QString& SubFolder(unsigned int index) const; ///< Returns the string from _SubFolderList at index

		static bool Debug() { return _Debug; } ///< returns debug flag
		static void SetDebug(const bool Debug = true) { _Debug = Debug; } ///< sets debug flag

		// Functions
		bool AppendChild(GuiFileListTreeModelItem *Child); ///< Returns true and appends child to child list if it's a proper child (type of child fits to type of this item)
		int FindSub( const QString& SearchTerm ) const; ///< Returns the index of the element which has this string as first string in _SubFolderList or -1 if none has
		int FollowSub( const QStringList& FollowList, int& i ) const; ///< Follows the FollowList beginning at i as long as the folders in it are consistent with the folders in _SubFolderList and increases i in the process. If not the complete list can be followed i is the index of the last identical entry in FollowList and the return value is the index in _SubFolderList. If the complete list can be followed i is the index of the element in FollowList that is the last of _SubFolderList and the return value is -1.
		GuiFileListTreeModelItem *Split( int index ); ///< Splits this folder type GuiFileListTreeModelItem into two GuiFileListTreeModelItem where the folder at index in the _SubFolderList is the first folder of the new GuiFileListTreeModelItem that is appended as a child to the old GuiFileListTreeModelItem, where the _SubFolderList is cut at index and _Data is adjusted (returns pointer to old GuiFileListTreeModelItem)
		std::vector<std::string>& GetFiles( std::vector<std::string>& FileList ) const; ///<  Appends all files for this parent to the end of FileList and returns it
	};

} // namespace rpwa

#endif /* GuiFileListTreeModelItem_H */
