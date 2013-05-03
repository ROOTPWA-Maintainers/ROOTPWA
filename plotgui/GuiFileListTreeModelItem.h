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

#include <QList>
#include <QVariant>
#include <QString>

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

		static bool _Debug; ///< if set to true, debug messages are printed

		// Functions

	public:
		// Constructors + Destructors
		GuiFileListTreeModelItem(const QString& Data, E_Types Type, GuiFileListTreeModelItem *Parent = 0); ///< Fills Data, Type and Parent into the member variables

		// Get && Set
		GuiFileListTreeModelItem *Child(int Row); ///< Returns child with index Index
		GuiFileListTreeModelItem *Parent(); ///< Returns the parent
		E_Types Type() const; ///< Returns the type of this item

		// Functions
		bool AppendChild(GuiFileListTreeModelItem *Child); ///< Returns true and appends child to child list if it's a proper child (type of child fits to type of this item)
	};

} // namespace rpwa

#endif /* GuiFileListTreeModelItem_H */
