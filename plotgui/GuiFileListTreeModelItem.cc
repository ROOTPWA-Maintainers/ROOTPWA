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
//      Code file for the GuiFileListTreeModelItem class that provides an
//		item for the model for QTreeview for a file list view
//
//
// Author List:
//      Stephan Schmeing          TUM            (original author)
//
//
//-------------------------------------------------------------------------

#include "reportingUtils.hpp"

#include "GuiFileListTreeModelItem.h"

using namespace std;
using namespace rpwa;

bool GuiFileListTreeModelItem::_Debug = false;

// Fills Data, Type and Parent into the member variables
GuiFileListTreeModelItem::GuiFileListTreeModelItem(const QString& Data, E_Types Type, GuiFileListTreeModelItem *Parent):
	GuiStringTreeModelItem( Data, Parent),
	_Type(Type){
}

// Returns child with index Index
GuiFileListTreeModelItem *GuiFileListTreeModelItem::Child(int Index){
	return static_cast<GuiFileListTreeModelItem *>( _Children.value(Index) );
}

// Returns the parent
GuiFileListTreeModelItem *GuiFileListTreeModelItem::Parent(){
	return static_cast<GuiFileListTreeModelItem *>( _Parent );
}

// Returns the type of this item
GuiFileListTreeModelItem::E_Types GuiFileListTreeModelItem::Type() const{
	return _Type;
}

// Returns true and appends child to child list if it's a proper child (type of child fits to type of this item)
bool GuiFileListTreeModelItem::AppendChild(GuiFileListTreeModelItem *Child){
	if( Child && ( ( Folder == Child->Type() ) || ( ( File == Child->Type() ) && ( Folder == _Type) ) ) ){ // Proper child?
		_Children.append(Child);
		return true;
	}
	else{
		if( _Debug ){
			printDebug << "GuiFileListTreeModelItem::AppendChild: Child not proper: (ChildPtr: "<<Child<<"; ChildType: "<<Child->Type()<<"; ParentType: "<<_Type<<")\n";
		}
		return false;
	}
}

