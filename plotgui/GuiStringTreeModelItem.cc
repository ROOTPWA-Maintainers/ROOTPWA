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
//      Code file for the GuiStringTreeModelItem class that provides an
//		item for the model for QTreeview for a string view
//
//
// Author List:
//      Stephan Schmeing          TUM            (original author)
//
//
//-------------------------------------------------------------------------

#include "reportingUtils.hpp"

#include "GuiStringTreeModelItem.h"

using namespace std;
using namespace rpwa;

bool GuiStringTreeModelItem::_Debug = false;

// Fills Data and Parent into the member variables
GuiStringTreeModelItem::GuiStringTreeModelItem(const QString& Data, GuiStringTreeModelItem *Parent):
	_Data(Data),
	_Parent(Parent){
}

// Deletes all children
GuiStringTreeModelItem::~GuiStringTreeModelItem(){
	qDeleteAll(_Children);
}

// Returns the number of children
int GuiStringTreeModelItem::NChildren() const{
	return _Children.count();
}

// Returns child with index Index
GuiStringTreeModelItem *GuiStringTreeModelItem::Child(int Index){
	return static_cast<GuiStringTreeModelItem *>( _Children.value(Index) );
}

// Returns the parent
GuiStringTreeModelItem *GuiStringTreeModelItem::Parent(){
	return static_cast<GuiStringTreeModelItem *>( _Parent );
}

// Returns the index where this element is listed at its parent
int GuiStringTreeModelItem::IndexAtParent() const{
	if (_Parent){
		return _Parent->_Children.indexOf( const_cast<GuiStringTreeModelItem*>(this) );
	}
	else{
		return 0;
	}
}

// Returns the data string
const QString& GuiStringTreeModelItem::Data() const{
	return _Data;
}

// Deletes child with index Index
void GuiStringTreeModelItem::DeleteChild(int Index){
	delete _Children[Index];

	_Children.removeAt(Index);
}

// Deletes all children of this item
void GuiStringTreeModelItem::DeleteChildren(){
	for( int i=0; i < _Children.size(); ++i){
		delete _Children[i];
	}

	_Children.clear();
}

// Returns the index of the element which has this string as data string or -1 if none has
int GuiStringTreeModelItem::Find( const QString& SearchTerm ) const{
	for( int i = 0; i < _Children.size(); ++i ){
		if( _Children[i]->Data() == SearchTerm ){
			return i;
		}
	}

	// Has not found a children with this data string
	return -1;
}
