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
//      Code file for the GuiWaveTreeModelItem class that provides an
//		item for the model for QTreeview for the wave view in GuiPwaMain
//
//
// Author List:
//      Stephan Schmeing          TUM            (original author)
//
//
//-------------------------------------------------------------------------

#include "reportingUtils.hpp"

#include "GuiWaveTreeModelItem.h"

using namespace std;
using namespace rpwa;

bool GuiWaveTreeModelItem::_Debug = false;

// Fills Data, Type and Parent into the member variables
GuiWaveTreeModelItem::GuiWaveTreeModelItem(const QString& Data, E_Types Type, GuiWaveTreeModelItem *Parent):
	_Data(Data),
	_Parent(Parent),
	_Type(Type){
}

// Deletes all children
GuiWaveTreeModelItem::~GuiWaveTreeModelItem(){
	qDeleteAll(_Children);
}

// Returns the number of children
int GuiWaveTreeModelItem::NChildren() const{
	return _Children.count();
}

// Returns child with index Index
GuiWaveTreeModelItem *GuiWaveTreeModelItem::Child(int Index){
	return _Children.value(Index);
}

// Returns the parent
GuiWaveTreeModelItem *GuiWaveTreeModelItem::Parent(){
	return _Parent;
}

// Returns the index where this element is listed at its parent
int GuiWaveTreeModelItem::IndexAtParent() const{
	if (_Parent){
		return _Parent->_Children.indexOf( const_cast<GuiWaveTreeModelItem*>(this) );
	}
	else{
		return 0;
	}
}

// Returns the type of this item
GuiWaveTreeModelItem::E_Types GuiWaveTreeModelItem::Type() const{
	return _Type;
}

// Returns the data string
const QString& GuiWaveTreeModelItem::Data() const{
	return _Data;
}

// Returns true and appends child to child list if it's a proper child (type of child fits to type of this item)
bool GuiWaveTreeModelItem::AppendChild(GuiWaveTreeModelItem *Child){
	if( Child && ( ( Child->Type() - 1 ) == _Type ) ){ // Proper child?
		_Children.append(Child);
		return true;
	}
	else{
		if( _Debug ){
			printDebug << "GuiWaveTreeModelItem::AppendChild: Child not proper: (ChildPtr: "<<Child<<"; ChildType: "<<Child->Type()<<"; ParentType: "<<_Type<<")\n";
		}
		return false;
	}
}

// Deletes all Children of this item
void GuiWaveTreeModelItem::DeleteChildren(){
	for( int i=0; i < _Children.size(); ++i){
		delete _Children[i];
	}

	_Children.clear();
}

