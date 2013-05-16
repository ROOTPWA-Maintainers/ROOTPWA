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
	GuiStringTreeModelItem( Data, Parent ),
	_Type(Type){
		if( Folder ==_Type ){
			_SubFolderList = Data.split('/');
			if( _SubFolderList.size() > 1){
				if( _SubFolderList[0].isEmpty() ){
					_SubFolderList[1].push_front( QChar('/') );
					_SubFolderList.removeAt(0);
				}
			}
		}
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

// Returns the string from _SubFolderList at index
const QString& GuiFileListTreeModelItem::SubFolder(unsigned int index) const{
	return _SubFolderList.at(index);
}

// Returns true and appends child to child list if it's a proper child (type of child fits to type of this item)
bool GuiFileListTreeModelItem::AppendChild(GuiFileListTreeModelItem *Child){
	if( Child && ( ( Folder == Child->Type() ) || ( ( File == Child->Type() ) && ( Folder == _Type) ) ) ){ // Proper child?
		_Children.append(Child);
		if( _Debug ){
			printDebug << "GuiFileListTreeModelItem::AppendChild: Child appended\n";
		}
		return true;
	}
	else{
		if( _Debug ){
			printDebug << "GuiFileListTreeModelItem::AppendChild: Child not proper: (ChildPtr: "<<Child<<"; ChildType: "<<Child->Type()<<"; ParentType: "<<_Type<<")\n";
		}
		return false;
	}
}

// Returns the index of the element which has this string as first string in _SubFolderList or -1 if none has
int GuiFileListTreeModelItem::FindSub( const QString& SearchTerm ) const{
	for( int i = 0; i < _Children.size(); ++i ){
		if( static_cast<GuiFileListTreeModelItem *>(_Children[i])->Type() == Folder ){
			if( static_cast<GuiFileListTreeModelItem *>(_Children[i])->SubFolder(0) == SearchTerm ){
				return i;
			}
		}
	}

	// Has not found a children with this data string
	return -1;
}

// Follows the FollowList beginning at i as long as the folders in it are consistent with the folders in _SubFolderList and increases i in the process. If not the complete list can be followed i is the index of the first non-identical entry in FollowList and the return value is the index in _SubFolderList. If the complete list can be followed i is the index of the element in FollowList after the last of _SubFolderList and the return value is -1.
int GuiFileListTreeModelItem::FollowSub( const QStringList& FollowList, int& i ) const{
	int n = 0;
	while( (n < _SubFolderList.size()) && (i < FollowList.size()) && (FollowList.at(i) == _SubFolderList.at(n)) ){
		++i;
		++n;
	}

	if( _SubFolderList.size() == n ){
		return -1;
	}
	else{
		return n;
	}
}

// Splits this folder type GuiFileListTreeModelItem into two GuiFileListTreeModelItem where the folder at index in the _SubFolderList is the first folder of the new GuiFileListTreeModelItem that is appended as a child to the old GuiFileListTreeModelItem, where the _SubFolderList is cut at index and _Data is adjusted (returns pointer to old GuiFileListTreeModelItem)
GuiFileListTreeModelItem *GuiFileListTreeModelItem::Split( int index ){
	// Create Child
	QStringList NewChildData = _SubFolderList.mid(index);
	GuiFileListTreeModelItem *NewChild = new GuiFileListTreeModelItem(NewChildData.join("/"), GuiFileListTreeModelItem::Folder, this);

	// Delete all folders that have been given to the child
	while ( _SubFolderList.size() > index ){
		_SubFolderList.removeLast();
	}
	_Data = _SubFolderList.join("/");

	// Move all children to the new child
	for( int i=0; i < _Children.size(); ++i ){
		_Children.at(i)->SetParent(NewChild);
		NewChild->AppendChild( static_cast<GuiFileListTreeModelItem *>(_Children.at(i)) );
	}
	_Children.clear();

	// Append the new child
	this->AppendChild(NewChild);

	return this;
}

///< Appends all files for this parent to the end of FileList and returns it
vector<string>& GuiFileListTreeModelItem::GetFiles( vector<string>& FileList ) const{
	if( _Debug ){
		printDebug << "GetFiles: Getting FileList of size: " << FileList.size() << '\n';
	}

	switch( _Type ){
	case Master:
	{
		for( int i=0; i < _Children.size(); ++i ){
			static_cast<GuiFileListTreeModelItem *>(_Children.at(i))->GetFiles( FileList );
		}

		break;
	}
	case Folder:
	{
		if( _Debug ){
			printDebug << "GetFiles: Accessing Folder: " << Data().toLatin1().constData() << '\n';
		}
		unsigned int n = FileList.size(); // Stores the first index where the folder has not been added to the String yet

		for( int i=0; i < _Children.size(); ++i ){
			static_cast<GuiFileListTreeModelItem *>(_Children.at(i))->GetFiles( FileList );

			if( _Debug ){
				printDebug << "GetFiles: Appending Folder: " << Data().toLatin1().constData() << '\n';
			}

			// Pushes the folder name to the front
			for( ; n < FileList.size(); ++n ){
				FileList[n].insert(0,"/");
				FileList[n].insert(0, Data().toLatin1().constData() );
				if( _Debug ){
					printDebug << "GetFiles: Resulting in: " << FileList[n] << '\n';
				}
			}
		}
		break;
	}
	case File:
	{
		if( _Debug ){
			printDebug << "GetFiles: Push back File: " << Data().toLatin1().constData() << '\n';
		}
		FileList.push_back( Data().toLatin1().constData() );
		break;
	}
	}

	if( _Debug ){
		printDebug << "GetFiles: Returning FileList of size: " << FileList.size() << '\n';
	}

	return FileList;
}
