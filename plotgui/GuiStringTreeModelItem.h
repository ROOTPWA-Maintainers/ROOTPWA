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
//      Header file for the GuiStringTreeModelItem class that provides an
//		item for the model for QTreeview for a string view
//
//
// Author List:
//      Stephan Schmeing          TUM            (original author)
//
//
//-------------------------------------------------------------------------

#ifndef GuiStringTreeModelItem_H
#define GuiStringTreeModelItem_H

#include <set>

#include <QList>
#include <QVariant>
#include <QString>

namespace rpwa{

	class GuiStringTreeModelItem{
	private:
		// Variables
		static bool _Debug; ///< if set to true, debug messages are printed

		// Functions

	protected:
		// Variables
		QList<GuiStringTreeModelItem*> _Children;
		QString _Data;
		GuiStringTreeModelItem *_Parent;

	public:
		// Constructors + Destructors
		GuiStringTreeModelItem(const QString& Data, GuiStringTreeModelItem *Parent = 0); ///< Fills Data and Parent into the member variables
		~GuiStringTreeModelItem(); ///< Deletes all children

		// Get && Set
		int NChildren() const; ///< Returns the number of children
		GuiStringTreeModelItem *Child(int Index); ///< Returns child with index Index
		GuiStringTreeModelItem *Parent(); ///< Returns the parent
		void SetParent(GuiStringTreeModelItem *Parent); ///< Sets the parent
		int IndexAtParent() const; ///< Returns the index where this element is listed at its parent
		const QString& Data() const; ///< Returns the data string

		static bool Debug() { return _Debug; } ///< returns debug flag
		static void SetDebug(const bool Debug = true) { _Debug = Debug; } ///< sets debug flag

		// Functions
		void DeleteChild(int Index); ///< Deletes child with index Index
		void DeleteChildren(); ///< Deletes all children of this item
		int Find( const QString& SearchTerm ) const; ///< Returns the index of the element which has this string as data string or -1 if none has
		void AddChildrenToSet( std::set<GuiStringTreeModelItem const*>& ChildrenSet ) const; ///< Adds pointers to all children (and subchildren and subsubchildren ...) to ChildrenSet
	};

} // namespace rpwa

#endif /* GuiStringTreeModelItem_H */
