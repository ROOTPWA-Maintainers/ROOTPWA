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
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//      Header file for the GuiWaveTreeModelItem class that provides an
//		item for the model for QTreeview for the wave view in GuiPwaMain
//
//
// Author List:
//      Stephan Schmeing          TUM            (original author)
//
//
//-------------------------------------------------------------------------

#ifndef GuiWaveTreeModelItem_H
#define GuiWaveTreeModelItem_H

#include <QList>
#include <QVariant>

namespace rpwa{

	class GuiWaveTreeModelItem{
	public:
		// Enums
		enum E_Types{
			Master,
			File,
			Tree,
			Wave
		};

	private:
		// Variables
	     QList<GuiWaveTreeModelItem*> _Children;
	     QVariant _Data;
	     GuiWaveTreeModelItem *_Parent;
	     E_Types _Type;

		static bool _Debug; ///< if set to true, debug messages are printed

		// Functions

	public:
		// Constructors + Destructors
		GuiWaveTreeModelItem(const QVariant& Data, E_Types Type, GuiWaveTreeModelItem *Parent = 0); ///< Fills Data, Type and Parent into the member variables
		~GuiWaveTreeModelItem(); ///< Deletes all children

		// Get && Set
		int NChildren() const; ///< Returns the number of children
		GuiWaveTreeModelItem *Child(int Row); ///< Returns child with index Index
		GuiWaveTreeModelItem *Parent(); ///< Returns the parent
		int IndexAtParent() const; ///< Returns the index where this element is listed at its parent
		E_Types Type() const; ///< Returns the type of this item
		const QVariant& Data() const; ///< Returns the data string

		// Functions
		bool AppendChild(GuiWaveTreeModelItem *Child); ///< Returns true and appends child to child list if it's a proper child (type of child fits to type of this item)
		void DeleteChildren(); ///< Deletes all Children of this item
	};

} // namespace rpwa

#endif /* GuiWaveTreeModelItem_H */
