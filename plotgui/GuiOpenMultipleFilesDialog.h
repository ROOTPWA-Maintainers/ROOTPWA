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
//      Header file for the GuiOpenMultipleFilesDialog class that provides
//		a dialog form where several files can be selected one after
//		another in different directories
//
//
// Author List:
//      Stephan Schmeing          TUM            (original author)
//
//
//-------------------------------------------------------------------------

#ifndef GuiOpenMultipleFilesDialog_H
#define GuiOpenMultipleFilesDialog_H

#include <QObject>
#include <QWidget>
#include <QFileDialog>
#include <QGridLayout>
#include <QListView>
#include <QModelIndex>
#include <QDir>

#include <ui_GuiSelectedFilesWidgetAdd.h>

#include "GuiFileListTreeModel.h"

namespace rpwa{

	class GuiOpenMultipleFilesDialog: public QFileDialog{
	Q_OBJECT
	private:
		class GuiSelectedFilesWidgetAdd: public QWidget, public Ui::GuiSelectedFilesWidgetAdd{
		private:
			// Variables
			GuiFileListTreeModel _SelectedFilesModel;
			QItemSelectionModel _SelectedFilesSelection;

			static bool _Debug; ///< if set to true, debug messages are printed
		public:
			// Variables

			// Constructors + Destructors
			GuiSelectedFilesWidgetAdd( QWidget * parent = 0 ); ///< Sets up the Widget for the selection view, which is added to the file dialog

			// Functions
			QModelIndexList SelectedFiles() const; ///< Returns a list of all selected files
			void AddFileString(const QString& Data); ///< Adds a file string to the model
			void AddFolderString(const QDir& Folder); ///< Adds a complete folder to the model
			void DeleteItem( const QModelIndex& Item ); ///< Deletes an item and all parents that are left without a child (except the _RootItem of course)
			void DeleteSelectedItems(); ///< Deletes all selected items and all parents that are left without a child (except the _RootItem of course)
		};

	private:
		// Variables
		GuiSelectedFilesWidgetAdd _SelectedFilesWidget;
		QListView *_FileDialogListView;

		static bool _Debug; ///< if set to true, debug messages are printed

		// Functions

	private slots:
		// Automatically connected slots

		// Manually connected slots
		void mon__BtAdd_clicked();
		void mon__BtAddAll_clicked();
		void mon__BtRemove_clicked();

	protected:
		// Variables

		// Functions

	public:
		// Constructors + Destructors
		GuiOpenMultipleFilesDialog( QWidget * parent = 0, const QString & caption = QString(), const QString & directory = QString(), const QString & filter = QString() ); ///< Initializes the file dialog

		// Get && Set

		// Functions
	};

} // namespace rpwa

#endif /* GuiOpenMultipleFilesDialog_H */
