#ifndef NGUI_SAVABLE_DOCUMENT_WINDOW_H
#define NGUI_SAVABLE_DOCUMENT_WINDOW_H

// K-3D
// Copyright (c) 1995-2004, Timothy M. Shead
//
// Contact: tshead@k-3d.com
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

/** \file
		\author Tim Shead (tshead@k-3d.com)
*/

#include "ui_component.h"
#include "unsaved_document.h"

#include <gtkmm/window.h>

namespace k3d { class idocument; }

namespace libk3dngui
{

class document_state;
	
///////////////////////////////////////////////////////////////////////////////
// savable_document_window

/// Provides a standard document-level window that is interactive (can participate in tutorials and macros)
class savable_document_window :
        public Gtk::Window,
	public ui_component,
	public unsaved_document
{
	typedef Gtk::Window base;

public:
	savable_document_window(document_state& Document, const std::string& Name);
	virtual ~savable_document_window();

	k3d::idocument& document();

	bool on_key_press_event(GdkEventKey* event);
	bool on_delete_event(GdkEventAny* event);

	/// Gives the user a chance to save changes, returns true iff all changes have been saved.
	const bool save_changes();
	/// Closes the window safely, giving the user a chance to save changes or cancel.
	void safe_close();
	/// Closes the window.
	void close();

private:
	unsaved_document* on_safe_close();

	/// Override in derived classes to handle cleanup when the window is closed
	virtual void on_close();

	document_state& m_document;
};

} // namespace libk3dngui

#endif // NGUI_SAVABLE_DOCUMENT_WINDOW_H

