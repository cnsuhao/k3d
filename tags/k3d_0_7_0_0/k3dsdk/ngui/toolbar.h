#ifndef NGUI_TOOLBAR_H
#define NGUI_TOOLBAR_H

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
#include <gtkmm/box.h>

namespace libk3dngui
{

namespace toolbar
{

/////////////////////////////////////////////////////////////////////////////
// control

/// Provides a standard toolbar control that is interactive (can participate in tutorials and macros)
class control :
        public Gtk::VBox,
	public ui_component
{
	typedef Gtk::VBox base;

public:
	control(k3d::icommand_node& Parent, const std::string& Name);

	Gtk::HBox& row(const unsigned int Row);

private:
	void pack_start();
	void pack_end();
//	void pack_start(Gtk::Widget& child, bool expand, bool fill, guint padding);
//	void pack_start(Gtk::Widget& child, Gtk::PackOptions options, guint padding);
//	void pack_end(Gtk::Widget& child, bool expand, bool fill, guint padding);
//	void pack_end(Gtk::Widget& child, Gtk::PackOptions options, guint padding);
};

} // namespace toolbar

} // namespace libk3dngui

#endif // NGUI_TOOLBAR_H


