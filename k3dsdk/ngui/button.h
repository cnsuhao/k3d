#ifndef K3DSDK_NGUI_BUTTON_H
#define K3DSDK_NGUI_BUTTON_H

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
#include <k3dsdk/signal_system.h>
#include <gtkmm/button.h>

namespace k3d
{

namespace ngui
{

namespace button
{

/////////////////////////////////////////////////////////////////////////////
// control

/// Provides a standard button control that is interactive (can participate in tutorials and macros)
class control :
        public Gtk::Button,
	public ui_component
{
	typedef Gtk::Button base;

public:
	control();
	control(const Glib::ustring& label, bool mnemonic = false);
	control(const Gtk::StockID& stock_id);
	control(const Glib::ustring& label, const Gtk::StockID& stock_id);
	control(Gtk::Widget& widget);

	void on_clicked();
};

} // namespace button

} // namespace ngui

} // namespace k3d

#endif // !K3DSDK_NGUI_BUTTON_H

