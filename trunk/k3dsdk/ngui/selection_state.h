#ifndef K3DSDK_NGUI_SELECTION_STATE_H
#define K3DSDK_NGUI_SELECTION_STATE_H

// K-3D
// Copyright (c) 1995-2009, Timothy M. Shead
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

#include <k3dsdk/ngui/selection.h>
#include <k3dsdk/selection.h>

#include <vector>

namespace k3d
{

class idocument;
class inode;

namespace ngui
{

namespace selection
{

/// Provides convenience methods for manipulating the selection state of a document.
class state
{
public:
	state(idocument& Document);

	/// Returns the current selection mode.
	const mode current_mode();
	/// Sets the current selection mode.
	void set_current_mode(const mode Mode);
	/// Connect a slot to a signal that will be emitted whenever the selection mode changes.
	sigc::connection connect_current_mode_changed_signal(const sigc::slot<void, ihint*>& Slot);

	/// Returns the current set of selected nodes.
	const std::vector<inode*> selected_nodes();

	/// Selects one node.
	void select(inode& Node);
	/// Selects a collection of components.
	void select(const k3d::selection::record& Selection);
	/// Selects a collection of components.
	void select(const k3d::selection::records& Selection);
	/// Selects all nodes or components (depending on current selection mode).
	void select_all();
	/// Selects all nodes.
	void select_all_nodes();

	/// Deselects one node.
	void deselect(inode& Node);
	/// Deselects a collection of components.
	void deselect(const k3d::selection::record& Selection);
	/// Deselects a collection of components.
	void deselect(const k3d::selection::records& Selection);
	/// Deselects all nodes or components (depending on current selection mode).
	void deselect_all();
	/// Deselects all nodes.
	void deselect_all_nodes();

private:
	class implementation;
	implementation& internal;
};

} // namespace selection

} // namespace ngui

} // namespace k3d

#endif // !K3DSDK_NGUI_SELECTION_STATE_H

