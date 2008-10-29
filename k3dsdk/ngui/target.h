#ifndef K3DSDK_NGUI_TARGET_H
#define K3DSDK_NGUI_TARGET_H

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

namespace libk3dngui
{

class document_state;
namespace viewport { class control; }

/// Centers the current selection in the given viewport by changing the viewport orientation
void aim_selection(document_state& DocumentState, viewport::control& Viewport);
/// Frames the current selection in the given viewport by changing the viewport position
void frame_selection(document_state& DocumentState, viewport::control& Viewport);

} // namespace libk3dngui

#endif // !K3DSDK_NGUI_TARGET_H
