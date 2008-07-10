#ifndef NGUI_KEYBOARD_H
#define NGUI_KEYBOARD_H

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
		\brief Declares the libk3dngui::entry control, which provides a standard MVC UI for string values
		\author Tim Shead (tshead@k-3d.com)
*/

#include <k3dsdk/keyboard.h>

namespace libk3dngui
{

/// Converts GDK keyboard modifiers to our native keyboard modifier type
const k3d::key_modifiers convert(const unsigned int Modifiers);
/// Converts our native keyboard modifier type to GDK keyboard modifiers
const unsigned int convert(const k3d::key_modifiers Modifiers);

} // namespace libk3dngui

#endif // !NGUI_KEYBOARD_H
