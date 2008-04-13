// K-3D
// Copyright (c) 1995-2006, Timothy M. Shead
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
	\author Timothy M. Shead (tshead@k-3d.com)
*/

#include <k3dsdk/module.h>

namespace libk3dvirtualglpainters
{

extern k3d::iplugin_factory& edge_painter_factory();
extern k3d::iplugin_factory& face_painter_factory();
extern k3d::iplugin_factory& point_painter_factory();
extern k3d::iplugin_factory& sds_edge_painter_factory();
extern k3d::iplugin_factory& sds_face_painter_factory();
extern k3d::iplugin_factory& sds_point_painter_factory();

} // namespace libk3dvirtualglpainters

K3D_MODULE_START(Registry)
	Registry.register_factory(libk3dvirtualglpainters::edge_painter_factory());
	Registry.register_factory(libk3dvirtualglpainters::face_painter_factory());
	Registry.register_factory(libk3dvirtualglpainters::point_painter_factory());
	Registry.register_factory(libk3dvirtualglpainters::sds_edge_painter_factory());
	Registry.register_factory(libk3dvirtualglpainters::sds_face_painter_factory());
	Registry.register_factory(libk3dvirtualglpainters::sds_point_painter_factory());
K3D_MODULE_END
