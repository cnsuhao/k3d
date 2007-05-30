#ifndef CONICS_DETAIL_H
#define CONICS_DETAIL_H

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
		\author Tim Shead (tshead@k-3d.com)
*/

#include <k3dsdk/bounded.h>
#include <k3dsdk/classes.h>
#include <k3dsdk/geometry.h>
#include <k3dsdk/drawable_gl.h>
#include <k3dsdk/imaterial.h>
#include <k3dsdk/material_client.h>
#include <k3dsdk/measurement.h>
#include <k3dsdk/module.h>
#include <k3dsdk/node.h>
#include <k3dsdk/node_change_signal.h>
#include <k3dsdk/persistent.h>
#include <k3dsdk/renderable_ri.h>
#include <k3dsdk/selection.h>
#include <k3dsdk/snappable.h>
#include <k3dsdk/snap_source.h>
#include <k3dsdk/snap_target.h>
#include <k3dsdk/transformable.h>

namespace libk3dconics
{

/////////////////////////////////////////////////////////////////////////////
// conic

class conic :
	public k3d::snappable<k3d::gl::drawable<k3d::ri::renderable<k3d::material_client<k3d::bounded<k3d::transformable<k3d::persistent<k3d::node_change_signal<k3d::node> > > > > > > >
{
	typedef k3d::snappable<k3d::gl::drawable<k3d::ri::renderable<k3d::material_client<k3d::bounded<k3d::transformable<k3d::persistent<k3d::node_change_signal<k3d::node> > > > > > > > base;

public:
	conic(k3d::iplugin_factory& Factory, k3d::idocument& Document) :
		base(Factory, Document)
	{
		add_snap_source(new k3d::snap_source(_("Center"), sigc::mem_fun(*this, &conic::center_source_position), sigc::mem_fun(*this, &conic::center_source_orientation)));
		add_snap_target(new k3d::snap_target(_("Center"), sigc::mem_fun(*this, &conic::center_target_position), sigc::mem_fun(*this, &conic::center_target_orientation)));
	}

private:
	const k3d::point3 center_source_position()
	{
		return k3d::point3();
	}

	bool center_source_orientation(k3d::vector3& Look, k3d::vector3& Up)
	{
		Look = k3d::vector3(0, 0, 1);
		Up = k3d::vector3(0, 1, 0);
		return true;
	}

	bool center_target_position(const k3d::point3& Position, k3d::point3& TargetPosition)
	{
		TargetPosition = k3d::point3();
		return true;
	}

	bool center_target_orientation(const k3d::point3& Position, k3d::vector3& Look, k3d::vector3& Up)
	{
		return false;
	}
};

} // namespace libk3dconics

#endif // CONICS_DETAIL_H
