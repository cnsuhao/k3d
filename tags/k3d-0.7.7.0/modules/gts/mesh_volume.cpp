// K-3D
// Copyright (c) 1995-2007, Timothy M. Shead
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
	\author Romain Behar (romainbehar@yahoo.com)
*/

#include <k3d-i18n-config.h>

#include <k3dsdk/document_plugin_factory.h>
#include <k3dsdk/measurement.h>
#include <k3dsdk/mesh_modifier.h>
#include <k3dsdk/node.h>

#include "helpers.h"

namespace module
{

namespace gts
{

/////////////////////////////////////////////////////////////////////////////
// mesh_volume

class mesh_volume :
	public k3d::mesh_modifier<k3d::node >
{
	typedef k3d::mesh_modifier<k3d::node > base;

public:
	mesh_volume(k3d::iplugin_factory& Factory, k3d::idocument& Document) :
		base(Factory, Document),
		m_volume(init_owner(*this) + init_name("volume") + init_label(_("Volume")) + init_description(_("Mesh volume")) + init_slot(sigc::mem_fun(*this, &mesh_volume::get_volume)))
	{
		m_input_mesh.changed_signal().connect(m_volume.make_reset_slot());
	}

	void on_create_mesh(const k3d::mesh& Input, k3d::mesh& Output)
	{
		Output = Input;
	}

	void on_update_mesh(const k3d::mesh& Input, k3d::mesh& Output)
	{
	}

	double get_volume()
	{
		if(k3d::mesh* const mesh = m_input_mesh.pipeline_value())
		{
			const gts_ptr<GtsSurface> gts_surface = convert(*mesh);
			return gts_surface_volume(gts_surface);
		}

		return 0;
	}

	static k3d::iplugin_factory& get_factory()
	{
		static k3d::document_plugin_factory<mesh_volume,
			k3d::interface_list<k3d::imesh_source,
			k3d::interface_list<k3d::imesh_sink > > > factory(
				k3d::uuid(0x905f24ef, 0x344f1c6f, 0x3242b79d, 0xeb8b0a57),
				"GTSMeshVolume",
				_("Calculates input mesh volume"),
				"Metrics",
				k3d::iplugin_factory::STABLE);

		return factory;
	}

private:
	k3d_data(double, immutable_name, change_signal, no_undo, computed_storage, no_constraint, read_only_property, no_serialization) m_volume;
};

/////////////////////////////////////////////////////////////////////////////
// mesh_volume_factory

k3d::iplugin_factory& mesh_volume_factory()
{
	return mesh_volume::get_factory();
}

} // namespace gts

} // namespace module
