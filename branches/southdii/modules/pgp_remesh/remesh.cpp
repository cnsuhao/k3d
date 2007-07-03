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
		\brief 
		\author Ian South-Dickinson (ian.southd@gmail.com)
*/
#include <k3dsdk/document_plugin_factory.h>
#include <k3dsdk/log.h>
#include <k3dsdk/module.h>
#include <k3dsdk/node.h>
#include <k3dsdk/material_client.h>
#include <k3dsdk/measurement.h>
#include <k3dsdk/mesh_modifier.h>
#include <k3dsdk/node.h>
#include <k3dsdk/persistent.h>
#include <k3dsdk/utility.h>
#include <iostream>
#include <map>
#include "mesh_info.h"
#include "diff_geom.h"

namespace libk3dquadremesh
{
	class pgp_remesh :
		public k3d::material_client<k3d::mesh_modifier<k3d::persistent<k3d::node> > >
		//public k3d::material_client<k3d::mesh_modifier<k3d::node > >
		//public k3d::node
	{
		typedef  k3d::material_client<k3d::mesh_modifier<k3d::persistent<k3d::node> > > base;
	public:
		pgp_remesh(k3d::iplugin_factory& Factory, k3d::idocument& Document) :
			base(Factory, Document),
			m_mean_coord(init_owner(*this) + init_name("use_mean") + init_label(_("Use Mean Coordinate Weights")) + init_description(_("Use Mean Coordinate Weights")) + init_value(true)),
			m_smooth(init_owner(*this) + init_name("use_smooth") + init_label(_("Smooth Curvature")) + init_description(_("Smooth Curvature")) + init_value(false)),
			m_symmetry(init_owner(*this) + init_name("smooth_4") + init_label(_("Smooth as 4-symmetry")) + init_description(_("Smooth as 4-symmetry")) + init_value(false))
		{
			m_mean_coord.changed_signal().connect(make_reset_mesh_slot());
			m_smooth.changed_signal().connect(make_reset_mesh_slot());
			m_symmetry.changed_signal().connect(make_reset_mesh_slot());
		}

		~pgp_remesh()
		{
		}

		void on_create_mesh(const k3d::mesh& InputMesh, k3d::mesh& OutputMesh) 
		{
			//k3d::log() << debug << "PGP: create mesh: " << std::endl;
			//detail::mesh_info m(InputMesh); 
			//k3d::log() << debug << "PGP: create mesh: diff geom" << std::endl;
			//detail::diff_geom diff(m);
			OutputMesh = InputMesh;
			//k3d::typed_array<k3d::vector3>* curv_p = new k3d::typed_array<k3d::vector3>;
			//k3d::log() << debug << "PGP: create mesh: fill diff geom" << std::endl;
			//diff.fill_diff_geom(OutputMesh);
			//// Will do this more efficiently later
		}
		void on_update_mesh(const k3d::mesh& InputMesh, k3d::mesh& OutputMesh)		  
		{
			detail::mesh_info m(InputMesh); 
			detail::diff_geom diff(m);
			OutputMesh = InputMesh;
			k3d::typed_array<k3d::vector3>* curv_p = new k3d::typed_array<k3d::vector3>;
			diff.fill_diff_geom(OutputMesh, m_mean_coord.value(), m_smooth.value(), m_symmetry.value());
		}

		static k3d::iplugin_factory& get_factory()
		{
		static k3d::document_plugin_factory<pgp_remesh,
			k3d::interface_list<k3d::imesh_source,
			k3d::interface_list<k3d::imesh_sink > > > factory(
			  k3d::uuid(0xc97aa4ce, 0x412c1ed2, 0x055044a4, 0xa151f085),
			  "PGP Remesh",
			  _("Quad remeshing using the PGP algorithm"),
			  "PGP",
			  k3d::iplugin_factory::EXPERIMENTAL);
			
			return factory;

		}
	private:
		k3d_data(bool, immutable_name, change_signal, with_undo, local_storage, no_constraint, writable_property, with_serialization) m_mean_coord;
		k3d_data(bool, immutable_name, change_signal, with_undo, local_storage, no_constraint, writable_property, with_serialization) m_smooth;
		k3d_data(bool, immutable_name, change_signal, with_undo, local_storage, no_constraint, writable_property, with_serialization) m_symmetry;

	};


	k3d::iplugin_factory& pgp_remesh_factory()
	{
		return pgp_remesh::get_factory();
	}
} // namespace pgp_module
