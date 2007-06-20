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
			base(Factory, Document)
		{
			k3d::log() << debug << "PGP Construct" << std::endl;
		}

		~pgp_remesh()
		{
			k3d::log() << debug << "PGP Deconstruct" << std::endl;
		}

		void on_create_mesh(const k3d::mesh& InputMesh, k3d::mesh& OutputMesh) 
		{
			detail::mesh_info m(InputMesh); 
			detail::diff_geom diff(m);
			OutputMesh = InputMesh;
			k3d::typed_array<k3d::vector3>* curv_p = new k3d::typed_array<k3d::vector3>;
			boost::shared_ptr<k3d::typed_array<k3d::vector3> > curv(curv_p);
			diff.fill_diff_geom(OutputMesh);
			curv->resize(OutputMesh.points->size());

			// Will do this more efficiently later
			for(int i = 0; i < curv->size(); i++) {
				curv->at(i).n[0] = diff.mean_curv[i][0];
				curv->at(i).n[1] = diff.mean_curv[i][1];
				curv->at(i).n[2] = diff.mean_curv[i][2];
			}

			OutputMesh.vertex_data["PGPMeanCurv"] = curv;

			k3d::log() << debug << "PGP: create mesh: " << curv.use_count() << " " << curv->size() << std::endl;
		}
		void on_update_mesh(const k3d::mesh& InputMesh, k3d::mesh& OutputMesh)		  
		{
			detail::mesh_info m(InputMesh); 
			detail::diff_geom diff(m);
			OutputMesh = InputMesh;
			k3d::typed_array<k3d::vector3>* curv_p = new k3d::typed_array<k3d::vector3>;
			boost::shared_ptr<k3d::typed_array<k3d::vector3> > curv(curv_p);
			diff.fill_diff_geom(OutputMesh);
			curv->resize(OutputMesh.points->size());

			// Will do this more efficiently later
			for(int i = 0; i < curv->size(); i++) {
				curv->at(i).n[0] = diff.mean_curv[i][0];
				curv->at(i).n[1] = diff.mean_curv[i][1];
				curv->at(i).n[2] = diff.mean_curv[i][2];
			}

			OutputMesh.vertex_data["PGPMeanCurv"] = curv;

			k3d::log() << debug << "PGP: update mesh" << std::endl;
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
	};


	k3d::iplugin_factory& pgp_remesh_factory()
	{
		return pgp_remesh::get_factory();
	}
} // namespace pgp_module
