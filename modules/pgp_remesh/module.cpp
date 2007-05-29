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

#include <iostream>
#include <map>

namespace libk3dquadremesh
{

namespace detail {
	typedef k3d::mesh::indices_t indices_t;

	// TODO: Add desc
	void calc_edge_face_adj(const k3d::mesh::polyhedra_t& Polyhedra, k3d::typed_array<indices_t>& adj) {
		for(int i = 0; i < Polyhedra.face_first_loops->size(); ++i) {
			int edge = Polyhedra.loop_first_edges->at(Polyhedra.face_first_loops->at(i));
			int first = edge;
			do {
				adj[edge] = i;				
				edge = Polyhedra.clockwise_edges->at(edge);
			} while( edge != first );
		}
		// error checking?
	}		

	// TODO: Add desc
	void calc_edge_companion_adj(const k3d::mesh::polyhedra_t& Polyhedra, k3d::typed_array<indices_t>& adj) {
		std::map<std::pair<indices_t, indices_t>, indices_t> segments;

		for(int i = 0; i < Polyhedra.edge_points->size(); ++i) {
			int v0 = Polyhedra.edge_points->at(i);
			int v1 = Polyhedra.edge_points->at(Polyhedra.clockwise_edges(i));

			segments[std::pair<v0,v1>] = i;
		}

		for(int i = 0; i < Polyhedra.edge_points->size(); ++i) {
			int v0 = Polyhedra.edge_points->at(i);
			int v1 = Polyhedra.edge_points->at(Polyhedra.clockwise_edges(i));

			int comp = segments[std::pair<v1,v0>];

			adj[i] = comp;
		}
		// assert correctness
	}

	// TODO: Add desc	
	void calc_vert_edge_adj(const k3d::mesh::polyhedra_t& Polyhedra, k3d::typed_array<indices_t>& adj) {
		for(int i = 0; i < Polyhedra.edge_points->size(); ++i) {
			int vert = Polyhedra.edge_points->at(i);
			adj[vert] = i;				
		}
	}

	// TODO: Add desc
	void calc_face_poly_adj(const k3d::mesh::polyhedra_t& Polyhedra, k3d::typed_array<indices_t>& adj) {
		for(int i = 0; i < Polyhedra.face_first_loops->size(); ++i) {
			
		}
	}

	class mesh_info 
	{
		mesh_info(const k3d::mesh& Mesh) {
		}
		
		k3d::typed_array<indices_t> edge_comp;
		k3d::typed_array<indices_t> edge_face;
		k3d::typed_array<indices_t> vert_edge;
		k3d::typed_array<indices_t> face_poly;
		//k3d::typed_array<> gaussian_curv;
		//k3d::typed_array<> mean_curv;
	};
};

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
			k3d::log() << debug << "PGP: create mesh: " << InputMesh.polyhedra->first_faces->size() << std::endl;
			k3d::log() << debug << "PGP: create mesh: " << InputMesh.polyhedra->first_faces->size() << std::endl;
			for(int i = 0; i < InputMesh.polyhedra->first_faces->size(); i++) {
				
			}
		}
		void on_update_mesh(const k3d::mesh& InputMesh, k3d::mesh& OutputMesh)		  
		{
			k3d::log() << debug << "PGP: update mesh" << std::endl;
		}

		static k3d::iplugin_factory& get_factory()
		{
			static k3d::document_plugin_factory<pgp_remesh> factory(
			  k3d::uuid(0xc97aa4ce, 0x412c1ed2, 0x055044a4, 0xa151f085),
			  "PGP Remesh",
			  "Quad remeshing using the PGP algorithm",
			  "PGP",
			  k3d::iplugin_factory::EXPERIMENTAL);
			
			return factory;

		}
	};

} // namespace pgp_module

K3D_MODULE_START(Registry)
	Registry.register_factory(libk3dquadremesh::pgp_remesh::get_factory());
K3D_MODULE_END
