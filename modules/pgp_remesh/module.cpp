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
#include <modules/qslim/MxMath.h>
#include <modules/qslim/MxTriangle.h>
#include <k3dsdk/utility.h>
#include <iostream>
#include <map>

namespace libk3dquadremesh
{

namespace detail {
	#include <modules/qslim/MxMath.h>
	#include <modules/qslim/MxTriangle.h>
	typedef k3d::mesh::indices_t indices_t;
	typedef size_t vert_t;
	typedef size_t edge_t;
	typedef size_t face_t;
	typedef size_t poly_t;

	/// TODO
	void calc_edge_face_adj(const k3d::mesh::polyhedra_t& Polyhedra, indices_t& adj) {
		for(size_t i = 0; i < Polyhedra.face_first_loops->size(); ++i) {
			size_t edge = Polyhedra.loop_first_edges->at(Polyhedra.face_first_loops->at(i));
			size_t first = edge;
			do {
				adj[edge] = i;				
				edge = Polyhedra.clockwise_edges->at(edge);
			} while( edge != first );
		}
		// error checking?
	}		

	/// TODO: Add desc
	void calc_edge_ccw_adj(const k3d::mesh::polyhedra_t& Polyhedra, indices_t& adj) {
		for(size_t i = 0; i < Polyhedra.loop_first_edges->size(); ++i) {
			size_t edge = Polyhedra.loop_first_edges->at(i);
			size_t first = edge;
			do {
				adj[Polyhedra.clockwise_edges->at(edge)] = edge;				
				edge = Polyhedra.clockwise_edges->at(edge);
			} while( edge != first );
		}
		// error checking?
	}		

	/// TODO: Add desc
	void calc_edge_companion_adj(const k3d::mesh::polyhedra_t& Polyhedra, indices_t& adj) {
		std::map<std::pair<size_t, size_t>, size_t> segments;

		for(size_t i = 0; i < Polyhedra.edge_points->size(); ++i) {
			size_t v0 = Polyhedra.edge_points->at(i);
			size_t v1 = Polyhedra.edge_points->at(Polyhedra.clockwise_edges->at(i));

			segments[std::pair<size_t, size_t>(v0,v1)] = i;
		}

		for(size_t i = 0; i < Polyhedra.edge_points->size(); ++i) {
			size_t v0 = Polyhedra.edge_points->at(i);
			size_t v1 = Polyhedra.edge_points->at(Polyhedra.clockwise_edges->at(i));

			size_t comp = segments[std::pair<size_t, size_t>(v1,v0)];

			adj[i] = comp;
		}
		// assert correctness
	}

	/// TODO: Add desc	
	void calc_vert_edge_adj(const k3d::mesh::polyhedra_t& Polyhedra, indices_t& adj) {
		for(size_t i = 0; i < Polyhedra.edge_points->size(); ++i) {
			size_t vert = Polyhedra.edge_points->at(i);
			adj.at(vert) = i;
		}
	}

	/// TODO: Add desc	
	void calc_vert_edge_ccw_adj(const k3d::mesh::polyhedra_t& Polyhedra, indices_t& adj) {
		indices_t comp;
		comp.resize(Polyhedra.edge_points->size());
		calc_edge_companion_adj(Polyhedra, comp);
		for(size_t i = 0; i < Polyhedra.edge_points->size(); ++i) {
			size_t vert = Polyhedra.edge_points->at(comp.at(i));
			adj.at(vert) = i;
		}
	}

	/// TODO: Add desc	
	void calc_edge_vert_ccw_adj(const k3d::mesh::polyhedra_t& Polyhedra, indices_t& adj) {
		for(size_t i = 0; i < Polyhedra.edge_points->size(); ++i) {
	//		adj[i] = 
		}
	}

	/// TODO: Add desc
	void calc_face_poly_adj(const k3d::mesh::polyhedra_t& Polyhedra, indices_t& adj) {
		indices_t comp;
		indices_t edge_face;
		
		comp.resize(Polyhedra.edge_points->size());
		edge_face.resize(Polyhedra.edge_points->size());

		calc_edge_companion_adj(Polyhedra, comp);
		calc_edge_face_adj(Polyhedra, edge_face);

		size_t empty = 0xFFFFFFFF;
		
		std::fill(adj.begin(), adj.end(), empty);
		std::vector<size_t> stack;

		// Depth first traversal of faces in polyhedra to label all the faces
		for(size_t i = 0; i < Polyhedra.first_faces->size(); ++i) {
			stack.push_back(Polyhedra.first_faces->at(i));
			while(!stack.empty()) {
				size_t f = stack.back();
				stack.pop_back();
				if(adj[f] != empty) {
					adj[f] = i;
					size_t edge = Polyhedra.loop_first_edges->at(Polyhedra.face_first_loops->at(f));
					size_t first = edge;

					do {
						stack.push_back(edge_face[comp[edge]]);
						edge = Polyhedra.clockwise_edges->at(edge);
					} while( edge != first );
				}
			}
		}
	}


	/// TODO: Add desc
	class mesh_info 
	{
	public:
		mesh_info(const k3d::mesh& Mesh) 
			: mesh(Mesh)
		{
			edge_comp.resize(Mesh.polyhedra->edge_points->size());
			edge_ccw.resize(Mesh.polyhedra->edge_points->size());
			edge_face.resize(Mesh.polyhedra->edge_points->size());

			face_poly.resize(Mesh.polyhedra->edge_points->size());
			vert_edge.resize(Mesh.points->size());

			mean_curv.resize(Mesh.points->size());
			curv_tens.resize(Mesh.points->size());
			gaus_curv.resize(Mesh.points->size());

			calc_edge_face_adj(*mesh.polyhedra, edge_face);
			calc_edge_ccw_adj(*mesh.polyhedra, edge_ccw);
			calc_edge_companion_adj(*mesh.polyhedra, edge_comp);
			calc_vert_edge_ccw_adj(*mesh.polyhedra, vert_edge);
			calc_edge_vert_ccw_adj(*mesh.polyhedra, vert_edge);
			calc_face_poly_adj(*mesh.polyhedra, face_poly);
		}

		void fill_diff_geom() {
			for(size_t i = 0; i < num_edges; ++i) {
				edge_cot[i] = cotangent(i);
			}
			for(size_t i = 0; i < num_verts; ++i) {
				mean_curv[i] = mean_curvature(i);
			}
		}

		Vec3 mean_curvature(vert_t vert) {
			edge_t edge = vert_edge[vert];
			double area = 0.0;
			Vec3 mean(0,0,0);

			do {
				area += area_mixed(edge_face[edge], edge);
				
				double w = edge_cot[edge] + edge_cot[edge_comp[edge]];
				Vec3 a(mesh.points->at(edge_vert[edge]).n);
				Vec3 b(mesh.points->at(edge_vert[edge_ccw[edge]]).n);
				
				a -= b;
				a *= w;

				mean += a;

				edge = edge_ccw[edge_comp[edge]];
			} while(edge != vert_edge[vert]);
			
			mean /= (2.0*area);
			return mean;
		}

		// Assumes edge_cot is filled
		double area_mixed(face_t face, edge_t edge) {
			Vec3 a(mesh.points->at(edge_vert[edge]).n);
			Vec3 b(mesh.points->at(edge_vert[edge_ccw[edge]]).n);
			Vec3 c(mesh.points->at(edge_vert[edge_ccw[edge_ccw[edge]]]).n);

			Vec3 AB = b-a;
			Vec3 BC = c-b;
			Vec3 CA = a-c;
			
			double AB_norm2 = norm2(AB);
			double CA_norm2 = norm2(CA);

			AB.Normalize();
			BC.Normalize();
			CA.Normalize();

			double cosA = - (AB * CA);
			double cosB = - (BC * AB);
			double cosC = - (CA * BC);

			double area = 0;

			// cosine is <= 0 if angle is obtuse/right ?
			if(cosA <= 0.0 || cosB <= 0.0 || cosC <= 0.0) {
				area = triangle_area(a,b,c);
				if(cosA <= 0.0)
					area /= 2.0;
				else
					area /= 4.0;
			} else {
				area = CA_norm2 * edge_cot[edge];
				area += AB_norm2 * edge_cot[edge_ccw[edge_ccw[edge]]];
				area *= 0.125;
			}

			return area;
		}

		// Returns cotangent of the angle of the vertex opposite of edge
		double cotangent(edge_t edge) {
			edge_t next = edge_ccw[edge];

			// get the 3d positions of the vertices
			Vec3 a(mesh.points->at(edge_vert[edge]).n);
			Vec3 b(mesh.points->at(edge_vert[edge_ccw[edge]]).n);
			Vec3 c(mesh.points->at(edge_vert[edge_ccw[edge_ccw[edge]]]).n);
			
			// a and b become the two edge vectors originating from vertex c
			a -= c;
			b -= c;

			// We can easily find the cosine of the angle, and the cotangent can be
			// expressed in terms of cosine. cos/sqrt(1-cos^2)
			// Can be further simplified by getting rid of the need to find 
			// the norms of a and b and removing 2 divisions 
			double ab = a * b;
			double aa = a * a;
			double bb = b * b;
			
			// Divison by zero would only occur for colinear triangles
			double denom2 = aa * bb - ab * ab;
			
			if(denom2 == 0.0) return 0.0;

			return ab/(sqrt(denom2));
		}

		size_t num_edges;
		size_t num_faces;
		size_t num_verts;

		const k3d::mesh& mesh;

		indices_t edge_comp;

		indices_t edge_face;
		indices_t face_edge;

		indices_t vert_edge;
		indices_t edge_vert;

		indices_t edge_ccw;

		indices_t face_poly;
		indices_t poly_face;

		std::vector<double> edge_cot;
		std::vector<double> gaus_curv;
		std::vector<Vec3> mean_curv;
		std::vector<Vec3> normal;
		std::vector<Vec3> curv_tens; // represents a b c values of tensor
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
			detail::mesh_info m(InputMesh); 
			k3d::log() << debug << "PGP: create mesh: " << InputMesh.polyhedra->first_faces->size() << std::endl;
			k3d::log() << debug << "PGP: create mesh: " << InputMesh.polyhedra->first_faces->size() << std::endl;

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
