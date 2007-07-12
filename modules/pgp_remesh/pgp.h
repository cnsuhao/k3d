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

#ifndef PGP_PGP_H
#define PGP_PGP_H

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
#include <utility>
#include "mesh_info.h"
#include "diff_geom.h"
#include <gmm/gmm.h>

namespace libk3dquadremesh
{

namespace detail {
	/// TODO: Add desc
	typedef diff_geom::tensor_t vec2;
	
	inline vec2 operator+(const vec2& a, const vec2& b) {
		return vec2(a.first+b.first, a.second+b.second);
	}
	inline vec2 operator-(const vec2& a, const vec2& b) {
		return vec2(a.first-b.first, a.second-b.second);
	}
	inline double operator*(const vec2& a, const vec2& b) {
		return a.first*b.first + a.second*b.second;
	}
	inline vec2 operator*(double a, const vec2& b) {
		return vec2(a*b.first, a*b.second);
	}
	inline vec2 operator*(const vec2& a, double b) {
		return vec2(b*a.first, b*a.second);
	}
	inline vec2 rotate90(const vec2& a, int n) {
		switch(n%4) {
			case 0:
				return vec2(a);
			case 1:
				return vec2(-a.second, a.first);
			case 2:
				return vec2(-a.first, -a.second);
			default:
				return vec2(a.second, -a.first);
		}
	}

	
	class PGP
	{
	public:
		PGP() 
		{
		}

		PGP(const mesh_info* Mesh, const diff_geom* Geom) 
			: mesh(Mesh), geom(Geom)
		{
		}

		~PGP() 
		{

		}

		void setup() 
		{
			face_data.resize(mesh->num_faces);
			vert_data.resize(mesh->num_verts);
			

			for(vert_t v = 0; v < mesh->num_verts; ++v) {
				double angle = 0.5*atan2(geom->tensor[v].second , geom->tensor[v].first);

				Vec3 temp_i = geom->tangent_basis_i[v];
				Vec3 temp_j = geom->tangent_basis_j[v];

				temp_i *= cos(angle);
				temp_j *= sin(angle);
				vert_data[v].dir = temp_i + temp_j;
			}

			// Fill face data with edge and vert indices, and project vf into triangle
			for(face_t f = 0; f < mesh->num_faces; ++f) {
				mesh_info::Face face = mesh->getFace(f);

				per_face &pf = face_data[f];

				pf.edge[0] = face[0]();
				pf.edge[1] = face[1]();
				pf.edge[2] = face[2]();

				pf.vert[0] = mesh->edge_vert[pf.edge[2]];
				pf.vert[1] = mesh->edge_vert[pf.edge[0]];
				pf.vert[2] = mesh->edge_vert[pf.edge[1]];

				mesh_info::Vert v0 = mesh->getVert(pf.vert[0]);
				mesh_info::Vert v1 = mesh->getVert(pf.vert[1]);
				mesh_info::Vert v2 = mesh->getVert(pf.vert[2]);
				mesh_info::Edge e0 = mesh->getEdge(pf.edge[0]);
				mesh_info::Edge e1 = mesh->getEdge(pf.edge[1]);
				mesh_info::Edge e2 = mesh->getEdge(pf.edge[2]);

				Vec3 norm = triangle_normal(v0.pos(), v1.pos(), v2.pos());
				Vec3 i = e0.dir();
				i.Normalize();
				Vec3 j = norm^i;
				j.Normalize();

				pf.e[0].first = e0.dir()*i;
				pf.e[0].second = 0.0;

				pf.e[1].first = e1.dir()*i;
				pf.e[1].second = e1.dir()*j;

				pf.e[2].first = e2.dir()*i;
				pf.e[2].second = e2.dir()*j;

				Vec3 temp = vert_data[pf.vert[0]].dir;
				temp = (temp*norm)*norm;
				temp = vert_data[pf.vert[0]].dir - temp;
				pf.vf[0].first  = temp*i;
				pf.vf[0].second = temp*j;

				temp = vert_data[pf.vert[1]].dir;
				temp = (temp*norm)*norm;
				temp = vert_data[pf.vert[1]].dir - temp;
				pf.vf[1].first  = temp*i;
				pf.vf[1].second = temp*j;
	
				temp = vert_data[pf.vert[2]].dir;
				temp = (temp*norm)*norm;
				temp = vert_data[pf.vert[2]].dir - temp;
				pf.vf[2].first  = temp*i;
				pf.vf[2].second = temp*j;

				std::vector<double> x(3);
				std::vector<double> b(3);
				b[0] = 1.0;
				b[1] = 1.0;
				b[2] = 0.0;

				gmm::dense_matrix<double> L(3,3);
				L(0,0) = pf.e[0].first*pf.e[0].first;
				L(0,1) = pf.e[1].first*pf.e[1].first;
				L(0,2) = pf.e[2].first*pf.e[2].first;
				L(1,0) = pf.e[0].second*pf.e[0].second;
				L(1,1) = pf.e[1].second*pf.e[1].second;
				L(1,2) = pf.e[2].second*pf.e[2].second;
				L(2,0) = 2.0*pf.e[0].first*pf.e[0].second;
				L(2,1) = 2.0*pf.e[1].first*pf.e[1].second;
				L(2,2) = 2.0*pf.e[2].first*pf.e[2].second;

				gmm::lu_solve(L, x, b);

				pf.lamda[0] = x[0];
				pf.lamda[1] = x[1];
				pf.lamda[2] = x[2];

				std::cout << "t = " << f << " L = " << "(" << x[0] << ", " << x[1] << ", " << x[2] << ")\n";

				double max1 = pf.vf[0]*pf.vf[1]; 
				double max2 = pf.vf[0]*pf.vf[2]; 

				double max_arg1 = 0;
				double max_arg2 = 0;

				for(int i = 1; i < 4; i++) {
					double dot1 = pf.vf[0]*rotate90(pf.vf[1], i);
					double dot2 = pf.vf[0]*rotate90(pf.vf[2], i);

					if(dot1 > max1) {
						max_arg1 = i;
						max1 = dot1;
					}
					if(dot2 > max2) {
						max_arg2 = i;
						max2 = dot2;
					}
				}

				pf.delta[0] = omega*0.5*((pf.vf[1] + pf.vf[2])*pf.e[0]);
				pf.delta[1] = omega*0.5*((pf.vf[0] + pf.vf[2])*pf.e[1]);
				pf.delta[2] = omega*0.5*((pf.vf[0] + pf.vf[1])*pf.e[2]);
			}
		}

		void solve() 
		{
			//gmm::row_matrix< gmm::wsvector<double> > A(mesh->num_edges, mesh->num_verts);

			//// compute

			//gmm::clean(A, 1E-10);
			//gmm::csr_matrix<double> B(mesh->num_edges, mesh->num_verts);
			//gmm::copy(A, B);
			//gmm::row_matrix< gmm::wsvector<double> > 

			//gmm::row_matrix< gmm::wsvector<double> > C
			//

			gmm::row_matrix< gmm::wsvector<double> > M(3,3);
			M(0,0) += 1;
			M(1,1) += 1;
			M(2,2) += 1;
			k3d::log() << M << std::endl;
		}

		void extract() 
		{
			
		}

		void remesh(k3d::mesh& OutputMesh) 
		{
			
		}


	private:
		const mesh_info* mesh; 
		const diff_geom* geom;
		struct per_face {
			edge_t edge[3];
			vert_t vert[3];
			int rot[3];
			double lamda[3];
			double delta[3];
			vec2 vf[3]; // vector field projected into triangle plane
			vec2 e[3]; // edge vectors in triangle space
		};

		struct per_vert {
			per_vert() : aligned(false) {}
			bool aligned;
			Vec3 dir;

		};
		std::vector<per_face> face_data;
		std::vector<per_vert> vert_data;
		double omega;
	};
};

}; // namespace pgp_module
#endif
