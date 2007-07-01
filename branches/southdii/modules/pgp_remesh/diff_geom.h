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

#ifndef PGP_DIFF_GEOM_H
#define PGP_DIFF_GEOM_H

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

namespace libk3dquadremesh
{

namespace detail {


	/// TODO: Add desc
	class diff_geom 
	{
	public:
		diff_geom(const mesh_info& Mesh) 
			: mesh(Mesh)
		{
		}

		void fill_diff_geom(k3d::mesh& OutputMesh);
		void smooth(int n, double h, int steps);
		Vec3 normal(vert_t vert);
		double principal_curve_tensor(vert_t vert, Vec3& curv_dir0,  Vec3& curv_dir1, Vec3& tens);
		double principal_curve_tensor2(vert_t vert, double radius, Vec3& curv_dir0,  Vec3& curv_dir1, Vec3& norm);		Vec3 isotropic_tensor(Vec3 tensor);
		void eigen(Vec3 tensor, double& e1, double& e2) ;
		double gaussian_curvature(vert_t vert);
		Vec3 mean_curvature(vert_t vert);
		/// Voronoi region of vertex on edge intersecting with a triangle
		double area_mixed(edge_t edge);

		/// Cotangent of the angle of the vertex opposite of edge
		double cotangent(edge_t edge);
		double mean_weight(edge_t edge);
		double mean_coord(edge_t e);		
		const mesh_info& mesh;
		double bb_diag;
		std::vector<double> edge_cot;
		std::vector<double> gaus_curv;
		std::vector<double> mean_coords;
		std::vector<double> mean_weights;
		std::vector<double> k1;
		std::vector<double> k2;

		std::vector<bool> edge_searched;
		std::vector<bool> face_searched;

		std::vector<double> rep_x;
		std::vector<double> rep_y;

		std::vector<Vec3> face_i_basis;
		std::vector<Vec3> face_j_basis;
	
		std::vector<Vec3> vert_i_basis;
		std::vector<Vec3> vert_j_basis;

		std::vector<Vec3> mean_curv;
		std::vector<Vec3> curv_tens; // represents a b c values of tensor

		bool leftOf(const Vec3& va, const Vec3& vb, const Vec3& vc, const Vec3& vd) {
			Vec3 a(va);
			Vec3 b(va);
			Vec3 c(va);
			a -= vd;
			b -= vd;
			c -= vd;
			
			double det =   a[0]*(b[1]*c[2] - c[1]*b[2]) 
						 + a[1]*(b[2]*c[0] - c[2]*b[0]) 
						 + a[2]*(b[0]*c[1] - c[0]*b[1]);

			return det >= 0.0;
		}
		
		// Signed dihedral angle of edge ab between oriented triangles abc and bad
		double signed_angle(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
			Vec3 t1_n = triangle_normal(a,b,c);
			Vec3 t2_n = triangle_normal(b,a,d);

			double angle = acos(t1_n*t2_n);
			if(leftOf(a,b,c,d))
				angle *= -1.0;

			return angle;
		}
		// 
		double segment_sphere_intersection_length(const Vec3& va, const Vec3& vb, double radius, const Vec3& center, bool a_in, Vec3& intersection) {
			double c_a = (vb-va)*(vb-va);
			double c_b = (va-center)*(vb-va);
			double c_c = (va-center)*(va-center) - radius*radius;
			double root = std::sqrt(c_b * c_b - 4.0 * c_a * c_c);
			double t = (-c_b + root)/(2*c_a);
			Vec3 temp;
			if(t < 0) 
				t = (-c_b - root)/(2*c_a);
			intersection = (vb-va);
			intersection *= t;
			intersection += va;

			if(a_in) {
				
				return (intersection-va).Length();
			} else {
				return (intersection-vb).Length();
			}
		}
	};
};

}; // namespace pgp_module
#endif
