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
#include <utility>
#include "mesh_info.h"

namespace libk3dquadremesh
{

namespace detail {
	#include <modules/qslim/MxMath.h>
	#include <modules/qslim/MxTriangle.h>
	Vec3 triangle_raw_normal(const Vec3& v1, const Vec3& v2, const Vec3& v3);
	double triangle_area(const Vec3& v1, const Vec3& v2, const Vec3& v3);
	Vec3 triangle_normal(const Vec3& v1, const Vec3& v2, const Vec3& v3);


	/// TODO: Add desc
	class diff_geom 
	{
	public:
		diff_geom() 
		{
			init = false;
			valid = false;
		}

		diff_geom(const mesh_info* Mesh) 
			: mesh(Mesh)
		{
			init = false;
			valid = true;
		}

		~diff_geom() 
		{
		}


		/// initialize vectors and precompute values
		void initialize();

		/// Sends curvature directions to output mesh for drawing
		void dump_draw_data(k3d::mesh& OutputMesh);

		/// Take smoothing steps with timestep h, stable for all values of h
		void smooth(double h, int steps, bool four_symm);
		
		/// Returns normal of tangent plane at vertex
		Vec3 normal(vert_t vert);
		
		/// Calculate principal curvature tensor, returns error from previously calculated gaussian curvature
		double principal_curve_tensor(vert_t vert, std::pair<double,double>& tens);
		
		/// Alternative principal curvature tensor
		double principal_curve_tensor2(vert_t vert, double radius, Vec3& curv_dir0,  Vec3& curv_dir1, Vec3& norm);		
		
		/// Converts symmetric tensor into isotropic tensor
		std::pair<double,double> isotropic_tensor(Vec3 tensor);
		
		/// Calculates eigenvalues of symmetric 2x2 matrix
		void eigen(Vec3 tensor, double& e1, double& e2) ;
		 
		/// Discrete gaussian curvature normal of vertex
		double gaussian_curvature(vert_t vert);
		
		/// Discrete mean curvature normal of vertex
		Vec3 mean_curvature(vert_t vert);

		/// Voronoi region of vertex on edge intersecting with a triangle
		double area_mixed(edge_t edge);

		/// Project the vector into the tangent plane at vertex vert
		Vec3 project(vert_t vert, const Vec3& x);
		
		/// Cotangent of the angle of the vertex opposite of edge
		double cotangent(edge_t edge);

		/// Mean weight of half edge
		double mean_weight(edge_t edge);

		/// Mean coordinates of half edge
		double mean_coord(edge_t e);		
		
		/// Mesh that is being worked with
		const mesh_info* mesh;
		
		/// Length of the bounding box diagonal
		double bb_diag;

		/// Stores precomputed edge cotangents
		std::vector<double> edge_cot;

		/// Stores precomputed gaussian curvature
		std::vector<double> gaus_curv;

		/// Stores precomputed mean coordinates
		std::vector<double> mean_coords;

		/// Stores precomputed mean weights
		std::vector<double> mean_weights;

		/// Max principal curvature
		std::vector<double> k1;

		/// Min principal curvature
		std::vector<double> k2;

		std::vector<bool> face_searched;

		/**
		 * Holds first column of each matrix in the isotropic curvature tensor field (curvature magnitudes have been normalize).
		 * Each vector is < cos(2*theta), sin(2*theta) >, where theta is the angle between the tensor and the
		 * x-axis in the local tangent plane.
		 */
		std::vector<std::pair<double,double> > tensor;

		bool mean, init, valid;
	
		/// X-axis basis for tangent plane at each vertex 
		std::vector<Vec3> tangent_basis_i;
		/// Y-axis basis for tangent plane at each vertex  
		std::vector<Vec3> tangent_basis_j;
		/// Z-axis basis for tangent plane at each vertex 
		std::vector<Vec3> tangent_basis_k;

		std::vector<Vec3> mean_curv;
		std::vector<Vec3> curv_tens;

		
		/// Is the point vd left of the plane specified by the oriented triangle with vertices va, vb, and vc? 
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
		
		/// Signed dihedral angle of edge ab between oriented triangles abc and bad
		double signed_angle(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
			Vec3 t1_n = triangle_normal(a,b,c);
			Vec3 t2_n = triangle_normal(b,a,d);

			double dot = t1_n*t2_n;

			if( dot <= -1.0 ) return k3d::pi();  
			if( dot >=  1.0 ) return 0;  
			double angle = acos(t1_n*t2_n);


			if(leftOf(a,b,c,d))
				angle *= -1.0;

			return angle;
		}
		/// Length of line segment va->vb within a sphere specified by center and radius
		double segment_sphere_intersection_length(const Vec3& va, const Vec3& vb, double radius, const Vec3& center, bool a_in, Vec3& intersection) {
			double c_a = (vb-va)*(vb-va);
			double c_b = (va-center)*(vb-va);
			double c_c = (va-center)*(va-center) - radius*radius;
			double root = c_b * c_b - 4.0 * c_a * c_c;
			if(root < 0) {
				if(a_in) intersection = vb;
				else intersection = va;
				return (vb-va).Length();
			}
			root = std::sqrt(root);
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
