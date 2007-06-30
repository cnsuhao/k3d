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
		double principal_curve_tensor2(vert_t vert, double radius, Vec3& curv_dir0,  Vec3& curv_dir1, Vec3& tens);		Vec3 isotropic_tensor(Vec3 tensor);
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
	};
};

}; // namespace pgp_module
#endif
