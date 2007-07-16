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
#include <k3dsdk/color.h>
#include <iostream>
#include <fstream>
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

		void setup(double omega);
		void solve();

		void extract(double omega);

		void remesh(k3d::mesh& OutputMesh);

	private:
		const mesh_info* mesh; 
		const diff_geom* geom;
		struct per_face {
			edge_t edge[3];
			vert_t vert[3];
			int rot[3];
			double lamda[3];
			double delta[3];
			double delta_p[3];
			vec2 vf[3]; // vector field projected into triangle plane
			vec2 e[3]; // edge vectors in triangle space
			double theta[3];
			double phi[3];
		};

		struct per_vert {
			per_vert() : aligned(false) {}
			bool aligned;
			Vec3 dir;

			double theta;
			double phi;
		};
		std::vector<per_face> face_data;
		std::vector<per_vert> vert_data;

		gmm::dense_matrix<double> M[3];
	};
};

}; // namespace pgp_module
#endif
