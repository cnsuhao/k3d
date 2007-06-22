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
#include "diff_geom.h"

namespace libk3dquadremesh
{

namespace detail {
	#include <modules/qslim/MxMath.h>
	#include <modules/qslim/MxTriangle.h>


	/// TODO: Add desc
	void diff_geom::fill_diff_geom(k3d::mesh& OutputMesh) 
	{
		k3d::typed_array<k3d::vector3>* curv_p = new k3d::typed_array<k3d::vector3>;
		boost::shared_ptr<k3d::typed_array<k3d::vector3> > curv(curv_p);

		k3d::typed_array<k3d::vector3>* p1_p = new k3d::typed_array<k3d::vector3>;
		boost::shared_ptr<k3d::typed_array<k3d::vector3> > p1(p1_p);

		k3d::typed_array<k3d::vector3>* p2_p = new k3d::typed_array<k3d::vector3>;
		boost::shared_ptr<k3d::typed_array<k3d::vector3> > p2(p2_p);

		k3d::typed_array<std::vector<edge_t> >* ring_p = new k3d::typed_array<std::vector<edge_t> >;
		boost::shared_ptr<k3d::typed_array<std::vector<edge_t> > > ring(ring_p);

		curv->resize(mesh.num_verts);
		p1->resize(mesh.num_verts);
		p2->resize(mesh.num_verts);
		ring->resize(mesh.num_verts);
		edge_cot.resize(mesh.num_edges);
		mean_curv.resize(mesh.num_edges);
		gaus_curv.resize(mesh.num_edges);

		k3d::log() << debug << "Start fill" << std::endl;
		for(size_t i = 0; i < mesh.num_edges; ++i) {
			edge_cot[i] = cotangent(i);
		}

		k3d::log() << debug << "Done cot" << std::endl;
		for(size_t i = 0; i < mesh.num_verts; ++i) {
			mean_curv[i] = mean_curvature(i);
			curv->at(i).n[0] = mean_curv[i][0];
			curv->at(i).n[1] = mean_curv[i][1];
			curv->at(i).n[2] = mean_curv[i][2];
		}
		k3d::log() << debug << "Done mean" << std::endl;
		for(size_t i = 0; i < mesh.num_verts; ++i) {
			gaus_curv[i] = gaussian_curvature(i);
		}
		avg_e1 = 0;
		avg_e2 = 0;
		Vec3 pc1,pc2, tens;
		k3d::log() << debug << "Done gaussian" << std::endl;
		for(size_t i = 0; i < mesh.num_verts; ++i) {
			principal_curve_tensor(i, pc1,  pc2, tens);
			p1->at(i).n[0] = pc1[0];
			p1->at(i).n[1] = pc1[1];
			p1->at(i).n[2] = pc1[2];
			
			p2->at(i).n[0] = pc2[0];
			p2->at(i).n[1] = pc2[1];
			p2->at(i).n[2] = pc2[2];
		}

		OutputMesh.vertex_data["PGPMeanCurv"] = curv;
		OutputMesh.vertex_data["PGPPrincCurv1"] = p1;
		OutputMesh.vertex_data["PGPPrincCurv2"] = p2;
		OutputMesh.vertex_data["PGPOneRing"] = ring;

		for(size_t i = 0; i < mesh.num_verts; ++i) {
			mesh_info::Vert v = mesh.getVert(i);
			mesh_info::Edge e = v.edge();
			edge_t first = e();
			do {
				ring->at(i).push_back(e.next().vert().index);
				e = e.comp().next();
			} while(e() != first);
		}
	}

	Vec3 diff_geom::normal(vert_t vert) 
	{
		return mean_curv[vert];	
	}

	//double diff_geom::principal_curve_tensor(vert_t vert, Vec3& curv_dir0,  Vec3& curv_dir1, Vec3& tens)
	//{
	//	// Choose basis for plane
	//	// For each edge adjacent to vert
	//	//		project edge into plane basis
	//	//		get mean curv for edge
	//	//		accumulate a and b coefficients
	//	//	Solve linear system to get a and b
	//	//  c = 2 * Kh - a
	//	//  check error = kg - (a*c - b*b)
	//	
	//	Vec3 ihat, jhat, khat;
	//	double d_x, d_y, w, k, len_square;
	//	double m[4] = {0,0,0,0};
	//	double x[2] = {0,0};
	//	double a,b,c;

	//	mesh_info::Vert v(mesh, vert);
	//	mesh_info::Edge e = v.edge();
	//	edge_t first = e();
	//	
	//	ihat = e.dir();
	//	ihat.Normalize();
	//	
	//	khat = normal(vert);
	//	jhat = khat ^ ihat;
	//	ihat = jhat ^ khat;

	//	ihat.Normalize();
	//	jhat.Normalize();
	//	khat.Normalize();
	//	Vec3 d;
	//	double len;

	//	do {
	//		d = khat;
	//		d *= (khat*e.dir());
	//		d = e.dir() - d;
	//		d.Normalize();

	//		d_x = d*ihat;
	//		d_y = d*jhat; 

	//		len = sqrt(d_x*d_x + d_y*d_y);

	//		d_x /= len;
	//		d_y /= len;
	//		
	//		len_square = e.dir() * e.dir(); 

	//		w = (edge_cot[e()] + edge_cot[e.comp()()]) * len_square;
	//		k = 2.0*(e.dir()*khat)/len_square;
	//		
	//		m[0] += w * d_x * d_x * d_x * d_x;
	//		m[1] += 2.0* w * d_x * d_x * d_x * d_y;
	//		m[2] += w * d_x * d_x * d_x * d_y;
	//		m[3] += w * d_x * d_x * d_y * d_y;

	//		x[0] += w * d_x * d_x * k;
	//		x[1] += w * d_x * d_y * k;

	//		e = e.comp().next();
	//	} while(e() != first);

	//	double det = m[0]*m[3] - m[1]*m[2];
	//	
	//	if(det == 0.0) {
	//		curv_dir0 = Vec3();
	//		curv_dir1 = Vec3();
	//		return 0;
	//	}
	//	det = 1/det;

	//	a = det*(m[3]*x[0] - m[1]*x[1]);
	//	b = det*(-1.0*m[2]*x[0] + m[0]*x[1]);
	//	c = mean_curv[vert].Length() - a;
	//	
	//	double error = gaus_curv[vert] - a*c + b*b;
	//	tens = Vec3(a,b,c);
	//	//Vec3 vi = Vec3(mesh.points->at(vert).n);
	//	//Vec3 vj = Vec3(mesh.points->at(vert_edge[edge_comp[edge]]).n);
	//	Vec3 iso = isotropic_tensor(tens);
	//	
	//	double *t = iso;
	//	double angle = 0.5*atan2(t[1], t[0]);
	//	double e1, e2;
	//	eigen(tens, e1, e2);

	//	Vec3 temp_i = ihat;
	//	Vec3 temp_j = jhat;

	//	temp_i *= e1*cos(angle);
	//	temp_j *= e1*sin(angle);
	//	curv_dir0 = temp_i + temp_j;

	//	temp_i = ihat;
	//	temp_j = jhat;
	//	temp_i *= e2*cos(angle + k3d::pi_over_2());
	//	temp_j *= e2*sin(angle + k3d::pi_over_2());
	//	curv_dir1 = temp_i + temp_j;

	//	return error;	
	//}
	

	double diff_geom::principal_curve_tensor(vert_t vert, Vec3& curv_dir0,  Vec3& curv_dir1, Vec3& tens)
	{
		// Choose basis for plane
		// For each edge adjacent to vert
		//		project edge into plane basis
		//		get mean curv for edge
		//		accumulate a and b coefficients
		//	Solve linear system to get a and b
		//  c = 2 * Kh - a
		//  check error = kg - (a*c - b*b)
		
		Vec3 ihat, jhat, khat;
		double d_x, d_y, w, k, len_square;
		double AtA[3] = {0,0,0}; // Symmetric 2x2
		double Atb[2] = {0,0}; 
		double a,b,c;
		double Kh2 = mean_curv[vert].Length();
		mesh_info::Vert v(mesh, vert);
		mesh_info::Edge e = v.edge().comp().next();
		edge_t first = e();
		
		bool print = (v.pos()[2] == 0.0) || (v() % 10 == 0);

		ihat = e.dir();
		ihat.Normalize();
		
		khat = normal(vert);
		jhat = khat ^ ihat;
		ihat = jhat ^ khat;

		ihat.Normalize();
		jhat.Normalize();
		khat.Normalize();

		if(print) {
			k3d::log() << "----" << v() << "----" << std::endl; 
			k3d::log() << "i = (" << ihat[0] << ", " << ihat[1] << ", " << ihat[2] << ")" << std::endl; 
			k3d::log() << "j = (" << jhat[0] << ", " << jhat[1] << ", " << jhat[2] << ")" << std::endl; 
			k3d::log() << "k = (" << khat[0] << ", " << khat[1] << ", " << khat[2] << ")" << std::endl; 
			k3d::log() << "Kh2 = (" << Kh2 << ")" << std::endl; 			
		}

		Vec3 d;

		double len;
		double area = 0;
		double valence = 0;

		do {
			area += area_mixed(e());
		} while(e() != first);

		area = 1.0/area;
		int n = 0;
		do {
			d = khat;
			d *= (khat*e.dir());
			d = e.dir() - d;
			d.Normalize();
			
			d_x = d*ihat;
			d_y = d*jhat; 
			
			len = sqrt(d_x*d_x + d_y*d_y);

			d_x /= len;
			d_y /= len;
			
			if(print)
				k3d::log() << n << ") d = (" << d_x << ", " << d_y << ")";
			n++;
			len_square = e.dir() * e.dir(); 

			w = area * (1.0/8.0) * (cotangent(e.index) + cotangent(e.comp().index)) * len_square;
			k = 2.0*(e.dir()*khat)/len_square;
			
			if(print)
				k3d::log()  << ", w = ("<<w<<")";

			double x2y2 = (d_x*d_x - d_y*d_y);
			double xy2  = 2.0*d_x*d_y;
			
			if(print)
				k3d::log() << ", k = (" << k << ")\n";
			
			AtA[0] += w * x2y2 * x2y2;
			AtA[1] += w * xy2 * x2y2;
			AtA[2] += w * xy2 * xy2;

			double sum = w * (k - Kh2 * d_y * d_y);
			Atb[0] += sum * x2y2;
			Atb[1] += sum * xy2;

			e = e.comp().next();
		} while(e() != first);

		if(print) {
			k3d::log() << "AtA = (" << AtA[0] << ", " << AtA[1] << ", " << AtA[2] << ")" << std::endl; 
			k3d::log() << "Atb = (" << Atb[0] << ", " << Atb[1] << ")" << std::endl; 
			k3d::log() << "K G = (" << gaus_curv[vert] << ")" << std::endl; 			
		}

		double det = AtA[0]*AtA[2] - AtA[1]*AtA[1];
		
		if(det == 0.0) {
			curv_dir0 = Vec3();
			curv_dir1 = Vec3();
			return 0;
		}
		det = 1/det;

		a = det*(AtA[2]*Atb[0] - AtA[1]*Atb[1]);
		b = det*(AtA[0]*Atb[1] - AtA[1]*Atb[0]);
		c = Kh2 - a;
		
		double error = gaus_curv[vert] - a*c + b*b;
		tens = Vec3(a,b,c);

		if(print) {
			k3d::log() << "Error = (" << error << ")" << std::endl; 
			k3d::log() << "abc = (" << a << ", " << b << ", " << c << ")" << std::endl;  
		}

		//Vec3 vi = Vec3(mesh.points->at(vert).n);
		//Vec3 vj = Vec3(mesh.points->at(vert_edge[edge_comp[edge]]).n);
		Vec3 iso = isotropic_tensor(tens);
		
		double *t = iso;
		double angle = 0.5*atan2(t[1], t[0]);
		double e1, e2;
		eigen(tens, e1, e2);

		if(print) {
			k3d::log() << "eigen = (" << e1 << ", " << e2 << ")" << std::endl;  
		}

		Vec3 temp_i = ihat;
		Vec3 temp_j = jhat;

		temp_i *= e1*cos(angle);
		temp_j *= e1*sin(angle);
		curv_dir0 = temp_i + temp_j;

		temp_i = ihat;
		temp_j = jhat;
		temp_i *= e2*cos(angle + k3d::pi_over_2());
		temp_j *= e2*sin(angle + k3d::pi_over_2());
		curv_dir1 = temp_i + temp_j;

		return error;	
	}
	
	Vec3 diff_geom::isotropic_tensor(Vec3 tensor) 
	{
		double *t = tensor;
		
		double lambda = 0.5*(t[0] + t[2]);
		Vec3 iso = Vec3(t[0] - lambda, t[1], 0);

		iso.Normalize();

		return iso;
	}

	void diff_geom::eigen(Vec3 tensor, double& e1, double& e2) 
	{
		double a = tensor[0];
		double b = tensor[1];
		double c = tensor[2];
		double root = sqrt(a*a + 4*b*b - 2*a*c + c*c );

		e1 = 0.5*(a + c - root);
		e2 = 0.5*(a + c + root);
	}

	double diff_geom::gaussian_curvature(vert_t vert) 
	{
		edge_t edge = mesh.vert_edge[vert];

		double curv = k3d::pi_times_2();
		double area = 0.0;
		do {
			area += area_mixed(edge);
			Vec3 a(mesh.mesh.points->at(mesh.edge_vert[edge]).n);
			Vec3 b(mesh.mesh.points->at(mesh.edge_vert[mesh.edge_ccw[edge]]).n);
			Vec3 c(mesh.mesh.points->at(mesh.edge_vert[mesh.edge_ccw[mesh.edge_ccw[edge]]]).n);
			
			Vec3 AB = b-a;
			Vec3 AC = c-a;
		
			AB.Normalize();
			AC.Normalize();

			curv -= acos(AB * AC);
			edge = mesh.edge_ccw[mesh.edge_comp[edge]];

		} while(edge != mesh.vert_edge[vert]);

		return curv/area;
	}

	Vec3 diff_geom::mean_curvature(vert_t vert) 
	{
		edge_t edge = mesh.vert_edge[vert];
		double area = 0.0;
		Vec3 mean(0,0,0);

		do {
			area += area_mixed(edge);
			
			double w = edge_cot[edge] + edge_cot[mesh.edge_comp[edge]];
			Vec3 a(mesh.mesh.points->at(mesh.edge_vert[edge]).n);
			Vec3 b(mesh.mesh.points->at(mesh.edge_vert[mesh.edge_ccw[edge]]).n);
			
			a -= b;
			a *= w;

			mean += a;

			edge = mesh.edge_ccw[mesh.edge_comp[edge]];
		} while(edge != mesh.vert_edge[vert]);
		
		mean /= (2.0*area);
		return mean;
	}

	/// Voronoi region of vertex on edge intersecting with a triangle
	double diff_geom::area_mixed(edge_t edge) 
	{
		Vec3 a(mesh.mesh.points->at(mesh.edge_vert[edge]).n);
		Vec3 b(mesh.mesh.points->at(mesh.edge_vert[mesh.edge_ccw[edge]]).n);
		Vec3 c(mesh.mesh.points->at(mesh.edge_vert[mesh.edge_ccw[mesh.edge_ccw[edge]]]).n);

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
		if(cosA < 0.0 || cosB <= 0.0 || cosC <= 0.0) {
			area = triangle_area(a,b,c);
			if(cosA <= 0.0)
				area /= 2.0;
			else
				area /= 4.0;
		} else {
			area = CA_norm2 * edge_cot[edge];
			area += AB_norm2 * edge_cot[mesh.edge_ccw[mesh.edge_ccw[edge]]];
			area *= 0.125;
		}

		return area;
	}

	/// Cotangent of the angle of the vertex opposite of edge
	double diff_geom::cotangent(edge_t edge) 
	{
		edge_t next = mesh.edge_ccw[edge];

		// get the 3d positions of the vertices
		Vec3 a(mesh.mesh.points->at(mesh.edge_vert[edge]).n);
		Vec3 b(mesh.mesh.points->at(mesh.edge_vert[mesh.edge_ccw[edge]]).n);
		Vec3 c(mesh.mesh.points->at(mesh.edge_vert[mesh.edge_ccw[mesh.edge_ccw[edge]]]).n);
		
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

};

}; // namespace pgp_module
