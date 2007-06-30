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
#include <gmm/gmm.h>

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

		vert_i_basis.resize(mesh.num_verts);
		vert_j_basis.resize(mesh.num_verts);
		face_i_basis.resize(mesh.num_faces);
		face_j_basis.resize(mesh.num_faces);
		rep_x.resize(mesh.num_verts);
		rep_y.resize(mesh.num_verts);
		mean_coords.resize(mesh.num_edges);
		mean_weights.resize(mesh.num_edges);

		k3d::log() << debug << "Start fill" << std::endl;
		for(size_t i = 0; i < mesh.num_edges; ++i) {
			edge_cot[i] = cotangent(i);
		}
		
		for(size_t i = 0; i < mesh.num_edges; ++i) {
			mean_weights[i] = mean_weight(i);
		}

		for(size_t i = 0; i < mesh.num_edges; ++i) {
			mean_coords[i] = mean_coord(i);
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


		//smooth(2, 1, 50);

		for(size_t i = 0; i < mesh.num_verts; ++i) {
			double angle = 0.5*atan2(rep_y[i] , rep_x[i]);

			Vec3 temp_i = vert_i_basis[i];
			Vec3 temp_j = vert_j_basis[i];

			temp_i *= cos(angle);
			temp_j *= sin(angle);
			pc1 = temp_i + temp_j;

			temp_i = vert_i_basis[i];
			temp_j = vert_j_basis[i];
			temp_i *= cos(angle + k3d::pi_over_2());
			temp_j *= sin(angle + k3d::pi_over_2());
			pc2 = temp_i + temp_j;

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

	// N-symmetry smoothing
	void diff_geom::smooth(int n, double h, int steps) {
		//std::vector<double> temp_x(mesh.num_verts), temp_y(mesh.num_verts);
		//std::vector<double> cos_Nangle(mesh.num_edges);
		//std::vector<double> sin_Nangle(mesh.num_edges);
		//int count = steps;
		//double error = 1;
		//temp_x.resize(mesh.num_verts);
		//temp_y.resize(mesh.num_verts);

		//cos_Nangle.resize(mesh.num_edges);
		//sin_Nangle.resize(mesh.num_edges);

		//for(edge_t e = 0; e < mesh.num_edges; e++) {
		//	Vec3 a = vert_i_basis[mesh.edge_vert[e]];
		//	Vec3 b = vert_i_basis[mesh.edge_vert[mesh.edge_ccw[e]]];
		//	double angle = 0.5*acos(a*b);
		//	cos_Nangle[e] = cos(n*angle);
		//	sin_Nangle[e] = sin(n*angle);
		//}
		
		std::vector<double> B(2*mesh.num_verts);
		std::vector<double> x(2*mesh.num_verts);
		
		gmm::row_matrix< gmm::wsvector<double> > M1(2*mesh.num_verts, 2*mesh.num_verts);
		for(vert_t v = 0; v < mesh.num_verts; v++) {
			mesh_info::Edge e = mesh.getVert(v).edge();
			edge_t first = e();
			
			M1(2*v, 2*v) = 1.0+h;
			M1(2*v+1, 2*v+1) = 1.0+h;
			
			do {
				vert_t w = e.comp().vert()();
				Vec3 geodesic = e.dir();

				Vec3 a = vert_i_basis[v];
				Vec3 b = vert_i_basis[w];
				double angle_i = acos(a*geodesic);
				double angle_j = acos(b*geodesic);
				double angle = angle_j - angle_i;
				double c = cos(n*angle);
				double s = sin(n*angle);

				M1(2*v, 2*w) = M1(2*v+1, 2*w+1) = -h*mean_coords[e()]*c;
				M1(2*v, 2*w+1) = h*mean_coords[e()]*s;
				M1(2*v+1, 2*w) = -h*mean_coords[e()]*s;			
			
				e = e.comp().next();
			} while(e() != first);		
		}

		gmm::csr_matrix<double> M;
		gmm::clean(M1, 1E-12);
		gmm::copy(M1, M);

		for(int v = 0; v < mesh.num_verts; v++) {
			B[2*v]   = rep_x[v];
			B[2*v+1] = rep_y[v];
		}

	
		gmm::diagonal_precond<gmm::csr_matrix<double> > PR(M);


		for(int i = 0; i < steps; ++i) {
			gmm::iteration iter(1E-7);
			iter.set_noisy(1);
//			gmm::bicgstab(M, x, B, PR, iter);
			double delta = 0;
			gmm::gmres(M, x, B, PR, 100, iter);
			for(vert_t v = 0; v < 2*mesh.num_verts; v++) {
				delta += (B[v]-x[v])*(B[v]-x[v]);
			}
			B.swap(x);

			if(delta < 1E-10) break;
		}

		for(vert_t v = 0; v < mesh.num_verts; v++) {
			rep_x[v] = B[2*v];
			rep_y[v] = B[2*v+1];
		}
		//for(vert_t v = 0; v < mesh.num_verts; v++) {
		//	double x = rep_x[v];
		//	double y = rep_y[v];
		//	x *= (1+h);
		//	y *= (1+h);

		//	mesh_info::Edge e = mesh.getVert(v).edge();
		//	edge_t first = e();
		//	do {
		//		double x_j, y_j;

		//		x_j = rep_x[e.next().vert()()];
		//		y_j = rep_y[e.next().vert()()];
		//		
		//		x -= mean_coords[e()]*h*(cos_Nangle[e()]*x_j - sin_Nangle[e()]*y_j);
		//		y -= mean_coords[e()]*h*(sin_Nangle[e()]*x_j + cos_Nangle[e()]*y_j);

		//		e = e.comp().next();
		//	} while(e() != first);

		//	temp_x[v] = x;
		//	temp_y[v] = y;

		//	error += (x-rep_x[v])*(x-rep_x[v]) + (y-rep_y[v])*(y-rep_y[v]);
		//}
	}


	Vec3 diff_geom::normal(vert_t vert) 
	{
		Vec3 mc = mean_curv[vert];
		if(mc*mc < 1E-10) {
			mesh_info::Edge e = mesh.getVert(vert).edge();
			edge_t first = e();
			mc[0] = 0;
			mc[1] = 0;
			mc[2] = 0;
			int count = 0;
			do {
				mc += triangle_normal(e.face()[0].start(), e.face()[1].start(), e.face()[2].start());
				count ++;
				e = e.comp().next();
			} while(e() != first);
			mc *= (1.0/count);
		}
		return mc;	
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
	

	double diff_geom::principal_curve_tensor2(vert_t vert, double radius, Vec3& curv_dir0,  Vec3& curv_dir1, Vec3& tens)
	{
		//gmm::dense_matrix<double> tensor(3,3);
		//gmm::clear(tensor);
		//gmm::dense_matrix<double> outer(3,3);
		//gmm::clear(outer);
		//gmm::dense_matrix<double> orient(3,3);
		//
		//mesh_info::Vert v = mesh.getVert(vert);
		//mesh_info::Edge e = mesh.getEdge(v());

		//Vec3 center = v.pos();

		//double area = 0;
		//double angle;
		//double intersection;
		//double r_sqr = radius*radius;

		//std::vector<edge_t> search;

		//search.push_back(e());

		//for(int i = 0; i < search.size(); i++) {
		//	e = mesh.getEdge(search[i]);

		//	if(face_searched[e.face()()]) {
		//		face_searched[e.face()()] = true;
		//		// add face area;
		//		// also add all edge matrices that have not been added yet
		//		// add all edge companions adjacent to faces that have not been searched yet

		//		Vec3 a = e.start();
		//		Vec3 b = e.next().start();
		//		Vec3 c = e.next().next().start();
		//		
		//		bool a_in,b_in,c_in;

		//		a_in = (a-center)*(a-center) < r_sqr;
		//		b_in = (b-center)*(b-center) < r_sqr;
		//		c_in = (c-center)*(c-center) < r_sqr;

		//		if(a_in && b_in) {
		//			
		//		} else {
		//			
		//		}
		//	}
		//}
		//
		//for(int i = 0; i < search.size(); i++ {
		//	// remove faces as having been searched
		//}
	}

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
		mesh_info::Edge e = v.edge();
		edge_t first = e();
				
		bool print = (v.pos()[2] == 0.0);// || (v() % 100 == 0);
		
		ihat = e.dir();
		ihat.Normalize();
		
		khat = normal(vert);
		jhat = khat ^ ihat;
		ihat = jhat ^ khat;

		ihat.Normalize();
		jhat.Normalize();
		khat.Normalize();

		vert_i_basis[vert] = ihat;
		vert_j_basis[vert] = jhat;

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
			e = e.comp().next();
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

			//w = (cotangent(e()) + cotangent(e.comp()())) * len_square;
			w = mean_coords[e()];
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
			rep_x[vert] = 0;
			rep_y[vert] = 0;

			curv_dir0 = Vec3();
			curv_dir1 = Vec3();
			return 0;
		}
		det = 1/det;

		a = det*(AtA[2]*Atb[0] - AtA[1]*Atb[1]);
		b = det*(AtA[0]*Atb[1] - AtA[1]*Atb[0]);
		c = Kh2 - a;
		
		double error = gaus_curv[vert] - a*c + b*b;
		//if(error > 1) {
		//	b = sqrt(std::abs(a*c - gaus_curv[vert]));
		//}

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
		rep_x[vert] = t[0];
		rep_y[vert] = t[1];

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

	double diff_geom::mean_weight(edge_t edge) 
	{
		mesh_info::Edge e(mesh.getEdge(edge));

		Vec3 e_v = e.dir();
		Vec3 enn_v = e.next().next().dir();
		Vec3 ecn_v = e.comp().next().dir();

		double inv_len = 1.0/e_v.Length();

		e_v *= inv_len;
		enn_v.Normalize();
		ecn_v.Normalize();

		double cos_a0 = -(e_v * enn_v);
		double cos_a1 = e_v * ecn_v;
		
		//double tan_half_a0 = std::sqrt((1-cos_a0)/(1+cos_a0));
		//double tan_half_a1 = std::sqrt((1-cos_a1)/(1+cos_a1));
		
		// Which one is faster/more accurate?
		 double tan_half_a0 = tan(0.5*acos(cos_a0));
		 double tan_half_a1 = tan(0.5*acos(cos_a1));
		
		return inv_len*(tan_half_a0 + tan_half_a1);
	}

	double diff_geom::mean_coord(edge_t edge)
	{
		mesh_info::Edge e = mesh.getEdge(edge);
		edge_t first = e();

		double sum = 0;

		do {
			sum += mean_weights[e()];
			e = e.comp().next();
		} while(e() != first);

		return mean_weights[first]/sum;
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
