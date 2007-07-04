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
	Vec3 triangle_raw_normal(const Vec3& v1, const Vec3& v2, const Vec3& v3)
	{
		Vec3 a = v2 - v1;
		Vec3 b = v3 - v1;
		return a^b;
	}

	double triangle_area(const Vec3& v1, const Vec3& v2, const Vec3& v3)
	{
		return 0.5 * norm(triangle_raw_normal(v1, v2, v3));
	}

	
	Vec3 triangle_normal(const Vec3& v1, const Vec3& v2, const Vec3& v3)
	{
		Vec3 n = triangle_raw_normal(v1, v2, v3);
		n.Normalize();
		
		return n;
	}


	
	void diff_geom::initialize() {
		edge_cot.resize(mesh->num_edges);
		mean_curv.resize(mesh->num_edges);
		gaus_curv.resize(mesh->num_edges);
		face_searched.resize(mesh->num_faces, false);

		k1.resize(mesh->num_verts);
		k2.resize(mesh->num_verts);
		tangent_basis_i.resize(mesh->num_verts);
		tangent_basis_j.resize(mesh->num_verts);
		tangent_basis_k.resize(mesh->num_verts);
		tensor.resize(mesh->num_verts);
		mean_coords.resize(mesh->num_edges);
		mean_weights.resize(mesh->num_edges);

		Vec3 bb_min, bb_max;
		bb_min[0] = bb_max[0] = mesh->getVert(0).pos()[0];
		bb_min[1] = bb_max[1] = mesh->getVert(0).pos()[1];
		bb_min[2] = bb_max[2] = mesh->getVert(0).pos()[2];

		for(size_t i = 1; i < mesh->num_verts; ++i) {
			mesh_info::Vert v = mesh->getVert(i);
			if(v.pos()[0] > bb_max[0]) bb_max[0] = v.pos()[0];
			if(v.pos()[1] > bb_max[1]) bb_max[1] = v.pos()[1];
			if(v.pos()[2] > bb_max[2]) bb_max[2] = v.pos()[2];

			if(v.pos()[0] < bb_min[0]) bb_min[0] = v.pos()[0];
			if(v.pos()[1] < bb_min[1]) bb_min[1] = v.pos()[1];
			if(v.pos()[2] < bb_min[2]) bb_min[2] = v.pos()[2];
		}

		bb_diag = (bb_min-bb_max).Length();

//		std::cout << " cot \n";
		for(size_t i = 0; i < mesh->num_edges; ++i) {
			edge_cot[i] = cotangent(i);
		}
		
//		std::cout << " gauss \n";
		for(size_t i = 0; i < mesh->num_verts; ++i) {
			gaus_curv[i] = gaussian_curvature(i);
		}

//		std::cout << " weight \n";
		for(size_t i = 0; i < mesh->num_edges; ++i) {
			mean_weights[i] = mean_weight(i);
		}

//		std::cout << " coordt \n";
		for(size_t i = 0; i < mesh->num_edges; ++i) {
			mean_coords[i] = mean_coord(i);
		}

//		std::cout << " mean curv \n";
		for(size_t i = 0; i < mesh->num_verts; ++i) {
			mean_curv[i] = mean_curvature(i);
		}

//		std::cout << " basis \n";
		for(size_t i = 0; i < mesh->num_verts; ++i) {
			mesh_info::Vert v = mesh->getVert(i);
			mesh_info::Edge e = v.edge();
	
			Vec3 ihat = e.dir();
			ihat.Normalize();
			
			Vec3 khat = normal(i);
			Vec3 jhat = khat ^ ihat;
			ihat = jhat ^ khat;

			ihat.Normalize();
			jhat.Normalize();
			khat.Normalize();

			tangent_basis_i[i] = ihat;
			tangent_basis_j[i] = jhat;
			tangent_basis_k[i] = khat;
		}


//		std::cout << " pc \n";
		for(size_t i = 0; i < mesh->num_verts; ++i) {
			principal_curve_tensor(i, tensor[i]);
		}

//		std::cout << " end \n";

	}
	
	
	void diff_geom::dump_draw_data(k3d::mesh& OutputMesh) 
	{
		k3d::typed_array<k3d::vector3>* curv_p = new k3d::typed_array<k3d::vector3>;
		boost::shared_ptr<k3d::typed_array<k3d::vector3> > curv(curv_p);
		k3d::typed_array<k3d::vector3>* p1_p = new k3d::typed_array<k3d::vector3>;
		boost::shared_ptr<k3d::typed_array<k3d::vector3> > p1(p1_p);

		k3d::typed_array<k3d::vector3>* p2_p = new k3d::typed_array<k3d::vector3>;
		boost::shared_ptr<k3d::typed_array<k3d::vector3> > p2(p2_p);

		curv->resize(mesh->num_verts);
		p1->resize(mesh->num_verts);
		p2->resize(mesh->num_verts);


		Vec3 pc1,pc2, tens;
//		std::cout << " dump 1 \n";
		for(size_t i = 0; i < mesh->num_verts; ++i) {
			if(std::abs(tensor[i].first) <= 1E-10 && std::abs(tensor[i].second) <= 1E-10) {
				p1->at(i).n[0] = 0;
				p1->at(i).n[1] = 0;
				p1->at(i).n[2] = 0;
				
				p2->at(i).n[0] = 0;
				p2->at(i).n[1] = 0;
				p2->at(i).n[2] = 0;

				curv->at(i).n[0] = mean_curv[i][0];
				curv->at(i).n[1] = mean_curv[i][1];
				curv->at(i).n[2] = mean_curv[i][2];
				continue;
			}
			double angle = 0.5*atan2(tensor[i].second , tensor[i].first);

			Vec3 temp_i = tangent_basis_i[i];
			Vec3 temp_j = tangent_basis_j[i];

			temp_i *= cos(angle);
			temp_j *= sin(angle);
			pc1 = temp_i + temp_j;

			temp_i = tangent_basis_i[i];
			temp_j = tangent_basis_j[i];
			temp_i *= cos(angle + k3d::pi_over_2());
			temp_j *= sin(angle + k3d::pi_over_2());
			pc2 = temp_i + temp_j;
			p1->at(i).n[0] = pc1[0];
			p1->at(i).n[1] = pc1[1];
			p1->at(i).n[2] = pc1[2];
			
			p2->at(i).n[0] = pc2[0];
			p2->at(i).n[1] = pc2[1];
			p2->at(i).n[2] = pc2[2];

			curv->at(i).n[0] = mean_curv[i][0];
			curv->at(i).n[1] = mean_curv[i][1];
			curv->at(i).n[2] = mean_curv[i][2];
		}

		OutputMesh.vertex_data["PGPMeanCurv"] = curv;
		OutputMesh.vertex_data["PGPPrincCurv1"] = p1;
		OutputMesh.vertex_data["PGPPrincCurv2"] = p2;
	}

	void diff_geom::smooth(double h, int steps, bool four_symm = false) {
		//std::vector<double> temp_x(mesh->num_verts), temp_y(mesh->num_verts);
		//std::vector<double> cos_Nangle(mesh->num_edges);
		//std::vector<double> sin_Nangle(mesh->num_edges);
		//int count = steps;
		//double error = 1;
		//temp_x.resize(mesh->num_verts);
		//temp_y.resize(mesh->num_verts);

		//cos_Nangle.resize(mesh->num_edges);
		//sin_Nangle.resize(mesh->num_edges);

		//for(edge_t e = 0; e < mesh->num_edges; e++) {
		//	Vec3 a = vert_i_basis[mesh->edge_vert[e]];
		//	Vec3 b = vert_i_basis[mesh->edge_vert[mesh->edge_ccw[e]]];
		//	double angle = 0.5*acos(a*b);
		//	cos_Nangle[e] = cos(n*angle);
		//	sin_Nangle[e] = sin(n*angle);
		//}
		if(steps == 0) return;
		double n = four_symm ? 4 : 2;		
		std::vector<double> B(2*mesh->num_verts);
		std::vector<double> x(2*mesh->num_verts);
		
		
		gmm::row_matrix< gmm::wsvector<double> > M1(2*mesh->num_verts, 2*mesh->num_verts);
		for(vert_t v = 0; v < mesh->num_verts; v++) {
			mesh_info::Edge e = mesh->getVert(v).edge();
			edge_t first = e();
			
			M1(2*v, 2*v) = 1.0+h;
			M1(2*v+1, 2*v+1) = 1.0+h;
			
			do {
				vert_t w = e.comp().vert()();

				Vec3 g_v = project(v, e.dir());
				Vec3 g_w = project(w, e.dir());

				double angle_i = atan2(g_v[1], g_v[0]);
				double angle_j = atan2(g_w[1], g_w[0]);
				while(angle_i < 0) angle_i += k3d::pi_times_2();
				while(angle_j < 0) angle_j += k3d::pi_times_2();
				double angle = angle_i - angle_j;
				while(angle < 0) angle += k3d::pi_times_2();

				double c = cos(n*angle);
				double s = sin(n*angle);

				M1(2*v, 2*w) = M1(2*v+1, 2*w+1) = -h*mean_coords[e()]*c;
				M1(2*v, 2*w+1) = h*mean_coords[e()]*s;
				M1(2*v+1, 2*w) = -h*mean_coords[e()]*s;			
			
				e = e.comp().next();
			} while(e() != first);		
		}

		gmm::csr_matrix<double> M; // ! pre compute this!@!!!
		gmm::clean(M1, 1E-8);
		gmm::copy(M1, M);

		if(four_symm) {
			for(int v = 0; v < mesh->num_verts; v++) {
				double angle = (0.5)*atan2(tensor[v].second, tensor[v].first);
				B[2*v]   = cos(n*angle);
				B[2*v+1] = sin(n*angle);
			} 
		} else {
			for(int v = 0; v < mesh->num_verts; v++) {
				B[2*v]   = tensor[v].first;
				B[2*v+1] = tensor[v].second;
			} 
		}


	
		gmm::diagonal_precond<gmm::csr_matrix<double> > PR(M);

		int i = 0;
		for(i = 0; i < steps; ++i) {
			gmm::iteration iter(1E-10);
			iter.set_noisy(0);
			//iter.set_maxiter(150);
			gmm::bicgstab(M, x, B, PR, iter);
			double delta = 0;
//			gmm::gmres(M, x, B, PR, 100, iter);
			for(vert_t v = 0; v < 2*mesh->num_verts; v++) {
				delta += (B[v]-x[v])*(B[v]-x[v]);
			}
			B.swap(x);

			if(delta < 1E-10) break;
		}
		std::cout << i << "/" << steps << std::endl;

		for(vert_t v = 0; v < mesh->num_verts; v++) {
			double angle = (1.0/n)*atan2(B[2*v+1], B[2*v]);
			tensor[v].first  = cos(2*angle);
			tensor[v].second = sin(2*angle);
		}

	}


	Vec3 diff_geom::normal(vert_t vert) 
	{
		Vec3 mc = mean_curv[vert];
		Vec3 n;

		mesh_info::Edge e = mesh->getVert(vert).edge();
		edge_t first = e();
		n[0] = 0;
		n[1] = 0;
		n[2] = 0;
		int count = 0;
		do {
			n += triangle_normal(e.face()[0].start(), e.face()[1].start(), e.face()[2].start());
			count ++;
			e = e.comp().next();
		} while(e() != first);
		n *= (1.0/count);

		if(mc*mc < 1E-8) {
			n.Normalize();		
			return n;
		}
		if(mc*n < 0) {
			mc *= -1.0;
		}
		mc.Normalize();
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
	//	//Vec3 vi = Vec3(mesh->points->at(vert).n);
	//	//Vec3 vj = Vec3(mesh->points->at(vert_edge[edge_comp[edge]]).n);
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
	
	Vec3 diff_geom::project(vert_t vert, const Vec3& x) {
		Vec3 dir = tangent_basis_k[vert];
		dir *= (tangent_basis_k[vert]*x);
		dir = x - dir;
		return Vec3(dir*tangent_basis_i[vert], dir*tangent_basis_j[vert], 0);
	}

	double diff_geom::principal_curve_tensor2(vert_t vert, double radius, Vec3& curv_dir0,  Vec3& curv_dir1, Vec3& norm)
	{
////		std::cout << radius;
//		gmm::dense_matrix<double> tensor(3,3);
//		gmm::clear(tensor);
//		
//		mesh_info::Vert v = mesh->getVert(vert);
//		mesh_info::Edge e = mesh->getEdge(v());
//
//		Vec3 center = v.pos();
//
//		double area = 0;
//		double r_sqr = radius*radius;
//
//		std::vector<edge_t> search;
//
//		search.push_back(e());
//		
////		k3d::log() << debug << "start loop"  << std::endl;
//
//		for(int i = 0; i < search.size(); i++) {
//			e = mesh->getEdge(search[i]);
//			if(!face_searched[e.face()()]) {
//				face_searched[e.face()()] = true;
//				// add face area;
//				// also add all edge matrices that have not been added yet
//				// add all edge companions adjacent to faces that have not been searched yet
//
//				Vec3 va = e.start();
//				Vec3 vb = e.next().start();
//				Vec3 vc = e.next().next().start();
//				
//				double p;
//				bool a_in,b_in,c_in;
//				double length,angle;
//				a_in = (va-center)*(va-center) < r_sqr;
//				b_in = (vb-center)*(vb-center) < r_sqr;
//				c_in = (vc-center)*(vc-center) < r_sqr;
//				Vec3 edge;
//				int j = 0;
//				Vec3 intersection[2];
//				double intersection_area;
//
//				if(!face_searched[e.comp().face()()]) {
////					k3d::log() << debug << i << " : check e1"  << std::endl;
//					search.push_back(e.comp()());
//					Vec3 edge = e.dir();
//					if(a_in && b_in) {
//						length = edge.Length();
//					} else if(a_in || b_in) {
//						length = segment_sphere_intersection_length(va, vb, radius, center, a_in, intersection[j]);
//						j++;
//					}
//					angle = signed_angle(va, vb, vc, e.comp().next().end());
//					p = angle*length;
////					k3d::log() << debug << i << " : check e1 " << angle << " " << length << std::endl;
//
//					tensor(0,0) += p*edge[0]*edge[0];
//					tensor(1,1) += p*edge[1]*edge[1];
//					tensor(2,2) += p*edge[2]*edge[2];
//
//					tensor(0,1) += p*edge[0]*edge[1];
//					tensor(0,2) += p*edge[0]*edge[2];
//					tensor(1,2) += p*edge[1]*edge[2];
//				} else if((a_in && !b_in) || (!a_in && b_in)) {
//					segment_sphere_intersection_length(va, vb, radius, center, a_in, intersection[j]);
//					j++;
//				}
//
//				if(!face_searched[e.next().comp().face()()]) {
//					search.push_back(e.next().comp()());
//					Vec3 edge = e.next().dir();
//					if(b_in && c_in) {
//						length = edge.Length();
//					} else if(b_in || c_in) {
//						length = segment_sphere_intersection_length(vb, vc, radius, center, b_in, intersection[j]);
//						j++;
//					}
//					angle = signed_angle(vb, vc, va, e.next().comp().next().end());
////					k3d::log() << debug << i << " : check e2 " << angle << " " << length << std::endl;
//					p = angle*length;
//					tensor(0,0) += p*edge[0]*edge[0];
//					tensor(1,1) += p*edge[1]*edge[1];
//					tensor(2,2) += p*edge[2]*edge[2];
//
//					tensor(0,1) += p*edge[0]*edge[1];
//					tensor(0,2) += p*edge[0]*edge[2];
//					tensor(1,2) += p*edge[1]*edge[2];
//				} else if((b_in && !c_in) || (!b_in && c_in)) {
//					segment_sphere_intersection_length(vb, vc, radius, center, b_in, intersection[j]);
//					j++;
//				}
//
//				if(!face_searched[e.next().next().comp().face()()]) {
//					search.push_back(e.next().next().comp()());
//					Vec3 edge = e.next().next().dir();
//					if(c_in && a_in) {
//						length = edge.Length();
//					} else if(c_in || a_in) {
//						// j should not be 2 here, because that would imply that a_in && c_in is true
//						length = segment_sphere_intersection_length(vc, va, radius, center, c_in, intersection[j]);
//						j++;
//					}
//					angle = signed_angle(vc, va, vb, e.next().comp().next().next().end());
////					k3d::log() << debug << i << " : check e3 " << angle << " " << length << std::endl;
//					p = angle*length;
//					tensor(0,0) += p*edge[0]*edge[0];
//					tensor(1,1) += p*edge[1]*edge[1];
//					tensor(2,2) += p*edge[2]*edge[2];
//
//					tensor(0,1) += p*edge[0]*edge[1];
//					tensor(0,2) += p*edge[0]*edge[2];
//					tensor(1,2) += p*edge[1]*edge[2];
//				} else if((c_in && !a_in) || (!c_in && a_in)) {
//					// j should not be 2 here, because that would imply that a_in && c_in is true
//					segment_sphere_intersection_length(vc, va, radius, center, c_in, intersection[j]);
//					j++;
//				}
//
//				if(j == 2) {
//					double temp_area = triangle_area(va,vb,vc);
//
//					if(!a_in && b_in && c_in) 
//						temp_area -= std::abs(triangle_area(va, intersection[0], intersection[1]));
//					if(a_in && !b_in && c_in) 
//						temp_area -= std::abs(triangle_area(vb, intersection[0], intersection[1]));
//					if(a_in && b_in && !c_in) 
//						temp_area -= std::abs(triangle_area(vb, intersection[0], intersection[1]));
//
//					if(a_in && !b_in && !c_in) 
//						temp_area = std::abs(triangle_area(va, intersection[0], intersection[1]));
//					if(!a_in && b_in && !c_in) 
//						temp_area = std::abs(triangle_area(vb, intersection[0], intersection[1]));
//					if(!a_in && !b_in && c_in) 
//						temp_area = std::abs(triangle_area(vb, intersection[0], intersection[1]));
//
//					area += temp_area;
//				} else {
//					area += triangle_area(va,vb,vc);
//				}
//			}
//		}
//
//		tensor(1,0) = tensor(0,1);
//		tensor(2,0) = tensor(0,2);
//		tensor(2,1) = tensor(1,2);
//		
////		std::cout << std::endl << area << std::endl;
//		gmm::scale(tensor, 1.0/area);
//		std::vector<double> eigval;
//		eigval.resize(3);
//		gmm::dense_matrix<double> eigvect(3,3);
////		std::cout << tensor << std::endl;
//		gmm::symmetric_qr_algorithm(tensor, eigval, eigvect);
//
////		std::cout << eigvect << std::endl;
////		std::cout << eigval << std::endl;
//
//		int max,min,n;
//		double ev[3];
//		ev[0] = std::abs(eigval[0]);
//		ev[1] = std::abs(eigval[1]);
//		ev[2] = std::abs(eigval[2]);
//
//		if(ev[0] > ev[1] && ev[1] > ev[2]) {
//			max = 0; min = 1; n = 2;
//		}
//		if(ev[2] > ev[1] && ev[1] > ev[0]) {
//			max = 2; min = 1; n = 0;
//		}
//
//		if(ev[1] > ev[2] && ev[2] > ev[0]) {
//			max = 1; min = 2; n = 0;
//		}
//		if(ev[0] > ev[2] && ev[2] > ev[1]) {
//			max = 0; min = 2; n = 1;
//		}
//
//		if(ev[2] > ev[0] && ev[0] > ev[1]) {
//			max = 2; min = 0; n = 1;
//		}
//		if(ev[1] > ev[0] && ev[0] > ev[2]) {
//			max = 1; min = 0; n = 2;
//		}
//		//std::cout << max << " " << min << " " << n << std::endl << "-----------" << std::endl;
//
//		curv_dir0[0] = eigvect(0,max);
//		curv_dir0[1] = eigvect(1,max);
//		curv_dir0[2] = eigvect(2,max);
//
//		curv_dir1[0] = eigvect(0,min);
//		curv_dir1[1] = eigvect(1,min);
//		curv_dir1[2] = eigvect(2,min);
//
//		norm[0] = eigvect(0,n);
//		norm[1] = eigvect(1,n);
//		norm[2] = eigvect(2,n);
//
//		// cleanup
//		for(int i = 0; i < search.size(); i++) {
//			e = mesh->getEdge(search[i]);
//			face_searched[e.face()()] = false;
//		}
	}

	double diff_geom::principal_curve_tensor(vert_t vert, std::pair<double,double>& iso)
	{
		// Choose basis for plane
		// For each edge adjacent to vert
		//		project edge into plane basis
		//		get mean curv for edge
		//		accumulate a and b coefficients
		//	Solve linear system to get a and b
		//  c = 2 * Kh - a
		//  check error = kg - (a*c - b*b)
		
		double d_x, d_y, w, k, len_square;

		//double AtA[3] = {0,0,0}; // Symmetric 2x2
		//double Atb[2] = {0,0}; 

		double a,b,c;
		bool print = false;
		mesh_info::Vert v = mesh->getVert(vert);
		mesh_info::Edge e = v.edge();
		vert_t first = e();

		double Kh2 = mean_curv[vert].Length();
		if(std::abs(Kh2) < 1E-8) {
			iso.first = 0;
			iso.second = 0;
			return 0;
		}

		Vec3 ihat = tangent_basis_i[vert];
		Vec3 jhat = tangent_basis_j[vert];
		Vec3 khat = tangent_basis_k[vert];

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
		int valence = 0;

		do {
			area += area_mixed(e());
			valence++;
			e = e.comp().next();
		} while(e() != first);

		gmm::dense_matrix<double> A(valence, 3);
		gmm::dense_matrix<double> AtA(3, 3);		
		std::vector<double> B(valence, 0.0);
		std::vector<double> AtB(3,0.0);
		std::vector<double> X(3,0.0);

		int i = 0;
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

//			w = area*(cotangent(e()) + cotangent(e.comp()())) * len_square;
			w = mean_coords[e()];
			k = 2.0*(e.dir()*khat)/len_square;
			
			if(print)
				k3d::log()  << ", w = ("<<w<<")";

			//double x2y2 = (d_x*d_x - d_y*d_y);
			//double xy2  = 2.0*d_x*d_y;
			
			if(print)
				k3d::log() << ", k = (" << k << ")\n";
			w = std::sqrt(std::abs(w));
			A(i, 0) = w* d_x * d_x;
			A(i, 1) = w* 2.0 *d_x * d_y;
			A(i, 2) = w* d_y * d_y;
			B[i] = w * k;

			//AtA[0] += w * x2y2 * x2y2;
			//AtA[1] += w * xy2 * x2y2;
			//AtA[2] += w * xy2 * xy2;

			//double sum = w * (k - Kh2 * d_y * d_y);
			//Atb[0] += sum * x2y2;
			//Atb[1] += sum * xy2;

			e = e.comp().next();
			i++;
		} while(e() != first);

		//if(print) {
		//	k3d::log() << "AtA = (" << AtA[0] << ", " << AtA[1] << ", " << AtA[2] << ")" << std::endl; 
		//	k3d::log() << "Atb = (" << Atb[0] << ", " << Atb[1] << ")" << std::endl; 
		//	k3d::log() << "K G = (" << gaus_curv[vert] << ")" << std::endl; 			
		//}

		//double det = AtA[0]*AtA[2] - AtA[1]*AtA[1];
		std::vector<int> pivot; pivot.resize(3);

		
		gmm::mult(gmm::transposed(A), B, AtB);
		gmm::mult(gmm::transposed(A), A, AtA);

		gmm::lu_factor(AtA, pivot);
		double det = gmm::lu_det(AtA);

		if(std::abs(det) < 1E-10) {
			iso.first = 0;
			iso.second = 0;

			return 0;
		}
		gmm::lu_solve(AtA, pivot, X, AtB);
		
		a = X[0];
		b = X[1];
		c = X[2];

		//det = 1/det;

		//a = det*(AtA[2]*Atb[0] - AtA[1]*Atb[1]);
		//b = det*(AtA[0]*Atb[1] - AtA[1]*Atb[0]);
		//c = Kh2 - a;
		
		double error = gaus_curv[vert] - a*c + b*b;
		//if(error > 1) {
		//	b = sqrt(std::abs(a*c - gaus_curv[vert]));
		//}

		Vec3 tens = Vec3(a,b,c);
		iso = isotropic_tensor(tens);
		if(print) {
			k3d::log() << "Error = (" << error << ")" << std::endl; 
			k3d::log() << "abc = (" << a << ", " << b << ", " << c << ")" << std::endl;  
			k3d::log() << "iso = (" << iso.first << ", " << iso.second << ")" << std::endl;  
		}


		//double e1, e2;
		//eigen(tens, e1, e2);
		
		return error;	
	}
	
	std::pair<double,double> diff_geom::isotropic_tensor(Vec3 tens) 
	{
		double lambda = 0.5*(tens[0] + tens[2]);
		Vec3 iso = Vec3(tens[0] - lambda, tens[1], 0);

		iso.Normalize();
		std::pair<double,double> ret;
		ret.first = iso[0];
		ret.second = iso[1];

		return ret;
	}

	void diff_geom::eigen(Vec3 tens, double& e1, double& e2) 
	{
		double a = tens[0];
		double b = tens[1];
		double c = tens[2];
		double root = sqrt(a*a + 4*b*b - 2*a*c + c*c );

		e1 = 0.5*(a + c - root);
		e2 = 0.5*(a + c + root);
	}

	double diff_geom::gaussian_curvature(vert_t vert) 
	{
		edge_t edge = mesh->vert_edge[vert];

		double curv = k3d::pi_times_2();
		double area = 0.0;
		do {
			area += area_mixed(edge);
			Vec3 a(mesh->mesh->points->at(mesh->edge_vert[edge]).n);
			Vec3 b(mesh->mesh->points->at(mesh->edge_vert[mesh->edge_ccw[edge]]).n);
			Vec3 c(mesh->mesh->points->at(mesh->edge_vert[mesh->edge_ccw[mesh->edge_ccw[edge]]]).n);
			
			Vec3 AB = b-a;
			Vec3 AC = c-a;
		
			AB.Normalize();
			AC.Normalize();
			double dot = AB*AC;

			if( dot <= -1.0 ) return k3d::pi();  
			if( dot >=  1.0 ) return 0;  
			double angle = 
	
			curv -= acos(dot);
			edge = mesh->edge_ccw[mesh->edge_comp[edge]];

		} while(edge != mesh->vert_edge[vert]);

		return curv/area;
	}

	Vec3 diff_geom::mean_curvature(vert_t vert) 
	{
		edge_t edge = mesh->vert_edge[vert];
		double area = 0.0;
		Vec3 mean(0,0,0);

		do {
			area += area_mixed(edge);
			
			double w = edge_cot[edge] + edge_cot[mesh->edge_comp[edge]];
			Vec3 a(mesh->mesh->points->at(mesh->edge_vert[edge]).n);
			Vec3 b(mesh->mesh->points->at(mesh->edge_vert[mesh->edge_ccw[edge]]).n);
			
			a -= b;
			a *= w;

			mean += a;

			edge = mesh->edge_ccw[mesh->edge_comp[edge]];
		} while(edge != mesh->vert_edge[vert]);
		
		mean /= (2.0*area);
		return mean;
	}

	/// Voronoi region of vertex on edge intersecting with a triangle
	double diff_geom::area_mixed(edge_t edge) 
	{
		Vec3 a(mesh->mesh->points->at(mesh->edge_vert[edge]).n);
		Vec3 b(mesh->mesh->points->at(mesh->edge_vert[mesh->edge_ccw[edge]]).n);
		Vec3 c(mesh->mesh->points->at(mesh->edge_vert[mesh->edge_ccw[mesh->edge_ccw[edge]]]).n);

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
			area += AB_norm2 * edge_cot[mesh->edge_ccw[mesh->edge_ccw[edge]]];
			area *= 0.125;
		}

		return area;
	}

	double diff_geom::mean_weight(edge_t edge) 
	{
		mesh_info::Edge e(mesh->getEdge(edge));

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
		mesh_info::Edge e = mesh->getEdge(edge);
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
		edge_t next = mesh->edge_ccw[edge];

		// get the 3d positions of the vertices
		Vec3 a(mesh->mesh->points->at(mesh->edge_vert[edge]).n);
		Vec3 b(mesh->mesh->points->at(mesh->edge_vert[mesh->edge_ccw[edge]]).n);
		Vec3 c(mesh->mesh->points->at(mesh->edge_vert[mesh->edge_ccw[mesh->edge_ccw[edge]]]).n);
		
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
