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
#include "pgp.h"
#include "mesh_info.h"
#include "diff_geom.h"
#include <gmm/gmm.h>

namespace libk3dquadremesh
{

namespace detail {
	void PGP::setup(double omega) 
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
			
			pf.lamda[0] = 1;
			pf.lamda[1] = 1;
			pf.lamda[2] = 1;

			//std::cout << "t = " << f << " L = " << "(" << x[0] << ", " << x[1] << ", " << x[2] << ")\n";

			double max1 = pf.vf[0]*pf.vf[1]; 
			double max2 = pf.vf[0]*pf.vf[2]; 

			int max_arg1 = 0;
			int max_arg2 = 0;

			//for(int i = 2; i < 4; i+=2) {
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

			pf.rot[0] = 0;
			pf.rot[1] = max_arg1;
			pf.rot[2] = max_arg2;


			for(int i = 0; i < 3; ++i) {
				int v0 = (i+1)%3;
				int v1 = (i+2)%3;

				pf.delta[i] = omega*0.5*((rotate90(pf.vf[v0], pf.rot[v0]) + rotate90(pf.vf[v1], pf.rot[v1]))*pf.e[i]);
				pf.delta_p[i] = omega*0.5*((rotate90(pf.vf[v0], pf.rot[v0]+1) + rotate90(pf.vf[v1], pf.rot[v1]+1))*pf.e[i]);
			}
		}
	}

	void PGP::solve() 
	{
		gmm::dense_matrix<double> D(4,4);
		gmm::dense_matrix<double> Temp1(4,4);
		gmm::dense_matrix<double> Temp2(4,4);
		M[0].resize(4,4);
		M[1].resize(4,4);
		M[2].resize(4,4);

		//M[0](0,2) = -1;
		//M[0](1,3) = 1;

		M[0](0,2) = 1;
		M[0](1,3) = -1;

		M[0](2,0) = 1;
		M[0](3,1) = 1;

		k3d::log() << -1 << ": mult a" <<std::endl;
		gmm::mult(M[0], M[0], M[1]);
		k3d::log() << -1 << ": mult b" <<std::endl;
		gmm::mult(M[0], M[1], M[2]);
		
		std::cout << D << std::endl;			
		std::cout << M[0] << std::endl;			
		std::cout << M[1] << std::endl;			
		std::cout << M[2] << std::endl;			
		
		std::vector<int> mapping(mesh->num_verts, 0);
		

		for(size_t p = 0; p < mesh->mesh->polyhedra->first_faces->size(); p++) {
			face_t f = mesh->mesh->polyhedra->first_faces->at(p);
			mapping[face_data[f].vert[0]] = -face_data[f].vert[0];
			k3d::log() << "poly constrain " << face_data[f].vert[0] <<std::endl;
		}
		
		int curr = 0;
		for(vert_t v = 0; v < mapping.size(); ++v) {
			if(mapping[v] == 0) {
				mapping[v] = curr++;
			}
		}
		gmm::row_matrix< gmm::wsvector<double> > A(4*curr, 4*curr);
		std::vector<double> X(4*curr, 0.0);
		std::vector<double> B(4*curr, 0.0);

		std::vector<double> sum_lamba(mesh->num_verts, 0.0);
		
		for(face_t f = 0; f < mesh->num_faces; ++f) {
			for(int e = 0; e < 3; ++e) {
				vert_t i = face_data[f].vert[(e+1)%3];
				vert_t j = face_data[f].vert[(e+2)%3];
				sum_lamba[i] += face_data[f].lamda[e];
				sum_lamba[j] += face_data[f].lamda[e];
			}								
		}
		
		for(vert_t v = 0; v < mapping.size(); ++v) {
			if(mapping[v] >= 0) {
				vert_t vm = 4*mapping[v];
				//std::cout << v << "->" << vm << " = " << sum_lamba[v] << std::endl;
				A(vm + 0, vm + 0) = sum_lamba[v];
				A(vm + 1, vm + 1) = sum_lamba[v];
				A(vm + 2, vm + 2) = sum_lamba[v];
				A(vm + 3, vm + 3) = sum_lamba[v];
			}
		}
					
		for(face_t f = 0; f < mesh->num_faces; ++f) {
			for(int e = 0; e < 3; ++e) {
				int v0 = (e+1)%3;
				int v1 = (e+2)%3;
				int i = 4*mapping[face_data[f].vert[v0]];
				int j = 4*mapping[face_data[f].vert[v1]];
				bool constrain0 = (mapping[face_data[f].vert[v0]] < 0);
				bool constrain1 = (mapping[face_data[f].vert[v1]] < 0);

				// Taken care of in seperate loop
				//if(!constrain0) {
				//	A(i + 0, i + 0) += face_data[f].lamda[e];
				//	A(i + 1, i + 1) += face_data[f].lamda[e];
				//	A(i + 2, i + 2) += face_data[f].lamda[e];
				//	A(i + 3, i + 3) += face_data[f].lamda[e];
				//}

				//if(!constrain1) {
				//	A(j + 0, j + 0) += face_data[f].lamda[e];
				//	A(j + 1, j + 1) += face_data[f].lamda[e];
				//	A(j + 2, j + 2) += face_data[f].lamda[e];
				//	A(j + 3, j + 3) += face_data[f].lamda[e];
				//}

				double cos_d = face_data[f].lamda[e]*std::cos(face_data[f].delta[e]);
				double sin_d = face_data[f].lamda[e]*std::sin(face_data[f].delta[e]);
				double cos_dt = face_data[f].lamda[e]*std::cos(face_data[f].delta_p[e]);
				double sin_dt = face_data[f].lamda[e]*std::sin(face_data[f].delta_p[e]);
				if(f == 0) {

//					std::cout << cos_d << " " << cos_dt << " " << sin_d << " "<< sin_dt << " " << std::endl;
//					std::cout << D << std::endl;
				}
				gmm::clear(D);
//				if(f == 0) std::cout << "0 " << D << std::endl;
				D(0,0) = -cos_d;
//				if(f == 0) std::cout << "1 " << D << std::endl;
				D(1,1) = -cos_d;
//				if(f == 0) std::cout << "2 " << D << std::endl;
				D(0,1) = sin_d;
//				if(f == 0) std::cout << "3 " << D << std::endl;
				D(1,0) = -sin_d;
//				if(f == 0) std::cout << "4 " << D << std::endl;

				D(2,2) = -cos_dt;
//				if(f == 0) std::cout << "5 " << D << std::endl;
				D(3,3) = -cos_dt;
//				if(f == 0) std::cout << "6 " << D << std::endl;
				D(2,3) = sin_dt;
//				if(f == 0) std::cout << "7 " << D << std::endl;
				D(3,2) = -sin_dt;
//				if(f == 0) std::cout << "8 " << D << std::endl;
//					std::cout << D << std::endl;
				if(f == 0) {
//					std::cout << "8 " << D << std::endl;
				}
				
				int r0 = face_data[f].rot[v0];
				int r1 = face_data[f].rot[v1];
				if((r1 % 2) == 1) r1 = (r1+2)%4;

				//k3d::log() << f << ": " << r0 << " " << r1 << std::endl;
				if(r0 == 0 && r1 == 0) {
					gmm::copy(D, Temp2);
				} else if (r0 == 0 && r1 != 0) {
					//k3d::log() << "mult 1" <<std::endl;
					gmm::mult(M[r1-1], D, Temp2);
				} else if (r0 != 0 && r1 == 0) {
					//k3d::log() << "mult 2" <<std::endl;
					gmm::mult(D, M[r0-1], Temp2);
				} else {
					//k3d::log() << "mult 3" <<std::endl;
					gmm::mult(D, M[r0-1], Temp1);
					//k3d::log() << "mult 4" <<std::endl;
					gmm::mult(M[r1-1], Temp1, Temp2);
				}
				if(f == 0) {
					std::cout << Temp2 << std::endl;
				}

				if(constrain0 && !constrain1) {
					//std::cout << "con0 |";
					double den = 1.0/sum_lamba[face_data[f].vert[v0]];
					//std::cout << "|";

					B[j + 0] -= (Temp2(0,0) + Temp2(0,2))*den;
					B[j + 1] -= (Temp2(1,0) + Temp2(1,2))*den;
					B[j + 2] -= (Temp2(2,0) + Temp2(2,2))*den;
					B[j + 3] -= (Temp2(3,0) + Temp2(3,2))*den;
					//std::cout << "|" <<std::endl;
					
				} else if(!constrain0 && constrain1) {
					//std::cout << "con1 |";
					double den = 1.0/sum_lamba[face_data[f].vert[v1]];
					//std::cout << "|";

					B[i + 0] -= (Temp2(0,0) + Temp2(2,0))*den;
					B[i + 1] -= (Temp2(0,1) + Temp2(2,1))*den;
					B[i + 2] -= (Temp2(0,2) + Temp2(2,2))*den;
					B[i + 3] -= (Temp2(0,3) + Temp2(2,3))*den;
					//std::cout << "|" <<std::endl;

				} else if(!constrain0 && !constrain1) {
					//std::cout << "- | " << i << " " << j << std::endl;
					A(j + 0, i + 0) += Temp2(0, 0);
					A(j + 0, i + 1) += Temp2(0, 1);
					A(j + 0, i + 2) += Temp2(0, 2);
					A(j + 0, i + 3) += Temp2(0, 3);

					A(j + 1, i + 0) += Temp2(1, 0);
					A(j + 1, i + 1) += Temp2(1, 1);
					A(j + 1, i + 2) += Temp2(1, 2);
					A(j + 1, i + 3) += Temp2(1, 3);

					A(j + 2, i + 0) += Temp2(2, 0);
					A(j + 2, i + 1) += Temp2(2, 1);
					A(j + 2, i + 2) += Temp2(2, 2);
					A(j + 2, i + 3) += Temp2(2, 3);

					A(j + 3, i + 0) += Temp2(3, 0);
					A(j + 3, i + 1) += Temp2(3, 1);
					A(j + 3, i + 2) += Temp2(3, 2);
					A(j + 3, i + 3) += Temp2(3, 3);

					///// 

					A(i + 0, j + 0) += Temp2(0, 0);
					A(i + 0, j + 1) += Temp2(1, 0);
					A(i + 0, j + 2) += Temp2(2, 0);
					A(i + 0, j + 3) += Temp2(3, 0);

					A(i + 1, j + 0) += Temp2(0, 1);
					A(i + 1, j + 1) += Temp2(1, 1);
					A(i + 1, j + 2) += Temp2(2, 1);
					A(i + 1, j + 3) += Temp2(3, 1);

					A(i + 2, j + 0) += Temp2(0, 2);
					A(i + 2, j + 1) += Temp2(1, 2);
					A(i + 2, j + 2) += Temp2(2, 2);
					A(i + 2, j + 3) += Temp2(3, 2);

					A(i + 3, j + 0) += Temp2(0, 3);
					A(i + 3, j + 1) += Temp2(1, 3);
					A(i + 3, j + 2) += Temp2(2, 3);
					A(i + 3, j + 3) += Temp2(3, 3);
					//std::cout << "-" << std::endl;
				}
			}
		}

		gmm::clean(A, 1E-10);
		gmm::csr_matrix<double> S(4*curr, 4*curr);
		gmm::copy(A, S);
		
/*			gmm::dense_matrix<double> print_out(4*curr, 4*curr);
		gmm::copy(A, print_out);
		
		std::ofstream fout("bblah.mat");
		for(vert_t v = 0; v < 4*curr; ++v) {
			for(vert_t w = 0; w < 4*curr; ++w) {
				fout << print_out(v,w) << ' ';
			}
			fout << std::endl;
		}
		fout.close();
*/		
		gmm::diagonal_precond<gmm::csr_matrix<double> > PR(S);
		//gmm::identity_matrix PR;
		gmm::iteration iter(1E-10);
		iter.set_noisy(0);
		gmm::cg(S, X, B, PR, iter);
		
		for(vert_t v = 0; v < mapping.size(); ++v) {
			if(mapping[v] < 0) {
				vert_data[v].theta = 0;
				vert_data[v].phi   = 0;
			} else {
				vert_t vm = 4*mapping[v];
				vert_data[v].theta = atan2(X[vm + 1], X[vm + 0]);
				vert_data[v].phi   = atan2(X[vm + 3], X[vm + 2]);
				if(vert_data[v].theta < 0) vert_data[v].theta += k3d::pi_times_2();
				if(vert_data[v].phi < 0) vert_data[v].phi += k3d::pi_times_2();
			}

			//std::cout << v << ": (" << vert_data[v].theta << ", " << vert_data[v].phi << ")" << std::endl;
			
		}
		

	}

	void PGP::extract(double omega) 
	{
		for(face_t f = 0; f < mesh->num_faces; ++f) {
			per_face &pf = face_data[f];

			pf.theta[0] = vert_data[pf.vert[0]].theta;
			pf.theta[1] = vert_data[pf.vert[1]].theta;
			pf.theta[2] = vert_data[pf.vert[2]].theta;

			pf.phi[0] = vert_data[pf.vert[0]].phi;
			pf.phi[1] = vert_data[pf.vert[1]].phi;
			pf.phi[2] = vert_data[pf.vert[2]].phi;

			for(int v = 0; v < 2; v++) {
				// propagate from v to (v+1)%3 over edge (v+2)%3
				int i = v%3;
				int j = (v+1)%3;
				int e = (v+2)%3;

				int r,s,t;

				vec2 n = pf.e[e];
				n = n*(1.0/sqrt(n*n));

				double max_dot = pf.vf[i]*pf.vf[j]; 
				r = 0;

				//for(int i = 2; i < 4; i+=2) {
				for(int k = 1; k < 4; k++) {
					double dot = pf.vf[i]*rotate90(pf.vf[j], k);

					if(dot > max_dot) {
						max_dot = dot;
						r = k;
					}
				}

				pf.vf[j] = rotate90(pf.vf[j], r);
				vec2 param = vec2(pf.theta[j], pf.phi[j]);
				param = rotate90(param, r);
				pf.theta[j] = param.first;
				pf.phi[j] = param.second;

				s = (int)floor(-1*(pf.theta[i] - (k3d::pi()/omega)*(n*(pf.vf[i] + pf.vf[j])) - pf.theta[j])/k3d::pi_times_2() + 0.5);
				t = (int)floor(-1*(pf.phi[i] - (k3d::pi()/omega)*(n*(rotate90(pf.vf[i], 1) + rotate90(pf.vf[j], 1))) - pf.phi[j])/k3d::pi_times_2() + 0.5);
				pf.theta[j] -= k3d::pi_times_2()*s;
				pf.phi[j] -= k3d::pi_times_2()*t;
			}
		}
	}

	void PGP::remesh(k3d::mesh& OutputMesh) 
	{
		k3d::typed_array<k3d::color>* c1_p = new k3d::typed_array<k3d::color>;
		boost::shared_ptr<k3d::typed_array<k3d::color> > c1(c1_p);

		k3d::typed_array<k3d::color>* c2_p = new k3d::typed_array<k3d::color>;
		boost::shared_ptr<k3d::typed_array<k3d::color> > c2(c2_p);

		k3d::typed_array<vec2>* tex_p = new k3d::typed_array<vec2>;
		boost::shared_ptr<k3d::typed_array<vec2> > tex(tex_p);

		c1->resize(mesh->num_edges);
		c2->resize(mesh->num_edges);
		tex->resize(mesh->num_edges);

		k3d::mesh::polyhedra_t* poly = new k3d::mesh::polyhedra_t(*OutputMesh.polyhedra);
		boost::shared_ptr<k3d::mesh::polyhedra_t> poly1(poly);
		k3d::mesh::named_arrays& a1 = poly->face_varying_data;

		a1["PGP_pre_theta_color"] = c1;
		a1["PGP_pre_phi_color"] = c2;
		a1["PGP_uv"] = tex;
		double angle;
	
		for(face_t f = 0; f < mesh->num_faces; ++f) {
			edge_t e0 = face_data[f].edge[0];
			edge_t e1 = face_data[f].edge[1];
			edge_t e2 = face_data[f].edge[2];
			//c1->at(e0) = k3d::color(k3d::basic_hsv((180.0/k3d::pi()) * vert_data[mesh->edge_vert[e1]].theta,1,1));
			//c1->at(e1) = k3d::color(k3d::basic_hsv((180.0/k3d::pi()) * vert_data[mesh->edge_vert[e2]].theta,1,1));
			//c1->at(e2) = k3d::color(k3d::basic_hsv((180.0/k3d::pi()) * vert_data[mesh->edge_vert[e0]].theta,1,1));
			//c2->at(e0) = k3d::color(k3d::basic_hsv((180.0/k3d::pi()) * vert_data[mesh->edge_vert[e1]].phi,1,1));
			//c2->at(e1) = k3d::color(k3d::basic_hsv((180.0/k3d::pi()) * vert_data[mesh->edge_vert[e2]].phi,1,1));
			//c2->at(e2) = k3d::color(k3d::basic_hsv((180.0/k3d::pi()) * vert_data[mesh->edge_vert[e0]].phi,1,1));
			angle = face_data[f].theta[0];
			while(angle < 0) angle += k3d::pi_times_2();
			while(angle > k3d::pi_times_2()) angle -= k3d::pi_times_2();
			c1->at(e1) = k3d::color(k3d::basic_hsv((180.0/k3d::pi()) * angle,1,1));
			tex->at(e1) = vec2(face_data[f].theta[0]/k3d::pi_times_2(), face_data[f].phi[0]/k3d::pi_times_2());

			angle = face_data[f].theta[1];
			while(angle < 0) angle += k3d::pi_times_2();
			while(angle > k3d::pi_times_2()) angle -= k3d::pi_times_2();
			c1->at(e2) = k3d::color(k3d::basic_hsv((180.0/k3d::pi()) * angle,1,1));
			tex->at(e2) = vec2(face_data[f].theta[1]/k3d::pi_times_2(), face_data[f].phi[1]/k3d::pi_times_2());

			angle = face_data[f].theta[2];
			while(angle < 0) angle += k3d::pi_times_2();
			while(angle > k3d::pi_times_2()) angle -= k3d::pi_times_2();
			c1->at(e0) = k3d::color(k3d::basic_hsv((180.0/k3d::pi()) * angle,1,1));
			tex->at(e0) = vec2(face_data[f].theta[2]/k3d::pi_times_2(), face_data[f].phi[2]/k3d::pi_times_2());

			angle = face_data[f].phi[0];
			while(angle < 0) angle += k3d::pi_times_2();
			while(angle > k3d::pi_times_2()) angle -= k3d::pi_times_2();
			c2->at(e1) = k3d::color(k3d::basic_hsv((180.0/k3d::pi()) * angle,1,1));

			angle = face_data[f].phi[1];
			while(angle < 0) angle += k3d::pi_times_2();
			while(angle > k3d::pi_times_2()) angle -= k3d::pi_times_2();
			c2->at(e2) = k3d::color(k3d::basic_hsv((180.0/k3d::pi()) * angle,1,1));

			angle = face_data[f].phi[2];
			while(angle < 0) angle += k3d::pi_times_2();
			while(angle > k3d::pi_times_2()) angle -= k3d::pi_times_2();
			c2->at(e0) = k3d::color(k3d::basic_hsv((180.0/k3d::pi()) * angle,1,1));
		}

		OutputMesh.polyhedra = poly1;
	}

};

}; // namespace pgp_module
