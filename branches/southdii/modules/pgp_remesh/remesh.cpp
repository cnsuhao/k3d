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
	/// TODO
	void calc_edge_face_adj(const k3d::mesh::polyhedra_t& Polyhedra, indices_t& adj) 
	{
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
	void calc_edge_ccw_adj(const k3d::mesh::polyhedra_t& Polyhedra, indices_t& adj) 
	{
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
	void calc_edge_companion_adj(const k3d::mesh::polyhedra_t& Polyhedra, indices_t& adj) 
	{
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
	void calc_vert_edge_adj(const k3d::mesh::polyhedra_t& Polyhedra, indices_t& adj) 
	{
		for(size_t i = 0; i < Polyhedra.edge_points->size(); ++i) {
			size_t vert = Polyhedra.edge_points->at(i);
			adj.at(vert) = i;
		}
	}

	/// TODO: Add desc	
	void calc_vert_edge_ccw_adj(const k3d::mesh::polyhedra_t& Polyhedra, indices_t& adj) 
	{
		indices_t comp;
		comp.resize(Polyhedra.edge_points->size());
		calc_edge_companion_adj(Polyhedra, comp);
		for(size_t i = 0; i < Polyhedra.edge_points->size(); ++i) {
			size_t vert = Polyhedra.edge_points->at(comp.at(i));
			adj.at(vert) = i;
		}
	}

	/// TODO: Add desc	
	void calc_edge_vert_ccw_adj(const k3d::mesh::polyhedra_t& Polyhedra, indices_t& adj) 
	{
		adj.resize(Polyhedra.edge_points->size());
		for(size_t i = 0; i < Polyhedra.edge_points->size(); ++i) {
			adj[i] = Polyhedra.edge_points->at(Polyhedra.clockwise_edges->at(i));
		}
	}

	/// TODO: Add desc
	void calc_face_poly_adj(const k3d::mesh::polyhedra_t& Polyhedra, indices_t& adj) 
	{
		adj.resize(Polyhedra.face_first_loops->size());
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
			num_edges = Mesh.polyhedra->edge_points->size();
			num_faces = Mesh.polyhedra->face_first_loops->size();
			num_verts = Mesh.points->size();
			edge_comp.resize(num_edges);
			edge_ccw.resize(num_edges);
			edge_face.resize(num_edges);

			face_poly.resize(num_faces);
			vert_edge.resize(num_verts);
			edge_vert.resize(num_verts);

			mean_curv.resize(num_verts);
			curv_tens.resize(num_verts);
			gaus_curv.resize(num_verts);
			curv_tens.resize(num_verts);

			edge_cot.resize(num_edges);
			mean_curv.resize(num_verts);
			gaus_curv.resize(num_verts);

			calc_edge_face_adj(*mesh.polyhedra, edge_face);
			calc_edge_ccw_adj(*mesh.polyhedra, edge_ccw);
			calc_edge_companion_adj(*mesh.polyhedra, edge_comp);
			calc_vert_edge_ccw_adj(*mesh.polyhedra, vert_edge);
			calc_edge_vert_ccw_adj(*mesh.polyhedra, edge_vert);
			calc_face_poly_adj(*mesh.polyhedra, face_poly);
		}
		
		class Face;
		class Edge;
		class Vert;

		class Vert 
		{
		public:
			Vert(const mesh_info& _info, vert_t _vert = -1) :
				info(_info), index(_vert) {}
			
			Vert& operator=(vert_t vert) 
			{ 
				index = vert; 
				return *this; 
			}

			vert_t operator() () const
			{
				return index;
			}

			Vec3 pos() const
			{
				return Vec3(info.mesh.points->at(index).n);
			}

			Edge edge() const {
				return Edge(info, info.vert_edge[index]);
			}

			const mesh_info& info;
			vert_t index;
		};

		class Face 
		{
		public:
			Face(const mesh_info& _info, face_t _face = -1) :
				info(_info), index(_face) {}
			
			Face& operator=(face_t face) 
			{ 
				index = face; 
				return *this; 
			}

			Edge operator[] (size_t e) const
			{
				edge_t edge = info.mesh.polyhedra->loop_first_edges->at(
					info.mesh.polyhedra->face_first_loops->at(index));
				
				for(int i = 0; i < e; i++) {
					edge = info.edge_ccw[edge]; 
				}

				return Edge(info, edge);
			}
			
			face_t operator() () const
			{
				return index;
			}

			const mesh_info& info;
			face_t index;
		};

		class Edge 
		{
		public:
			Edge(const mesh_info& _info, edge_t _edge = -1) :
				info(_info), index(_edge) {}
			
			Edge(const Edge& e) :
				info(e.info), index(e.index) {}

			Edge& operator=(edge_t edge) 
			{ 
				index = edge; 
				return *this; 
			}
			
			Edge& operator=(const Edge& edge) 
			{ 
				index = edge.index; 
				return *this; 
			}

			// prefix counterclockwise traversal 
			Edge& operator++ ()
			{
				index = info.edge_ccw[index];
				return *this;
			}
			
			// postfix counterclockwise traversal
			Edge  operator++ (int) 
			{
				edge_t old = index;
				index = info.edge_ccw[index];
				return Edge(info, old);
			}
	
			// prefix clockwise traversal
			Edge& operator-- ()     
			{
				index = info.mesh.polyhedra->clockwise_edges->at(index);
				return *this;
			}
			
			// postfix clockwise traversal 
			Edge  operator-- (int) 
			{
				edge_t old = index;
				index = info.mesh.polyhedra->clockwise_edges->at(index);
				return Edge(info, old);
			}
	
			Edge next() const 
			{
				return Edge(info, info.edge_ccw[index]);
			}

			Edge prev() const 
			{
				return Edge(info, info.mesh.polyhedra->clockwise_edges->at(index));
			}

			Face face() const 
			{
				return Face(info, info.edge_face[index]);
			}

			Vert vert() const 
			{
				return Vert(info, info.edge_vert[index]);
			}

			Edge comp() const
			{
				return Edge(info, info.edge_comp[index]);
			}

			edge_t operator() () const
			{
				return index;
			}

			Vec3 start() const
			{
				return  Vec3(info.mesh.points->at(info.edge_vert[index]).n);
			}

			Vec3 end() const
			{
				return  Vec3(info.mesh.points->at(info.edge_vert[info.edge_ccw[index]]).n);
			}	

			Vec3 dir() const
			{
				return end()-start();
			}	

			const mesh_info& info;
			edge_t index;
		};




		void fill_diff_geom(k3d::mesh& OutputMesh) 
		{
			k3d::typed_array<k3d::vector3>* curv_p = new k3d::typed_array<k3d::vector3>;
			boost::shared_ptr<k3d::typed_array<k3d::vector3> > curv(curv_p);

			k3d::typed_array<k3d::vector3>* p1_p = new k3d::typed_array<k3d::vector3>;
			boost::shared_ptr<k3d::typed_array<k3d::vector3> > p1(p1_p);

			k3d::typed_array<k3d::vector3>* p2_p = new k3d::typed_array<k3d::vector3>;
			boost::shared_ptr<k3d::typed_array<k3d::vector3> > p2(p2_p);
			curv->resize(num_verts);
			p1->resize(num_verts);
			p2->resize(num_verts);

			k3d::log() << debug << "Start fill" << std::endl;
			for(size_t i = 0; i < num_edges; ++i) {
				edge_cot[i] = cotangent(i);
			}
			k3d::log() << debug << "Done cot" << std::endl;
			for(size_t i = 0; i < num_verts; ++i) {
				mean_curv[i] = mean_curvature(i);
				curv->at(i).n[0] = mean_curv[i][0];
				curv->at(i).n[1] = mean_curv[i][1];
				curv->at(i).n[2] = mean_curv[i][2];
			}
			k3d::log() << debug << "Done mean" << std::endl;
			for(size_t i = 0; i < num_verts; ++i) {
				gaus_curv[i] = gaussian_curvature(i);
			}
			Vec3 pc1,pc2;
			for(size_t i = 0; i < num_verts; ++i) {
				principal_curve_tensor(i, pc1,  pc2);
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



			k3d::log() << debug << "Done gaussian" << std::endl;
		}

		Vec3 normal(vert_t vert) {
			return mean_curv[vert];	
		}

		double principal_curve_tensor(vert_t vert, Vec3& curv_dir0,  Vec3& curv_dir1)
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
			double m[4] = {0,0,0,0};
			double x[2] = {0,0};
			double a,b,c;

			Vert v(*this, vert);
			Edge e = v.edge();
			edge_t first = e();
			

			ihat = e.dir();
			ihat.Normalize();
			
			khat = normal(vert);
			jhat = khat ^ ihat;
			ihat = jhat ^ khat;

			ihat.Normalize();
			jhat.Normalize();

			do {
				d_x = e.dir() * ihat;		
				d_y = e.dir() * jhat;
				len_square = e.dir() * e.dir(); 

				w = (edge_cot[e()] + edge_cot[e.comp()()]) * len_square;
				k = 2.0*(e.dir()*khat)/len_square;
				
				m[0] += w * d_x * d_x * d_x * d_x;
				m[1] += 2.0* w * d_x * d_x * d_x * d_y;
				m[2] += w * d_x * d_x * d_x * d_y;
				m[3] += w * d_x * d_x * d_y * d_y;

				x[0] += w * d_x * d_x * k;
				x[1] += w * d_x * d_y * k;

				e = e.comp().next();
			} while(e() != first);

			double det = m[0]*m[3] - m[1]*m[2];
			
			if(det == 0.0) {
				curv_dir0 = Vec3();
				curv_dir1 = Vec3();
				return 0;
			}
			det = 1/det;

			a = det*(m[3]*x[0] - m[1]*x[1]);
			b = det*(-1.0*m[2]*x[0] + m[0]*x[1]);
			c = 0.5*mean_curv[vert].Length() - a;
			
			double error = gaus_curv[vert] - a*c + b*b;
			
			//Vec3 vi = Vec3(mesh.points->at(vert).n);
			//Vec3 vj = Vec3(mesh.points->at(vert_edge[edge_comp[edge]]).n);
			Vec3 iso = isotropic_tensor(Vec3(a,b,c));
			
			double *t = iso;
			double angle = atan2(t[1], t[0])*0.5;

			Vec3 temp_i = ihat;
			Vec3 temp_j = jhat;

			temp_i *= cos(angle);
			temp_j *= sin(angle);
			curv_dir0 = temp_i + temp_j;

			temp_i = ihat;
			temp_j = jhat;
			temp_i *= cos(angle + k3d::pi_over_2());
			temp_j *= sin(angle + k3d::pi_over_2());
			curv_dir1 = temp_i + temp_j;
	
			return error;	
		}
		
		Vec3 isotropic_tensor(Vec3 tensor) 
		{
			double *t = tensor;
			
			double lambda = 0.5*(t[0] + t[2]);
			Vec3 iso = Vec3(t[0] - lambda, t[1], 0);

			iso.Normalize();

			return iso;
		}

		double gaussian_curvature(vert_t vert) 
		{
			edge_t edge = vert_edge[vert];

			double curv = k3d::pi_times_2();
			double area = 0.0;
			do {
				area += area_mixed(edge_face[edge], edge);
				Vec3 a(mesh.points->at(edge_vert[edge]).n);
				Vec3 b(mesh.points->at(edge_vert[edge_ccw[edge]]).n);
				Vec3 c(mesh.points->at(edge_vert[edge_ccw[edge_ccw[edge]]]).n);
				
				Vec3 AB = b-a;
				Vec3 AC = c-a;
			
				AB.Normalize();
				AC.Normalize();

				curv -= acos(AB * AC);
				edge = edge_ccw[edge_comp[edge]];

			} while(edge != vert_edge[vert]);

			return curv/area;
		}

		Vec3 mean_curvature(vert_t vert) 
		{
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

		/// Voronoi region of vertex on edge intersecting with a triangle
		double area_mixed(face_t face, edge_t edge) 
		{
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
			if(cosA < 0.0 || cosB <= 0.0 || cosC <= 0.0) {
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

		/// Cotangent of the angle of the vertex opposite of edge
		double cotangent(edge_t edge) 
		{
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

		std::vector<Vec3> face_i_basis;
		std::vector<Vec3> face_j_basis;

		std::vector<Vec3> mean_curv;
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
			OutputMesh = InputMesh;
			k3d::typed_array<k3d::vector3>* curv_p = new k3d::typed_array<k3d::vector3>;
			boost::shared_ptr<k3d::typed_array<k3d::vector3> > curv(curv_p);
			m.fill_diff_geom(OutputMesh);
			curv->resize(OutputMesh.points->size());

			// Will do this more efficiently later
			for(int i = 0; i < curv->size(); i++) {
				curv->at(i).n[0] = m.mean_curv[i][0];
				curv->at(i).n[1] = m.mean_curv[i][1];
				curv->at(i).n[2] = m.mean_curv[i][2];
			}

			OutputMesh.vertex_data["PGPMeanCurv"] = curv;

			k3d::log() << debug << "PGP: create mesh: " << curv.use_count() << " " << curv->size() << std::endl;
		}
		void on_update_mesh(const k3d::mesh& InputMesh, k3d::mesh& OutputMesh)		  
		{
			detail::mesh_info m(InputMesh); 
			OutputMesh = InputMesh;
			k3d::typed_array<k3d::vector3>* curv_p = new k3d::typed_array<k3d::vector3>;
			boost::shared_ptr<k3d::typed_array<k3d::vector3> > curv(curv_p);
			m.fill_diff_geom(OutputMesh);

			curv->resize(OutputMesh.points->size());

			// Will do this more efficiently later
			for(int i = 0; i < curv->size(); i++) {
				curv->at(i).n[0] = m.mean_curv[i][0];
				curv->at(i).n[1] = m.mean_curv[i][1];
				curv->at(i).n[2] = m.mean_curv[i][2];
			}

			OutputMesh.vertex_data["PGPMeanCurv"] = curv;

			k3d::log() << debug << "PGP: update mesh" << std::endl;
		}

		static k3d::iplugin_factory& get_factory()
		{
		static k3d::document_plugin_factory<pgp_remesh,
			k3d::interface_list<k3d::imesh_source,
			k3d::interface_list<k3d::imesh_sink > > > factory(
			  k3d::uuid(0xc97aa4ce, 0x412c1ed2, 0x055044a4, 0xa151f085),
			  "PGP Remesh",
			  _("Quad remeshing using the PGP algorithm"),
			  "PGP",
			  k3d::iplugin_factory::EXPERIMENTAL);
			
			return factory;

		}
	};


	k3d::iplugin_factory& pgp_remesh_factory()
	{
		return pgp_remesh::get_factory();
	}
} // namespace pgp_module
