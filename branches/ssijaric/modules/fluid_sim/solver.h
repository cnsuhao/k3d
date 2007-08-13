#ifndef FLUID_SIM_SOLVER_H
#define FLUID_SIM_SOLVER_H

#include <k3dsdk/document_plugin_factory.h>
#include <k3dsdk/imaterial.h> 
#include <k3dsdk/node.h>
#include <k3dsdk/persistent.h> 
#include <k3dsdk/material.h>
#include <k3dsdk/material_client.h>
#include <k3dsdk/measurement.h>
#include <k3dsdk/mesh_source.h>
#include <k3dsdk/mesh_modifier.h>
#include <k3dsdk/mesh.h>
#include <k3dsdk/module.h>
#include <k3dsdk/property.h>

#include <gmm/gmm.h>


#include "voxel_grid.h"
#include <list>

namespace fluid_sim
{

class solver : public k3d::node
{
public:
	solver(k3d::iplugin_factory& Factory, k3d::idocument& Document);
	static k3d::iplugin_factory& get_factory();

protected:
	typedef gmm::row_matrix< gmm::wsvector<float> > sparse_matrix;
	k3d_data(voxel_grid*, immutable_name, change_signal, no_undo, node_storage, no_constraint, node_property, no_serialization) m_voxel_grid;
	k3d_data(float, immutable_name, change_signal, with_undo, local_storage, no_constraint, read_only_property, no_serialization) m_timestep;
	k3d_data(int, immutable_name, change_signal, with_undo, local_storage, no_constraint, read_only_property, no_serialization) m_steps;
	k3d_data(float, immutable_name, change_signal, with_undo, local_storage, no_constraint, read_only_property, no_serialization) m_viscosity;


	int m_xfaces, m_yfaces, m_zfaces; // get from the voxel_grid
	int m_xvox, m_yvox, m_zvox;
	float m_vox_width;
	float m_dt;


	sigc::slot<void, iunknown*> start_solver_slot();

	void diffuse_velocities(voxel_grid& new_grid, const voxel_grid& old_grid);
	void project(voxel_grid& new_grid, const voxel_grid& old_grid);
	void add_force(voxel_grid& u, const array3d_i& forcesx, const array3d_i& forcesy, const array3d_i& forcesz);
	k3d::point3 trace_particle(const k3d::point3& p);


	struct idx {
		idx(int i,int j,int k) : m_i(i), m_j(j), m_k(k) { }
		int m_i, m_j, m_k;
	};

	void setup_diffusion_velocity_bc(array3d_i& faces);
	void setup_and_solve_diffusion(const std::list<idx>& faces, array3d_i& faces_r, voxel_grid& u, const voxel_grid& w, voxel_grid::velocity_type vtype);





};

k3d::iplugin_factory& solver_factory();

}

#endif
