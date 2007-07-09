#ifndef FLUID_SIM_VOXEL_GRID_H
#define FLUID_SIM_VOXEL_GRID_H

#include <k3dsdk/document_plugin_factory.h>
#include <k3dsdk/module.h>
#include <k3dsdk/node.h>
#include <k3dsdk/measurement.h>
#include <k3dsdk/persistent.h>
#include <k3dsdk/i18n.h>
#include <k3dsdk/mesh_source.h>

#include <boost/multi_array.hpp>

namespace fluid_sim
{

struct particle 
{
	k3d::point3 m_pos;
	k3d::vector3 m_u; // velocity
	double m_phi;
	double m_radius;
};

// each cell contains velocities on its sides - since 3 of the velocities are shared with the cell's neighbours, only 3 are 
// needed to be stored for each cell
class cell 
{
public:
	cell() : m_num_particles(0) { }
	double velocity_x() { return m_vx; }
	double velocity_y() { return m_vy; }
	double velocity_z() { return m_vz; }
private:
	double m_phi; // stored at the center
	double m_num_particles; // number of particles in the cell
	// velocities at faces of each cell
	double m_vx; 
	double m_vy;
	double m_vz;

	std::list<particle*> m_particles; // list of particles that this cell contains (will most likely be changed)
};

// the voxel grid is specifed in similar fashion as the k3d bounding box - with (nx, ny, nz) and (px, py, pz)
class voxel_grid : public k3d::node
{
public:
	voxel_grid(k3d::iplugin_factory& Factory, k3d::idocument& Document); 
	static k3d::iplugin_factory& get_factory();
	sigc::signal1<void, k3d::iunknown*> voxel_grid_changed_signal; 
	void vgrid_modified_x(k3d::iunknown* Hint);
	void vgrid_modified_y(k3d::iunknown* Hint);
	void vgrid_modified_z(k3d::iunknown* Hint);
	void vgrid_modified(k3d::iunknown* Hint);

	sigc::slot<void, k3d::iunknown*> vgrid_modified_slot_x();
	sigc::slot<void, k3d::iunknown*> vgrid_modified_slot_y();
	sigc::slot<void, k3d::iunknown*> vgrid_modified_slot_z();
	sigc::slot<void, k3d::iunknown*> vgrid_modified_slot();

	k3d::vector3 velocity(const k3d::point3& pos); // trilinearly interpolated velocity at pos
	const cell& get_cell(const k3d::point3& pos);

	
private:
	typedef boost::multi_array<cell, 3> array_type;


	// many of the properties will be removed - only voxel width will be modifiable by the user
	k3d_data(double, immutable_name, change_signal, with_undo, local_storage, no_constraint, writable_property, no_serialization) m_nx;
        k3d_data(double, immutable_name, change_signal, with_undo, local_storage, no_constraint, writable_property, no_serialization) m_ny;
	k3d_data(double, immutable_name, change_signal, with_undo, local_storage, no_constraint, writable_property, no_serialization) m_nz;

	k3d_data(double, immutable_name, change_signal, with_undo, local_storage, no_constraint, writable_property, no_serialization) m_px;
        k3d_data(double, immutable_name, change_signal, with_undo, local_storage, no_constraint, writable_property, no_serialization) m_py;
	k3d_data(double, immutable_name, change_signal, with_undo, local_storage, no_constraint, writable_property, no_serialization) m_pz;


	k3d_data(double, immutable_name, change_signal, with_undo, local_storage, no_constraint, writable_property, no_serialization) m_vox_width;
	
	k3d_data(int, immutable_name, change_signal, with_undo, local_storage, no_constraint, read_only_property, no_serialization) m_rows;
	k3d_data(int, immutable_name, change_signal, with_undo, local_storage, no_constraint, read_only_property, no_serialization) m_cols;
	k3d_data(int, immutable_name, change_signal, with_undo, local_storage, no_constraint, read_only_property, no_serialization) m_slices;

	double m_orig_nx;
	double m_orig_ny;
	double m_orig_nz;

	double m_orig_px;
	double m_orig_py;
	double m_orig_pz;

	k3d::vector3 m_norigin; // <nx, ny, nz>
	k3d::vector3 m_porigin; // <px, py, pz>
	
	// voxel grid
	array_type* m_grid;


	double lerp(double fx, double v0, double v1) { return v0 + fx*(v1-v0); }

	// (nx, ny, nz) of the voxel (i,j,k)
	k3d::vector3 n_coord(int i, int j, int k) { return m_norigin + k3d::vector3(i, j, k) * m_vox_width.value(); }

	const cell& get_voxel(int i, int j, int k) { return (*m_grid)[i][j][k]; };

	// boost's [i][j][k] performs range checking by default - can be disabled
	// convenience function for getting velocities stored in individual voxels
	// more efficient to just call get_voxel, and use cell's accessor methods to access velocities
	double get_velocity_i(int i, int j, int k) { return ((*m_grid)[i][j][k]).velocity_x(); }
	double get_velocity_j(int i, int j, int k) { return ((*m_grid)[i][j][k]).velocity_y(); }
	double get_velocity_k(int i, int j, int k) { return ((*m_grid)[i][j][k]).velocity_z(); }

	// return the location of the specified velocity component for voxel (i,j,k)
	k3d::vector3 get_velocity_i_pos(int i, int j, int k) { }
	k3d::vector3 get_velocity_j_pos(int i, int j, int k) { }
	k3d::vector3 get_velocity_k_pos(int i, int j, int k) { }

	





};

k3d::iplugin_factory& voxel_grid_factory();
	
}

#endif
