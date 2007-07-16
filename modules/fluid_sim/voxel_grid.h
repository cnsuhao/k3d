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

// The grid doesn't use multi_array to represent voxels anymore
// rather, each velocity component is stored in a float* array, so that it can be passed onto pois3d to handle
// dissipation and projection steps


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
	//const cell& get_cell(const k3d::point3& pos);
	int get_i(const k3d::point3& pos) { return (int)(pos[0] - m_nx.value())/m_num_cols; }
	int get_j(const k3d::point3& pos) { return (int)(pos[1] - m_ny.value())/m_num_rows; }
	int get_k(const k3d::point3& pos) { return (int)(pos[2] - m_nz.value())/m_num_slices; }

	// should NOT be doing this, but arrays are needed in the fluid solver in order to use POIS3D from FISHPAK
	float* get_velocity_x() { return m_velocity_x; }
	float* get_velocity_y() { return m_velocity_y; }
	float* get_velocity_z() { return m_velocity_z; }

	float nx() { return m_nx.value(); }
	float ny() { return m_ny.value(); }
	float nz() { return m_nz.value(); }

	float px() { return m_px.value(); }
	float py() { return m_py.value(); }
	float pz() { return m_pz.value(); }

	/* i = rows, j = columns, k = slice */
	float get_velocity_i(int i, int j, int k) { return m_velocity_x[m_num_rows*m_num_cols*k + m_num_cols*i + j]; }
	float get_velocity_j(int i, int j, int k) { return m_velocity_y[m_num_rows*m_num_cols*k + m_num_cols*i + j]; }
	float get_velocity_k(int i, int j, int k) { return m_velocity_y[m_num_rows*m_num_cols*k + m_num_cols*i + j]; }


	// return the location of the specified velocity component for voxel (i,j,k)
	k3d::point3 get_velocity_i_pos(int i, int j, int k) { }
	k3d::point3 get_velocity_j_pos(int i, int j, int k) { }
	k3d::point3 get_velocity_k_pos(int i, int j, int k) { }

	int number_of_voxels() { return m_num_rows * m_num_cols * m_num_slices; }
	int num_rows() { return m_num_rows; }
	int num_cols() { return m_num_cols; }
	int num_slices() { return m_num_slices; }

	int indexof(int i, int j, int k) { 
		if (i < 0) i = 0;
		else if (i > m_num_rows) i = m_num_rows;
		if (j < 0) j = 0;
		else if (j > m_num_cols) j = m_num_cols;
		if (k < 0) k = 0;
		else if (k > m_num_slices) k = m_num_slices;

		return m_num_rows*m_num_cols*k + m_num_cols*i + j; 
	}

	// interpolate different velocity components
	float velocity_x(const k3d::point3& pos);
	float velocity_y(const k3d::point3& pos);
	float velocity_z(const k3d::point3& pos);

	float voxel_width() { return m_vox_width.value(); }

private:
	typedef boost::multi_array<cell, 3> array_type;


	// many of the properties will be removed - only voxel width will be modifiable by the user
	k3d_data(double, immutable_name, change_signal, with_undo, local_storage, no_constraint, read_only_property, no_serialization) m_nx;
        k3d_data(double, immutable_name, change_signal, with_undo, local_storage, no_constraint, read_only_property, no_serialization) m_ny;
	k3d_data(double, immutable_name, change_signal, with_undo, local_storage, no_constraint, read_only_property, no_serialization) m_nz;

	k3d_data(double, immutable_name, change_signal, with_undo, local_storage, no_constraint, read_only_property, no_serialization) m_px;
        k3d_data(double, immutable_name, change_signal, with_undo, local_storage, no_constraint, read_only_property, no_serialization) m_py;
	k3d_data(double, immutable_name, change_signal, with_undo, local_storage, no_constraint, read_only_property, no_serialization) m_pz;


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
	//array_type* m_grid;

	// voxel grid  - must be changed again to just float* m_velocity which has the size of m_rows*m_cols*m_slices in order to work with
	// FISHPAK's pois3d function when converted from Fortran to C with f2c
	float* m_velocity_x;
	float* m_velocity_y;
	float* m_velocity_z;

	float* m_density;

	float m_visc;

	int m_num_rows;
	int m_num_cols;
	int m_num_slices;


	double lerp(double fx, double v0, double v1) { return v0 + fx*(v1-v0); }

	// (nx, ny, nz) of the voxel (i,j,k)
	k3d::vector3 n_coord(int i, int j, int k) { return m_norigin + k3d::vector3(i, j, k) * m_vox_width.value(); }

};

k3d::iplugin_factory& voxel_grid_factory();
	
}

#endif
