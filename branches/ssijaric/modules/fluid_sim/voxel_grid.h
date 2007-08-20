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

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include <boost/generator_iterator.hpp>


#include "types.h"


namespace fluid_sim
{
	typedef boost::multi_array<float, 3> array3d_f; // float
	typedef boost::multi_array<int, 3> array3d_i;
	typedef boost::minstd_rand base_generator_type;


// the voxel grid is specifed in similar fashion as the k3d bounding box - with (nx, ny, nz) and (px, py, pz)
class voxel_grid : public k3d::node
{
public:
	voxel_grid(k3d::iplugin_factory& Factory, k3d::idocument& Document); 
	voxel_grid(const voxel_grid& grid);
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

	float nx() { return m_nx.value(); }
	float ny() { return m_ny.value(); }
	float nz() { return m_nz.value(); }

	float px() { return m_px.value(); }
	float py() { return m_py.value(); }
	float pz() { return m_pz.value(); }

	int number_of_voxels() const { return m_xvox*m_yvox*m_zvox; }

	int xvoxels() const { return m_xvox; }
	int yvoxels() const { return m_yvox; }
	int zvoxels() const { return m_zvox; }

	int xfaces() const { return m_xfaces; }
	int yfaces() const { return m_yfaces; }
	int zfaces() const { return m_zfaces; }

	float vx(int i, int j, int k) const { return (*m_grid_vx)[i][j][k]; }
	float vy(int i, int j, int k) const { return (*m_grid_vy)[i][j][k]; }
	float vz(int i, int j, int k) const { return (*m_grid_vz)[i][j][k]; }

	float& vx(int i, int j, int k) { return (*m_grid_vx)[i][j][k]; }
	float& vy(int i, int j, int k) { return (*m_grid_vy)[i][j][k]; }
	float& vz(int i, int j, int k) { return (*m_grid_vz)[i][j][k]; }


	float vx(const idx& l) const { return (*m_grid_vx)[l.i][l.j][l.k]; }
	float vy(const idx& l) const { return (*m_grid_vy)[l.i][l.j][l.k]; }
	float vz(const idx& l) const { return (*m_grid_vz)[l.i][l.j][l.k]; }


	float& vx(const idx& l) { return (*m_grid_vx)[l.i][l.j][l.k]; }
	float& vy(const idx& l) { return (*m_grid_vy)[l.i][l.j][l.k]; }
	float& vz(const idx& l) { return (*m_grid_vz)[l.i][l.j][l.k]; }

	k3d::point3 lower_voxel_corner(int i, int j, int k) { return k3d::point3(m_nx.value() + m_voxel_width*i, m_ny.value() + m_voxel_width*j,
									     m_nz.value() + m_voxel_width*k); }

	k3d::point3 loc_vx(int i, int j, int k) { return k3d::point3(m_nx.value() + m_voxel_width*i, m_ny.value() + m_voxel_width*(j+0.5),
			                                            m_nz.value() + m_voxel_width*(k+0.5)); }

	k3d::point3 loc_vy(int i, int j, int k) { return k3d::point3(m_nx.value() + m_voxel_width*(i+0.5), m_ny.value() + m_voxel_width*j,
			                                             m_nz.value() + m_voxel_width*(k+0.5)); }

	k3d::point3 loc_vz(int i, int j, int k) { return k3d::point3(m_nx.value() + m_voxel_width*(i+0.5), m_ny.value() + m_voxel_width*(j + 0.5),
			                                             m_nz.value() + m_voxel_width*k); }





	void make_obstacle(int i, int j, int k) {
		(*m_vox_type)[i][j][k] = OBSTACLE;
		(*m_grid_vx)[i][j][k] = 0;
		(*m_grid_vy)[i][j][k] = 0;
		(*m_grid_vz)[i][j][k] = 0;
	}

	k3d::point3 random_location_in_cell(int i, int j, int k);	
	void add_particle(int i, int j, int k);

	// interpolate different velocity components
	
	float interpolate_vx(const k3d::point3& pos);
	float interpolate_vy(const k3d::point3& pos);
	float interpolate_vz(const k3d::point3& pos);

	float interpolate_vx(float x, float y, float z);
	float interpolate_vy(float x, float y, float z);
	float interpolate_vz(float x, float y, float z);

	void setup_outer_voxel_layer();

	void mark_cell_as_fluid(const k3d::point3& p);
	void fluid_to_empty();


	float voxel_width() const { return m_voxel_width; }

	enum voxel_type {
		FLUID, 
		OBSTACLE,
		EMPTY // air
	};

	enum velocity_type {
		VX,
		VY,
		VZ
	};

	float interpolate(float x, float y, float z, velocity_type vtype);

	
	//typedef boost::multi_array<float, 3> array3d_f; // float
	//typedef boost::multi_array<int, 3> array3d_i;
	typedef boost::multi_array<voxel_type,3> array3d_type; // type

	bool is_solid(int i, int j, int k) const { return ((*m_vox_type)[i][j][k] == OBSTACLE) ? true : false; }
	bool is_fluid(int i, int j, int k) const { return ((*m_vox_type)[i][j][k] == FLUID) ? true : false; }
	bool is_air(int i, int j, int k) const { return ((*m_vox_type)[i][j][k] == EMPTY) ? true : false; }

	float& pressure(int i, int j, int k) { return (*m_pressure)[i][j][k]; }
	float& density(int i, int j, int k) { return (*m_density)[i][j][k]; }

	float pressure(int i, int j, int k) const { return (*m_pressure)[i][j][k]; }
	float density(int i, int j, int k) const { return (*m_density)[i][j][k]; }

	voxel_type vox_type(int i, int j, int k) const { return (*m_vox_type)[i][j][k]; }
	voxel_type& vox_type(int i, int j, int k) { return (*m_vox_type)[i][j][k]; }

private:
	// many of the properties will be removed - only voxel width will be modifiable by the user
	k3d_data(float, immutable_name, change_signal, with_undo, local_storage, no_constraint, read_only_property, no_serialization) m_nx;
        k3d_data(float, immutable_name, change_signal, with_undo, local_storage, no_constraint, read_only_property, no_serialization) m_ny;
	k3d_data(float, immutable_name, change_signal, with_undo, local_storage, no_constraint, read_only_property, no_serialization) m_nz;

	k3d_data(float, immutable_name, change_signal, with_undo, local_storage, no_constraint, read_only_property, no_serialization) m_px;
        k3d_data(float, immutable_name, change_signal, with_undo, local_storage, no_constraint, read_only_property, no_serialization) m_py;
	k3d_data(float, immutable_name, change_signal, with_undo, local_storage, no_constraint, read_only_property, no_serialization) m_pz;


	k3d_data(float, immutable_name, change_signal, with_undo, local_storage, no_constraint, writable_property, no_serialization) m_vox_width;
	
	k3d_data(int, immutable_name, change_signal, with_undo, local_storage, no_constraint, read_only_property, no_serialization) m_rows;
	k3d_data(int, immutable_name, change_signal, with_undo, local_storage, no_constraint, read_only_property, no_serialization) m_cols;
	k3d_data(int, immutable_name, change_signal, with_undo, local_storage, no_constraint, read_only_property, no_serialization) m_slices;

	float m_orig_nx;
	float m_orig_ny;
	float m_orig_nz;

	float m_orig_px;
	float m_orig_py;
	float m_orig_pz;

	/**
	float m_cur_nx;
	float m_cur_ny;
	float m_cur_nz;

	float m_cur_px;
	float m_cur_py;
	float m_cut_pz;
	**/

	// right now, widthx = widthy = widthz, but this may change in the future
	float m_vox_width_x;
	float m_vox_width_y;
	float m_vox_width_z;

	float m_voxel_width;

	k3d::vector3 m_norigin; // <nx, ny, nz>
	k3d::vector3 m_porigin; // <px, py, pz>

	array3d_f* m_grid_vx;
	array3d_f* m_grid_vy;
	array3d_f* m_grid_vz;

	array3d_f* m_density;
	array3d_f* m_pressure;
	array3d_type* m_vox_type;

	float m_visc;

	// number of components in each direction
	int m_xfaces; 
	int m_yfaces; 
	int m_zfaces; 

	int m_xvox;
	int m_yvox;
	int m_zvox;
};

k3d::iplugin_factory& voxel_grid_factory();
	
}

#endif
