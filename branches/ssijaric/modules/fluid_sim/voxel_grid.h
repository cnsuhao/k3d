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

#include "types.h"


namespace fluid_sim
{

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

	float nx() { return m_nx.value(); }
	float ny() { return m_ny.value(); }
	float nz() { return m_nz.value(); }

	float px() { return m_px.value(); }
	float py() { return m_py.value(); }
	float pz() { return m_pz.value(); }

	int number_of_voxels() { return (m_xcomps-1) * (m_ycomps-1) * (m_zcomps-1); }

	// interpolate different velocity components
	float interpolate_vx(const k3d::point3& pos);
	float interpolate_vy(const k3d::point3& pos);
	float interpolate_vz(const k3d::point3& pos);

	float interpolate_vx(float x, float y, float z);
	float interpolate_vy(float x, float y, float z);
	float interpolate_vz(float x, float y, float z);

	float voxel_width() { return m_vox_width.value(); }

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


	// right now, widthx = widthy = widthz, but this may change in the future
	float m_vox_width_x;
	float m_vox_width_y;
	float m_vox_width_z;

	k3d::vector3 m_norigin; // <nx, ny, nz>
	k3d::vector3 m_porigin; // <px, py, pz>

	array3d_f* m_grid_vx;
	array3d_f* m_grid_vy;
	array3d_f* m_grid_vz;

	array3d_f* m_density;

	float m_visc;

	// number of components in each direction
	int m_xcomps; // components along x axis = columns
	int m_ycomps; // components along y axis = slices
	int m_zcomps; // components along z axis = rows
};

k3d::iplugin_factory& voxel_grid_factory();
	
}

#endif
