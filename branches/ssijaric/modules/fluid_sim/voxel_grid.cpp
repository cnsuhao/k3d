#include <k3dsdk/document_plugin_factory.h>
#include <k3dsdk/module.h>
#include <k3dsdk/node.h>

#include "voxel_grid.h"


// many of the properties will be removed...only voxel width will be modifiable by the user
namespace fluid_sim
{
	voxel_grid::voxel_grid(k3d::iplugin_factory& Factory, k3d::idocument& Document) :
		k3d::node(Factory, Document),
		m_nx(init_owner(*this) + init_name("nx") + init_label(_("Nx")) + init_description(_("Nx")) + init_value(-5) + init_step_increment(1) + init_units(typeid(k3d::measurement::scalar))),
		m_ny(init_owner(*this) + init_name("ny") + init_label(_("Ny")) + init_description(_("Ny")) + init_value(-5) + init_step_increment(1) + init_units(typeid(k3d::measurement::scalar))),
		m_nz(init_owner(*this) + init_name("nz") + init_label(_("Nz")) + init_description(_("Nz")) + init_value(-5) + init_step_increment(1) + init_units(typeid(k3d::measurement::scalar))),
		m_px(init_owner(*this) + init_name("px") + init_label(_("Px")) + init_description(_("Px")) + init_value(5) + init_step_increment(1) + init_units(typeid(k3d::measurement::scalar))),
		m_py(init_owner(*this) + init_name("py") + init_label(_("Py")) + init_description(_("Py")) + init_value(5) + init_step_increment(1) + init_units(typeid(k3d::measurement::scalar))),
		m_pz(init_owner(*this) + init_name("pz") + init_label(_("Pz")) + init_description(_("Pz")) + init_value(5) + init_step_increment(1) + init_units(typeid(k3d::measurement::scalar))),
		m_vox_width(init_owner(*this) + init_name("voxel_width") + init_label(_("Voxel Side Length")) + init_description(_("Voxel Width")) + init_value(0.5) + init_step_increment(1) + init_units(typeid(k3d::measurement::distance))),
		m_rows(init_owner(*this) + init_name("rows") + init_label(_("Rows")) + init_description(_("Rows")) + init_value(20) + init_step_increment(1) + init_units(typeid(k3d::measurement::scalar))),
		m_cols(init_owner(*this) + init_name("cols") + init_label(_("Cols")) + init_description(_("Cols")) + init_value(20) + init_step_increment(1) + init_units(typeid(k3d::measurement::scalar))),
		m_slices(init_owner(*this) + init_name("slices") + init_label(_("Slices")) + init_description(_("Slices")) + init_value(20) + init_step_increment(1) + init_units(typeid(k3d::measurement::scalar)))



	{
		m_vox_width.changed_signal().connect(vgrid_modified_slot());

		m_orig_nx = m_nx.value();
		m_orig_ny = m_ny.value();
		m_orig_nz = m_nz.value();

		m_orig_px = m_px.value();
		m_orig_py = m_py.value();
		m_orig_pz = m_pz.value();

		m_norigin[0] = m_nx.value();
		m_norigin[1] = m_ny.value();
		m_norigin[2] = m_nz.value();

		m_porigin[0] = m_px.value();
		m_porigin[1] = m_py.value();
		m_porigin[2] = m_pz.value();

	
		m_xcomps = m_cols.value() + 1;
		m_ycomps = m_slices.value() + 1 ;
		m_zcomps = m_rows.value() + 1;

		m_grid_vx = new array3d_f(m_xcomps, m_ycomps, m_zcomps);
		m_grid_vy = new array3d_f(m_xcomps, m_ycomps, m_zcomps);
		m_grid_vz = new array3d_f(m_xcomps, m_ycomps, m_zcomps);
		m_density = new array3d_f(m_xcomps-1, m_ycomps-1, m_zcomps-1);

		
	}

	k3d::iplugin_factory& voxel_grid::get_factory()
 	{
		static k3d::document_plugin_factory<voxel_grid, k3d::interface_list<k3d::inode> > factory(
				k3d::uuid(0xe418f0ad, 0x534ede85, 0x810d43a1, 0x7ee5fa11),
				"VoxelGridPlugin",
				"Voxel Grid",
				"Fluid",
				k3d::iplugin_factory::EXPERIMENTAL);
		return factory;
	}

	sigc::slot<void, k3d::iunknown*> voxel_grid::vgrid_modified_slot()
	{
		return sigc::mem_fun(*this, &voxel_grid::vgrid_modified);
	}


	sigc::slot<void, k3d::iunknown*> voxel_grid::vgrid_modified_slot_x()
        {
                return sigc::mem_fun(*this, &voxel_grid::vgrid_modified_x);
        }

	sigc::slot<void, k3d::iunknown*> voxel_grid::vgrid_modified_slot_y()
        {
                return sigc::mem_fun(*this, &voxel_grid::vgrid_modified_y);
        }

	sigc::slot<void, k3d::iunknown*> voxel_grid::vgrid_modified_slot_z()
        {
                return sigc::mem_fun(*this, &voxel_grid::vgrid_modified_z);
        }



	// this should really be split into three callback frunctions -one for change in x, y, or z.  That way, not
	// all of m_cols, m_slices and m_rows need to be computed, except for when voxel_width changes
	void voxel_grid::vgrid_modified_x(k3d::iunknown* Hint)
	{
		float diff = std::fabs(m_orig_px - m_orig_nx);
		float rem = fmod(diff, m_vox_width.value());
		if (rem != 0) {
			float offset = rem/2.0;
			m_nx.set_value(m_orig_nx - offset);
			m_px.set_value(m_orig_px + offset);
		}
		m_cols.set_value((int)((std::fabs(m_px.value() - m_nx.value()))/m_vox_width.value()));
		
	}

	void voxel_grid::vgrid_modified_y(k3d::iunknown* Hint) {
		float diff = std::fabs(m_orig_py - m_orig_ny);
		float rem = fmod(diff, m_vox_width.value());
		if (rem != 0) {
			float offset = rem/2.0;
			m_ny.set_value(m_orig_ny - offset);
			m_py.set_value(m_orig_py + offset);
		}

		m_rows.set_value((int)((std::fabs(m_py.value() - m_ny.value()))/m_vox_width.value()));
	}

	void voxel_grid::vgrid_modified_z(k3d::iunknown* Hint) {
		float diff = std::fabs(m_orig_pz - m_orig_nz);
		float rem = fmod(diff, m_vox_width.value());
		if (rem != 0) {
			float offset = rem/2.0;
			m_nz.set_value(m_orig_nz - offset);
			m_pz.set_value(m_orig_pz + offset);
		}

		m_slices.set_value((int)((std::fabs(m_pz.value() - m_nz.value()))/m_vox_width.value()));
	}


	void voxel_grid::vgrid_modified(k3d::iunknown* Hint) {
		vgrid_modified_x(Hint);
		vgrid_modified_y(Hint);
		vgrid_modified_z(Hint);
		m_vox_width_x = m_vox_width.value();
		m_vox_width_y = m_vox_width.value();
		m_vox_width_z = m_vox_width.value();
		voxel_grid_changed_signal.emit(Hint);
	}


	k3d::iplugin_factory& voxel_grid_factory() {
		return voxel_grid::get_factory();
	}

	// tri-linear interpolation for different velcoity components
	
	float voxel_grid::interpolate_vx(const k3d::point3& pos)
	{
		return interpolate_vx(pos[0], pos[1], pos[2]);
	}

	float voxel_grid::interpolate_vy(const k3d::point3& pos)
	{
		return interpolate_vy(pos[0], pos[1], pos[2]);
	}

	float voxel_grid::interpolate_vz(const k3d::point3& pos)
	{
		return interpolate_vz(pos[0], pos[1], pos[2]);
	}

	float voxel_grid::interpolate_vx(float x, float y, float z)
	{
		assert(x >= m_nx.value() && x <= m_px.value() && y >= m_ny.value() && y <= m_py.value() && z >= m_nz.value() && z <= m_pz.value());

		float nx = m_nx.value();
		float ny = m_ny.value();
		float nz = m_nz.value();
		float vox_width = m_vox_width.value();

		float i_diff = x - nx;
		float j_diff = y - ny + 0.5*vox_width;
		float k_diff = z - nz + 0.5*vox_width;

		int i = (int)(i_diff/vox_width);
		int j = (int)(j_diff/vox_width);
		int k = (int)(k_diff/vox_width);

		std::cout << "i = " << i << ", j = " << j << ", k = " << k << std::endl;

		// find the location of v000
		float x0 = nx + i*vox_width;
		float y0 = ny + (j + 0.5)*vox_width;
		float z0 = nz + (k + 0.5)*vox_width;

		std::cout << "(x0, y0, z0) = (" << x0 << ", " << y0 << ", " << z0 << ")" << std::endl;

		float dx = (x - x0)/vox_width;
		float dy = (y - y0)/vox_width;
		float dz = (z - z0)/vox_width;


		if (j_diff < 0 &&  k_diff < 0) { // interpolate over 2 points along i - (i,0,0) and (i+1,0,0)
			std::cout << "j < 0, k < 0" << std::endl;

			float v0 = (*m_grid_vx)(i,0,0);
			float v1 = (*m_grid_vx)(i+1,0,0);

			return v0*(1-dx) + v1*dx;
		}
		else if (j_diff < 0) {
			std::cout << "j < 0" << std::endl;

			float v00 = (*m_grid_vx)(i,0,k);
			float v10 = (*m_grid_vx)(i+1,0,k);
			float v11 = (*m_grid_vx)(i+1,0,k+1);
			float v01 = (*m_grid_vx)(i,0,k+1);

			return v00*(1-dx)*(1-dz) + v10*dx*(1-dz) + v11*dx*dz + v01*(1-dx)*dz;

		}	
		else if (k_diff < 0) {
			std::cout << "k < 0" << std::endl;

			float v00 = (*m_grid_vx)(i,j,0);
			float v10 = (*m_grid_vx)(i+1,j,0);
			float v01 = (*m_grid_vx)(i,j+1,0);
			float v11 = (*m_grid_vx)(i+1,j+1,0); // possibly wrong - check!!

			return v00*(1-dx)*(1-dy) + v10*dx*(1-dy) + v01*(1-dx)*dy + v11*dx*dy;

		}
		else {
			float v000 = (*m_grid_vx)(i,j,k);
			float v001 = (*m_grid_vx)(i,j,k+1);
			float v101 = (*m_grid_vx)(i+1,j,k+1);
			float v100 = (*m_grid_vx)(i+1,j,k);

			float v010 = (*m_grid_vx)(i,j+1,k);
			float v011 = (*m_grid_vx)(i,j+1,k+1);
			float v111 = (*m_grid_vx)(i+1,j+1,k+1);
			float v110 = (*m_grid_vx)(i+1,j+1,k);


			std::cout << v000 << " " << v001 << " " << v101 << " " << v100 << " " << v010 << " " <<
				v011 << " " << v111 << " " << v110 << std::endl;

			return v000*(1-dx)*(1-dy)*(1-dz) +
				v100*dx*(1-dy)*(1-dz) +
				v010*(1-dx)*dy*(1-dz) +
				v001*(1-dx)*(1-dy)*dz +
				v101*dx*(1-dy)*dz +
				v011*(1-dx)*dy*dz +
				v110*dx*dy*(1-dz) +
				v111*dx*dy*dz;
		}

	}

	float voxel_grid::interpolate_vy(float x, float y, float z)
	{
		//assert(x >= m_nx && x <= m_px && y >= m_ny && y <= m_py && z >= m_nz && z <= m_pz);
		//
		float nx = m_nx.value();
		float ny = m_ny.value();
		float nz = m_nz.value();
		float vox_width = m_vox_width.value();


		float i_diff = x - (nx + 0.5*vox_width);
		float j_diff = y - ny;
		float k_diff = z - (nz + 0.5*vox_width);

		int i = (int)(i_diff/vox_width);
		int j = (int)(j_diff/vox_width);
		int k = (int)(k_diff/vox_width);

		std::cout << "i = " << i << ", j = " << j << ", k = " << k << std::endl;

		// find the location of v000
		float x0 = nx + (i + 0.5)*vox_width;
		float y0 = ny + j*vox_width;
		float z0 = nz + (k + 0.5)*vox_width;

		std::cout << "(x0, y0, z0) = (" << x0 << ", " << y0 << ", " << z0 << ")" << std::endl;

		float dx = (x - x0)/vox_width;
		float dy = (y - y0)/vox_width;
		float dz = (z - z0)/vox_width;


		if (i_diff < 0 &&  k_diff < 0) { // interpolate over 2 points along i - (i,0,0) and (i+1,0,0)
			std::cout << "j < 0, k < 0" << std::endl;
			float v0 = (*m_grid_vy)[(0,j,0)];
			float v1 = (*m_grid_vy)[(0,j+1,0)];

			return v0*(1-dy) + v1*dy;
		}
		else if (i_diff < 0) {
			std::cout << "i < 0" << std::endl;

			float v00 = (*m_grid_vy)(0,j,k);
			float v10 = (*m_grid_vy)(0,j+1,k);
			float v11 = (*m_grid_vy)(0,j+1,k+1);
			float v01 = (*m_grid_vy)(0,j,k+1);

			return v00*(1-dy)*(1-dz) + v10*dy*(1-dz) + v11*dy*dz + v01*(1-dy)*dz;

		}	
		else if (k_diff < 0) {
			std::cout << "k < 0" << std::endl;

			float v00 = (*m_grid_vy)(i,j,0);
			float v10 = (*m_grid_vy)(i+1,j,0);
			float v01 = (*m_grid_vy)(i,j+1,0);
			float v11 = (*m_grid_vy)(i+1,j+1,0);

			return v00*(1-dx)*(1-dy) + v10*dx*(1-dy) + v01*(1-dx)*dy + v11*dx*dy;

		}
		else {
			float v000 = (*m_grid_vy)(i,j,k);
			float v001 = (*m_grid_vy)(i,j,k+1);
			float v101 = (*m_grid_vy)(i+1,j,k+1);
			float v100 = (*m_grid_vy)(i+1,j,k);

			float v010 = (*m_grid_vy)(i,j,k+1);
			float v011 = (*m_grid_vy)(i,j+1,k+1);
			float v111 = (*m_grid_vy)(i+1,j+1,k+1);
			float v110 = (*m_grid_vy)(i+1,j+1,k);


			std::cout << v000 << " " << v001 << " " << v101 << " " << v100 << " " << v010 << " " <<
				v011 << " " << v111 << " " << v110 << std::endl;

			return v000*(1-dx)*(1-dy)*(1-dz) +
				v100*dx*(1-dy)*(1-dz) +
				v010*(1-dx)*dy*(1-dz) +
				v001*(1-dx)*(1-dy)*dz +
				v101*dx*(1-dy)*dz +
				v011*(1-dx)*dy*dz +
				v110*dx*dy*(1-dz) +
				v111*dx*dy*dz;
		}

	}

	float voxel_grid::interpolate_vz(float x, float y, float z)
	{
		//assert(x >= m_nx && x <= m_px && y >= m_ny && y <= m_py && z >= m_nz && z <= m_pz);

		float nx = m_nx.value();
		float ny = m_ny.value();
		float nz = m_nz.value();
		float vox_width = m_vox_width.value();


		float i_diff = x - (nx + 0.5*vox_width);
		float j_diff = y - (ny + 0.5*vox_width);
		float k_diff = z - nz;

		int i = (int)(i_diff/vox_width);
		int j = (int)(j_diff/vox_width);
		int k = (int)(k_diff/vox_width);

		std::cout << "i = " << i << ", j = " << j << ", k = " << k << std::endl;

		// find the location of v000
		float x0 = nx + (i + 0.5)*vox_width;
		float y0 = ny + (j + 0.5)*vox_width;
		float z0 = nz + k*vox_width; 

		std::cout << "(x0, y0, z0) = (" << x0 << ", " << y0 << ", " << z0 << ")" << std::endl;

		float dx = (x - x0)/vox_width;
		float dy = (y - y0)/vox_width;
		float dz = (z - z0)/vox_width;


		if (i_diff < 0 &&  j_diff < 0) { // interpolate over 2 points along i - (i,0,0) and (i+1,0,0)
			std::cout << "j < 0, k < 0" << std::endl;
			float v0 = (*m_grid_vz)(0,0,k);
			float v1 = (*m_grid_vz)(0,0,k+1);

			return v0*(1-dz) + v1*dz;
		}
		else if (i_diff < 0) {
			std::cout << "i < 0" << std::endl;
			float v00 = (*m_grid_vz)(0,j,k);
			float v10 = (*m_grid_vz)(0,j+1,k);
			float v11 = (*m_grid_vz)(0,j+1,k+1);
			float v01 = (*m_grid_vz)(0,j,k+1);

			return v00*(1-dy)*(1-dz) + v10*dy*(1-dz) + v11*dy*dz + v01*(1-dy)*dz;

		}	
		else if (j_diff < 0) {
			std::cout << "k < 0" << std::endl;
			float v00 = (*m_grid_vz)(i,0,k);
			float v10 = (*m_grid_vz)(i+1,0,k);
			float v01 = (*m_grid_vz)(i,0,k+1);
			float v11 = (*m_grid_vz)(i+1,0,k+1);

			return v00*(1-dx)*(1-dz) + v10*dx*(1-dz) + v01*(1-dx)*dz + v11*dx*dz;

		}
		else {
			float v000 = (*m_grid_vz)(i,j,k);
			float v001 = (*m_grid_vz)(i,j,k+1);
			float v101 = (*m_grid_vz)(i+1,j,k+1);
			float v100 = (*m_grid_vz)(i+1,j,k);

			float v010 = (*m_grid_vz)(i,j+1,k);
			float v011 = (*m_grid_vz)(i,j+1,k+1);
			float v111 = (*m_grid_vz)(i+1,j+1,k+1);
			float v110 = (*m_grid_vz)(i+1,j+1,k);

			std::cout << v000 << " " << v001 << " " << v101 << " " << v100 << " " << v010 << " " <<
				v011 << " " << v111 << " " << v110 << std::endl;

			return v000*(1-dx)*(1-dy)*(1-dz) +
				v100*dx*(1-dy)*(1-dz) +
				v010*(1-dx)*dy*(1-dz) +
				v001*(1-dx)*(1-dy)*dz +
				v101*dx*(1-dy)*dz +
				v011*(1-dx)*dy*dz +
				v110*dx*dy*(1-dz) +
				v111*dx*dy*dz;
		}
	}



}
