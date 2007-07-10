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

		m_grid = new array_type(boost::extents[m_rows.value()][m_cols.value()][m_slices.value()]);
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
		double diff = std::fabs(m_orig_px - m_orig_nx);
		double rem = fmod(diff, m_vox_width.value());
		if (rem != 0) {
			double offset = rem/2.0;
			m_nx.set_value(m_orig_nx - offset);
			m_px.set_value(m_orig_px + offset);
		}
		m_cols.set_value((int)((std::fabs(m_px.value() - m_nx.value()))/m_vox_width.value()));
		
	}

	void voxel_grid::vgrid_modified_y(k3d::iunknown* Hint) {
		double diff = std::fabs(m_orig_py - m_orig_ny);
		double rem = fmod(diff, m_vox_width.value());
		if (rem != 0) {
			double offset = rem/2.0;
			m_ny.set_value(m_orig_ny - offset);
			m_py.set_value(m_orig_py + offset);
		}

		m_rows.set_value((int)((std::fabs(m_py.value() - m_ny.value()))/m_vox_width.value()));
	}

	void voxel_grid::vgrid_modified_z(k3d::iunknown* Hint) {
		double diff = std::fabs(m_orig_pz - m_orig_nz);
		double rem = fmod(diff, m_vox_width.value());
		if (rem != 0) {
			double offset = rem/2.0;
			m_nz.set_value(m_orig_nz - offset);
			m_pz.set_value(m_orig_pz + offset);
		}

		m_slices.set_value((int)((std::fabs(m_pz.value() - m_nz.value()))/m_vox_width.value()));
	}


	void voxel_grid::vgrid_modified(k3d::iunknown* Hint) {
		vgrid_modified_x(Hint);
		vgrid_modified_y(Hint);
		vgrid_modified_z(Hint);
		voxel_grid_changed_signal.emit(Hint);
	}


	k3d::iplugin_factory& voxel_grid_factory() {
		return voxel_grid::get_factory();
	}

	// tri-linear interpolation - adapted from Graphics Gems IV
	k3d::vector3 voxel_grid::velocity(const k3d::point3& pos)
	{
		typedef k3d::point3 point3;
		typedef k3d::vector3 vector3;

		vector3 vel;
		double v000, v100, v010, v110, v001, v101, v011, v111;
		double vi00, vi01, vi10, vi11, vij0, vij1, vijk;
		double fi, fj, fk;
		double vox_width = m_vox_width.value();

		// cell location
		point3 loc = (pos - m_norigin)/vox_width;

		int i = (int)loc[0];
		int j = (int)loc[1];
		int k = (int)loc[2];

		vector3 p0 = n_coord(i,j,k);


		// we need four voxels in order to  sample 8 velocities, so that we we can
		// trilinearly interpolate the velocity at pos
		// ie.  each voxel has 2 velocities (ie. one velocity is shared with its neighbour)


		// first interpolate the i component of the velocity (ie. the component along the x-axis)
		// 4 cases
		
		// locations of j and k velocity components on the voxel
		double posk = p0[3] + 0.5*vox_width;
		double posj = p0[2] + 0.5*vox_width;

		if (pos[3] > posk) {
			if (pos[2] > posj) {
				v000 = get_velocity_i(i,j,k);
				v100 = get_velocity_i(i+1,j,k);
				v110 = get_velocity_i(i+1,j+1,k);
				v010 = get_velocity_i(i,j+1,k);

				v001 = get_velocity_i(i,j,k+1);
				v101 = get_velocity_i(i+1,j,k+1);
				v111 = get_velocity_i(i+1,j+1,k+1);
				v011 = get_velocity_i(i,j+1,k+1);

			}
			else {
				v000 = get_velocity_i(i,j-1,k);
				v100 = get_velocity_i(i+1,j-1,k);
				v110 = get_velocity_i(i+1,j,k);
				v010 = get_velocity_i(i,j,k);

				v001 = get_velocity_i(i,j-1,k+1);
				v101 = get_velocity_i(i+1,j-1,k+1);
				v111 = get_velocity_i(i+1,j,k+1);
				v011 = get_velocity_i(i,j,k+1);
			}
		}
		else { // if pos[3] <= posk
			if (pos[2] > posj) {
				v000 = get_velocity_i(i,j,k-1);
				v100 = get_velocity_i(i+1,j,k-1);
				v110 = get_velocity_i(i+1,j+1,k-1);
				v010 = get_velocity_i(i,j+1,k-1);

				v001 = get_velocity_i(i,j,k);
				v101 = get_velocity_i(i+1,j,k);
				v111 = get_velocity_i(i+1,j+1,k);
				v011 = get_velocity_i(i,j+1,k);
			}
			else {
				v000 = get_velocity_i(i,j-1,k-1);
				v100 = get_velocity_i(i+1,j-1,k-1);
				v110 = get_velocity_i(i+1,j,k-1);
				v010 = get_velocity_i(i,j,k-1);

				v001 = get_velocity_i(i,j-1,k);
				v101 = get_velocity_i(i+1,j-1,k);
				v111 = get_velocity_i(i+1,j,k);
				v011 = get_velocity_i(i,j,k);
			}
		}



		// fi, fj, and fk (ie. fractions of the distance) must be computed!
		vi00 = lerp(fi, v000, v100);
		vi01 = lerp(fi, v001, v101);
		vi10 = lerp(fi, v010, v110);
		vi11 = lerp(fi, v011, v111);

		vij0 = lerp(fj, vi00, vi10);
		vij1 = lerp(fj, vi01, vi11);

		vijk = lerp(fk, vij0, vij1); // this is the x-component (ie. the interpolated i component of the evelocity);
		vel[1] = vijk;

		// DO THE SAME for j, k components


	





		
	}

}
