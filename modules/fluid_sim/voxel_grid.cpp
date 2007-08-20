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


		// number of faces in each direction	
		m_xfaces = m_cols.value() + 1;
		m_yfaces = m_slices.value() + 1 ;
		m_zfaces = m_rows.value() + 1;

		m_xvox = m_cols.value();
		m_yvox = m_slices.value();
		m_zvox = m_rows.value();

		m_grid_vx = new array3d_f(boost::extents[m_xfaces][m_yfaces][m_zfaces]);
		m_grid_vy = new array3d_f(boost::extents[m_xfaces][m_yfaces][m_zfaces]);
		m_grid_vz = new array3d_f(boost::extents[m_xfaces][m_yfaces][m_zfaces]);
		m_density = new array3d_f(boost::extents[m_xvox][m_yvox][m_zvox]);
		m_vox_type = new array3d_type(boost::extents[m_xvox][m_yvox][m_zvox]);
		m_voxel_width = m_vox_width.value();

		setup_outer_voxel_layer();
	}

	voxel_grid::voxel_grid(const voxel_grid& grid) : k3d::node((const_cast<voxel_grid&>(grid)).factory(), (const_cast<voxel_grid&>(grid)).document()),
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


		// number of faces in each direction	
		m_xfaces = m_cols.value() + 1;
		m_yfaces = m_slices.value() + 1 ;
		m_zfaces = m_rows.value() + 1;

		m_xvox = m_cols.value();
		m_yvox = m_slices.value();
		m_zvox = m_rows.value();

		m_grid_vx = new array3d_f(boost::extents[m_xfaces][m_yfaces][m_zfaces]);
		m_grid_vy = new array3d_f(boost::extents[m_xfaces][m_yfaces][m_zfaces]);
		m_grid_vz = new array3d_f(boost::extents[m_xfaces][m_yfaces][m_zfaces]);
		m_density = new array3d_f(boost::extents[m_xvox][m_yvox][m_zvox]);
		m_vox_type = new array3d_type(boost::extents[m_xvox][m_yvox][m_zvox]);
		m_voxel_width = m_vox_width.value();

		setup_outer_voxel_layer();

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
		m_voxel_width = m_vox_width.value();
		voxel_grid_changed_signal.emit(Hint);
	}


	k3d::iplugin_factory& voxel_grid_factory() {
		return voxel_grid::get_factory();
	}


	// the outer voxel layer conists of four solid walls
	void voxel_grid::setup_outer_voxel_layer()
	{
		for (int i = 0; i < m_xvox; ++i) {
			for (int j = 0; j < m_yvox; ++j) {
				(*m_vox_type)[i][j][m_zvox-1] = OBSTACLE;
				(*m_vox_type)[i][j][0] = OBSTACLE;
			}
		}

		for (int i = 0; i < m_xvox; ++i) {
			for (int k = 0; k < m_zvox; ++k) {
				(*m_vox_type)[i][m_yvox-1][k] = OBSTACLE;
				(*m_vox_type)[i][0][k] = OBSTACLE;
			}
		}

		for (int j = 0; j < m_yvox; ++j) {
			for (int k = 0; k < m_zvox; ++k) {
				(*m_vox_type)[m_xvox-1][j][k] = OBSTACLE;
				(*m_vox_type)[0][j][k] = OBSTACLE;
			}
		}


		// set velocities 
		for (int i = 0; i < m_xfaces; ++i) {
			for (int j = 0; j < m_yfaces; ++j) {
				(*m_grid_vx)[i][j][m_zfaces-1] = 0;
				(*m_grid_vy)[i][j][m_zfaces-1] = 0;
				(*m_grid_vz)[i][j][m_zfaces-1] = 0;
			}
		}

		for (int i = 0; i < m_xfaces; ++i) {
			for (int k = 0; k < m_zfaces; ++k) {
				(*m_grid_vx)[i][m_yfaces-1][k] = 0;
				(*m_grid_vy)[i][m_yfaces-1][k] = 0;
				(*m_grid_vz)[i][m_yfaces-1][k] = 0;
			}
		}

		for (int j = 0; j < m_yfaces; ++j) {
			for (int k = 0; k < m_zfaces; ++k) {
				(*m_grid_vx)[m_xfaces-1][j][k] = 0;
				(*m_grid_vy)[m_xfaces-1][j][k] = 0;
				(*m_grid_vz)[m_xfaces-1][k][k] = 0;
			}

		}

		// initially set all other voxels as air
		for (int i = 1; i < m_xvox-1; ++i) {
			for (int j = 1; j < m_yvox-1; ++j) {
				for (int k = 1; k < m_zvox-1; ++k) {
					(*m_vox_type)[i][j][k] = EMPTY;	
				}
			}
		}
		
	}

	k3d::point3 voxel_grid::random_location_in_cell(int i, int j, int k)
	{
		base_generator_type generator(42u);
		boost::uniform_real<> uni_dist(0,1);
		boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);

		// generate three random numbers, one for each coordinate
		float x = uni();
		float y = uni();
		float z = uni();

		return lower_voxel_corner(i,j,k) + k3d::point3(x,y,z)*m_voxel_width;
	}

	// set all fluid cell to air (not the obstacles, however)
	void voxel_grid::fluid_to_empty()
	{
		for (int i = 0; i < m_xvox; ++i) {
			for (int j = 0; j < m_yvox; ++j) {
				for (int k = 0; k < m_zvox; ++k) {
					if ((*m_vox_type)[i][j][k] == FLUID) 
						(*m_vox_type)[i][j][k] = EMPTY;
				}
			}
		}

	}

	void voxel_grid::mark_cell_as_fluid(const k3d::point3& p)
	{
		float nx = m_nx.value();
		float ny = m_ny.value();
		float nz = m_nz.value();
		float vox_width = m_vox_width.value();

		float i_diff = p[0] - nx;
		float j_diff = p[1] - ny;
		float k_diff = p[2] - nz;

		int i = (int)(i_diff/vox_width);
		int j = (int)(j_diff/vox_width);
		int k = (int)(k_diff/vox_width);

		(*m_vox_type)[i][j][k] = FLUID;
	}


	// tri-linear interpolation for different velcoity components

	float voxel_grid::interpolate_vx(const k3d::point3& pos)
	{
		return interpolate(pos[0], pos[1], pos[2], VX);
	}

	float voxel_grid::interpolate_vy(const k3d::point3& pos)
	{
		return interpolate(pos[0], pos[1], pos[2], VY);
	}

	float voxel_grid::interpolate_vz(const k3d::point3& pos)
	{
		return interpolate(pos[0], pos[1], pos[2], VZ);
	}

	float voxel_grid::interpolate_vx(float x, float y, float z)
	{
		return interpolate(x,y,z,VX);
	}

	float voxel_grid::interpolate_vy(float x, float y, float z)
	{
		return interpolate(x,y,z,VY);
	}

	float voxel_grid::interpolate_vz(float x, float y, float z)
	{
		return interpolate(x,y,z,VZ);

	}

	float voxel_grid::interpolate(float x, float y, float z, velocity_type vtype)
	{
		//assert(x >= m_nx && x <= m_px && y >= m_ny && y <= m_py && z >= m_nz && z <= m_pz);

		array3d_f* grid_v;
		if (vtype == VX)
			grid_v = m_grid_vx;
		else if (vtype == VY)
			grid_v = m_grid_vy;
		else
			grid_v = m_grid_vz;

		float nx = m_nx.value();
		float ny = m_ny.value();
		float nz = m_nz.value();
		float vox_width = m_vox_width.value();


		float i_diff = x - nx;
		float j_diff = y - ny;
		float k_diff = z - nz;

		int i = (int)(i_diff/vox_width);
		int j = (int)(j_diff/vox_width);
		int k = (int)(k_diff/vox_width);


		// find the location of v000
		float x0, y0, z0;
		if (vtype == VX) {
			x0 = nx + i*vox_width;
			y0 = ny + (j + 0.5)*vox_width;
			z0 = nz + (k + 0.5)*vox_width;
		}
		else if (vtype == VY) {		
			x0 = nx + (i + 0.5)*vox_width;
			y0 = ny + j*vox_width;
			z0 = nz + (k + 0.5)*vox_width;
		}
		else {
			x0 = nx + (i + 0.5)*vox_width;
			y0 = ny + (j + 0.5)*vox_width;
			z0 = nz + k*vox_width; 
		}


		float dx = (x - x0)/vox_width;
		float dy = (y - y0)/vox_width;
		float dz = (z - z0)/vox_width;


		float v000 = (*grid_v)[i][j][k];
		float v001 = (*grid_v)[i][j][k+1];
		float v101 = (*grid_v)[i+1][j][k+1];
		float v100 = (*grid_v)[i+1][j][k];

		float v010 = (*grid_v)[i][j+1][k];
		float v011 = (*grid_v)[i][j+1][k+1];
		float v111 = (*grid_v)[i+1][j+1][k+1];
		float v110 = (*grid_v)[i+1][j+1][k];


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
