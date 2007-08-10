#include "solver.h"

// semi-largangian method, as it applies to water simulation
//
// non-periodic boundary conditions

namespace fluid_sim
{

	solver::solver(k3d::iplugin_factory& Factory, k3d::idocument& Document) :
		k3d::node(Factory, Document),
		m_voxel_grid(init_owner(*this) + init_name("voxel_grid") + init_label(_("Voxel Grid")) + init_description(_("Voxel Grid")) + init_value<voxel_grid*>(0)),
		m_timestep(init_owner(*this) + init_name("timestep") + init_label(_("Timestep")) + init_description(_("Timestep")) + init_value(0.1) + init_step_increment(0.1) + init_units(typeid(k3d::measurement::scalar))),
		m_steps(init_owner(*this) + init_name("steps") + init_label(_("Steps")) + init_description(_("Steps")) + init_value(2) + init_step_increment(1) + init_units(typeid(k3d::measurement::scalar))),
		m_viscosity(init_owner(*this) + init_name("timestep") + init_label(_("Timestep")) + init_description(_("Timestep")) + init_value(0.1) + init_step_increment(0.1) + init_units(typeid(k3d::measurement::scalar)))


		
	{
		m_voxel_grid.changed_signal().connect(start_solver_slot());	
	}


	sigc::slot<void, k3d::iunknown*> solver::start_solver_slot()		
	{
		
		//return sigc::mem_fun(*this, &solver::start_solver);
	}



	k3d::iplugin_factory& solver::get_factory()
 	{
		static k3d::document_plugin_factory<solver, k3d::interface_list<k3d::inode> > factory(
				k3d::uuid(0x1ac5578c, 0xa944a032, 0x30384fbd, 0x3aab8ab4),
				"FluidSolverPlugin",
				"Fluid Solver",
				"Fluid",
				k3d::iplugin_factory::EXPERIMENTAL);
		return factory;
	}

	k3d::iplugin_factory& solver_factory() {
		return solver::get_factory();
	}


	// Projection step to make the field divergence free
	// Some notes: For all solid cells (ie. obstacles), the pressure values have the Neumann BC of 0 (ie. no change in pressure)
	// For this to work, we also set the velocities at cell faces of solid objects to 0
	// Why foes this work? ---> the projection step involves the subtraction of the pressure gradient from the velocity,
	// which gives the new velocity
	// In the case of the solid cell, the velocity at the face is 0, and the pressure gradient is 0 -> so the new velocity is
	// also 0
	void solver::project(voxel_grid& new_grid, const voxel_grid& old_grid)
	{
		// easier to use the u, w notation
		voxel_grid& u = new_grid;
		const voxel_grid& w = old_grid;
		int xvoxels = w.xvoxels();
		int yvoxels = w.yvoxels();
		int zvoxels = w.zvoxels();
		array3d<float> b(xvoxels, yvoxels, zvoxels);


		// since velocities are defined at each face, that means that there is one more velocity component in each
		// direction than there are pressure and density components (that are defined at the center of ech voxel)
		// Also, boundary conditions are defined at indices 0 and xvoxels, etc.
		//
		// Note that we negate both b and A in order to make A positive definite
		for (int i = 0; i < xvoxels; ++i) {
			for (int j = 0; j < yvoxels; ++j) {
				for (int k = 0; k < zvoxels; ++k) {
					b(i,j,k) = -(w.vx(i+1,j,k) - w.vx(i,j,k) +
						    w.vy(i,j+1,k) - w.vy(i,j,k) +
						    w.vz(i,j,k+1) - w.vz(i,j,k))/m_vox_width;
						
				}
			}
		}
		
		gmm::row_matrix< gmm::wsvector<float> > A(w.number_of_voxels(), w.number_of_voxels());

		for (int i = 0; i < xvoxels; ++i) {
			for (int j = 0; j < yvoxels; ++j) {
				for (int k = 0; k < zvoxels; ++k) {
					A((i,j,k), (i,j,k)) = 6;

					if (i > 0) {
						if (w.is_solid(i-1,j,k)) { A((i,j,k), (i,j,k)) -= 1; u.pressure(i-1,j,k) = u.pressure(i,j,k); A((i,j,k),(i-1,j,k)) = 0; }
						else A((i,j,k),(i-1,j,k)) = -1;
					}
					else {
						A((i,j,k), (i,j,k)) -= 1;
					}

					if (j > 0) {
						if (w.is_solid(i,j-1,k)) { A((i,j,k), (i,j,k)) -= 1; u.pressure(i,j-1,k) = u.pressure(i,j,k); A((i,j,k),(i,j-1,k)) = 0; }
						else A((i,j,k),(i,j-1,k)) = -1;
					}
					else {
						A((i,j,k), (i,j,k)) -= 1;
					}

					if (k > 0) {
						if (w.is_solid(i,j,k-1)) { A((i,j,k), (i,j,k)) -= 1; u.pressure(i,j,k-1) = u.pressure(i,j,k); A((i,j,k),(i,j,k-1)) = 0; }
						else A((i,j,k),(i,j,k-1)) = -1;
					}
					else {
						A((i,j,k), (i,j,k)) -= 1;
					}

					if (i < m_xvox) {
						if (w.is_solid(i+1,j,k)) { A((i,j,k), (i,j,k)) -= 1; u.pressure(i+1,j,k) = u.pressure(i,j,k); A((i,j,k),(i+1,j,k)) = 0; }
						else A((i,j,k),(i+1,j,k)) = -1;
					}
					else {
						A((i,j,k), (i,j,k)) -= 1;
					}

					if (j < m_yvox) {
						if (w.is_solid(i,j+1,k)) { A((i,j,k), (i,j,k)) -= 1; u.pressure(i,j+1,k) = u.pressure(i,j,k); A((i,j,k),(i,j+1,k)) = 0; }
						else A((i,j,k),(i,j+1,k)) = -1;
					}
					else {
						A((i,j,k), (i,j,k)) -= 1;
					}

					if (k < m_zvox) {
						if (w.is_solid(i,j,k+1)) { A((i,j,k), (i,j,k)) -= 1; u.pressure(i,j,k+1) = u.pressure(i,j,k); A((i,j,k),(i,j,k+1)) = 0; }
						else A((i,j,k),(i,j,k+1)) = -1;
					}
					else {
						A((i,j,k), (i,j,k)) -= 1;
					}

				}
			}
		}

		// now we have -Ax = -b, which is sparse, symmetric and positive definite...solve using CG




	}

	void solver::diffuse(array3d<float>& new_field, const array3d<float>& old_field)
	{
		// assert(new_field.size() == old_field.size());
		gmm::row_matrix< gmm::wsvector<double> > A(old_field.size(), old_field.size());
		float beta = m_viscosity.value()*m_timestep.value()/(m_vox_width*m_vox_width);



		// i = 0
		for (int j = 1; j < m_xvox; ++j) {
			for (int k = 1; k < m_zvox; ++k) {
				A((0,j,k), (0,j-1,k)) = -beta;
				A((0,j,k), (0,j,k-1)) = -beta;
				A((0,j,k), (0,j,k)) = 1 + 6*beta;
				A((0,j,k), (0,j,k+1)) = -beta;
				A((0,j,k), (0,j+1,k)) = -beta;
				A((0,j,k), (1,j,k)) = -beta;
			}
		}

		/*
		int i = m_xcomps-1;
		for (int j = 1; j < m_ycomps-1; ++j) {
			for (int k = 1; k < m_zcomps-1; ++k) {
				A((i,j,k), (i-1,j,k)) = -beta;
				A((i,j,k), (i,j-1,k)) = -beta;
				A((i,j,k), (i,j,k-1)) = -beta;
				A((i,j,k), (i,j,k)) = 1 + 6*beta;
				A((i,j,k), (i,j,k+1)) = -beta;
				A((i,j,k), (i,j+1,k)) = -beta;
			}
		}

		// j = 0
		//int j = 0;
		for (int i = 1; i < m_xcomps-1; ++i) {
			for (int k = 1; k < m_zcomps-1; ++k) {
				A((i,0,k), (i-1,0,k)) = -beta;
				A((i,0,k), (i,0,k-1)) = -beta;
				A((i,0,k), (i,0,k)) = 1 + 6*beta;
				A((i,0,k), (i,0,k+1)) = -beta;
				A((i,0,k), (i,1,k)) = -beta;
				A((i,0,k), (i+1,0,k)) = -beta;
			}
		}

		int j = m_ycomps-1;
		for (int i = 1; i < m_xcomps-1; ++i) {
			for (int k = 1; k < m_zcomps-1; ++k) {
				A((i,j,k), (i-1,j,k)) = -beta;
				A((i,j,k), (i,j-1,k)) = -beta;
				A((i,j,k), (i,j,k-1)) = -beta;
				A((i,j,k), (i,j,k)) = 1 + 6*beta;
				A((i,j,k), (i,j,k+1)) = -beta;
				A((i,j,k), (i+1,j,k)) = -beta;

			}
		}

		for (int i = 1; i < m_xcomps-1; ++i) {
			for (int j = 1; j < m_ycomps-1; ++j) {
				A((i,j,0), (i-1,j,0)) = -beta;
				A((i,j,0), (i,j-1,0)) = -beta;
				A((i,j,0), (i,j,0)) = 1 + 6*beta;
				A((i,j,0), (i,j,1)) = -beta;
				A((i,j,0), (i,j+1,0)) = -beta;
				A((i,j,0), (i+1,j,0)) = -beta;
			}
		}

		int k = m_zcomps-1;	
		for (int i = 1; i < m_xcomps-1; ++i) {
			for (int j = 1; j < m_ycomps-1; ++j) {
				A((i,j,k), (i-1,j,k)) = -beta;
				A((i,j,k), (i,j-1,k)) = -beta;
				A((i,j,k), (i,j,k-1)) = -beta;
				A((i,j,k), (i,j,k)) = 1 + 6*beta;
				A((i,j,k), (i,j+1,k)) = -beta;
				A((i,j,k), (i+1,j,k)) = -beta;

			}
		}

		for (int i = 1; i < m_xcomps-1; ++i) {
			for (int j = 1; j < m_ycomps-1; ++j) {
				for (int k = 1; k < m_zcomps-1; ++k) {
					A((i,j,k), (i-1,j,k)) = -beta;
					A((i,j,k), (i,j-1,k)) = -beta;
					A((i,j,k), (i,j,k-1)) = -beta;
					A((i,j,k), (i,j,k)) = 1 + 6*beta;
					A((i,j,k), (i,j,k+1)) = -beta;
					A((i,j,k), (i,j+1,k)) = -beta;
					A((i,j,k), (i+1,j,k)) = -beta;
				}
			}
		}
		*/

		// solve here using CG from gmm++

	}

	void solver::add_force(array3d<float>& field, const array3d<float>& forces)
	{
		float timestep = m_timestep.value();
		for (int i = 0; i < field.size(); ++i) {
			field[i] += timestep*forces[i];
		}
	}


}
