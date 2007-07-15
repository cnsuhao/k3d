#include "solver.h"

// semi-largangian, so far, as it applies to water simulation

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
		return sigc::mem_fun(*this, &solver::start_solver);
	}

	void solver::start_solver(k3d::iunknown* const Hint) 
	{
		simulate();
	}

	void solver::simulate()
	{
		voxel_grid* vgrid = m_voxel_grid.value();
		int num_voxels = vgrid->number_of_voxels();

		float* u0x = vgrid->get_velocity_x(); // this is wrong - don't need u0x, etc...just u0 
						      // which has the size num_rows*num_cols*num_slices
						      // chnage voxel_grid to reflect that
		float* u0y = vgrid->get_velocity_y();
		float* u0z = vgrid->get_velocity_z();

		// second grid
		float* work_space_ux = new float[num_voxels];
		float* work_space_uy = new float[num_voxels];
		float* work_space_uz = new float[num_voxels];

		float* u1x = work_space_ux;
		float* u1y = work_space_uy;
		float* u1z = work_space_uz;

		float viscosity = m_viscosity.value();

		float* forces; // = make it point to forces

		for (int i = 0; i < m_steps.value(); ++i) {
			swap(u0x, u1x); 
			vstep(u1x, u0x, viscosity, forces);

			swap(u0y, u1y);
			vstep(u1y, u0y, viscosity, forces);

			swap(u0z, u1z);
			vstep(u1z, u0z, viscosity, forces);

			// sstep(s1x, s0x, viscosity, densities); 
			// NOTE: don't need sstep for water animation - only the velocities are updated using the semi-Lagrangian method
		}

		delete work_space_ux;
		delete work_space_uy;
		delete work_space_uz;
	}

	void solver::vstep(float* u1, float* u0, float viscosity, float* forces)
	{
		// first, add forces
		// addForce...
		transport(u1, u0, u0);
		diffuse(u0, u1);
		project(u1, u0);
	}


	void solver::transport(float* u1, float* u0, float* ufield)
	{
		voxel_grid* vgrid = m_voxel_grid.value();
		float x0[vgrid->number_of_voxels()];
		k3d::vector3 location;

		for (int i = 0; i < vgrid->num_rows(); ++i) {
			for (int j = 0; j < vgrid->num_cols(); ++j) {
				for (int k = 0; k < vgrid->num_slices(); ++k) {
					// must do for velocity component, since they are stored on the faces
					location = m_voxel_grid.value()->get_velocity_i_pos(i,j,k);
					trace_particle(location, ufield, -m_timestep.value(), x0);
					//u1x = lerp(location, u0);

					location = m_voxel_grid.value()->get_velocity_j_pos(i,j,k);
					trace_particle(location, ufield, -m_timestep.value(), x0);
					// u1y = lerp(location, u0);

					location = m_voxel_grid.value()->get_velocity_k_pos(i,j,k);
					trace_particle(location, ufield, -m_timestep.value(), x0);
					// u1z = lerp(location, u0);

	
				}
			}
		}


	}

	// Runge-Kutta 2nd order - implement the modified Euler method
	void solver::trace_particle(const k3d::vector3& x, float* ufield, float timestep, float* result)
	{

	}


	// figure out pois3d from fishpak, and use the values from appendix B in SF for diffuse and project
	void solver::diffuse(float* u1, float* u0)
	{
				
	}

	void solver::project(float* u1, float* u0)
	{

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





}
