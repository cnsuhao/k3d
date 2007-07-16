#include "solver.h"

// semi-largangian method, as it applies to water simulation
// this module is perhaps better off when moved together with voxel grid (ie. voxel grid belongs in the solver plugin, not seperately)
//
// uses pois3d solver from FISHPAK - this works if there are no internal boundaries
// if internal boundaries are present, a multigrid method should be used or a relaxation routine
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

		m_l = vgrid->num_rows(); // i
		m_m = vgrid->num_cols(); // j
		m_n = vgrid->num_slices(); // k
		m_ldimf = m_l;
		m_mdimf = m_m;
		m_ierror = 0;


		int L = m_l;
		int M = m_m;
		int N = m_n;

		int maximum = max(m_l, m_m);
		maximum = max(maximum, m_n);

		m_A = new real[vgrid->num_slices()];
		m_B = new real[vgrid->num_slices()];
		m_C = new real[vgrid->num_slices()];

		// workspace
		m_W = new real[30 + L + M + 2*N + maximum + 7*((L+1)/2 + (M+1)/2)];

		// non-periodic boundary conditions
		m_lperod = 1; // l[0] = l[n] = 0
		m_mperod = 1;
		m_nperod = 1;
	
		// u0* should always hold the previous solution of the velocity field
		float* u0x = vgrid->get_velocity_x(); 
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
			// need to advect the field by itsef for each component of the velcity, since we are using a MAC grid 
			swap(u0x, u1x); 
			vstep(u1x, u0x, u0x, u0y, u0z, vel_x);

			
			swap(u0y, u1y);
			vstep(u1y, u0y, u0x, u0y, u0z, vel_y);

			swap(u0z, u1z);
			vstep(u1z, u0z, u0x, u0y, u0z, vel_z);

			// the last step is the projection ---> divergence free field

			project(u1x, u1y, u1z, u0x, u0y, u0z);
			

			// sstep(s1x, s0x, viscosity, densities); 
			// NOTE: don't need sstep for water simulation - only the velocities are updated using the semi-Lagrangian method,
			// Densities are updates using other methods
		}

		delete work_space_ux;
		delete work_space_uy;
		delete work_space_uz;

		delete m_A;
		delete m_B;
		delete m_C;
		delete m_W;
	}

	void solver::setup_diffusion_constants()
	{
		float dt = m_timestep.value();
		float visc = m_viscosity.value();
		float vox_size = m_voxel_grid.value()->voxel_width();
		int k = m_voxel_grid.value()->num_slices();
		float vox_size_sq = vox_size*vox_size;
		
		// our voxels have the same width, height, and length
		m_K1 = -dt*visc/vox_size_sq;
		m_K2 = m_K1;

		for (int i = 0; i < k; ++i) {
			m_A[i] = m_K1; // only since the height == length == width
			m_C[i] = m_K1;
			m_B[i] = 1.0 + 2.0*dt*visc*vox_size_sq;
		}
	}

	void solver::setup_projection_constants()		
	{
		float dt = m_timestep.value();
		float visc = m_viscosity.value();
		float vox_size = m_voxel_grid.value()->voxel_width();
		int k = m_voxel_grid.value()->num_slices();
		float vox_size_sq = vox_size*vox_size;


		m_K1 = 1.0/vox_size_sq;
		m_K2 = m_K1;
		for (int i = 0; i < k; ++i) {
			m_A[i] = m_K1; // only since the height == length == width
			m_C[i] = m_K1;
			m_B[i] = -2.0/vox_size_sq;
		}


	}

	void solver::vstep(float* new_vfield, float* old_vfield, float *old_x, float* old_y, float *old_z, velocity_component vc)
	{
		// first, add forces
		// addForce...
		//
		// follow the convention from Stable Fluids with u0, u1 to avoid swapping
		float* u0 = old_vfield;
		float *u1 = new_vfield;

		transport(u1, u0, old_x, old_y, old_z, vc);
		setup_diffusion_constants();
		diffuse(u1);

		// proejction step needs to be in the simulate routine
		//setup_projection_constants();
		//project(u1, old_x, old_y, old_z);
	}


	// find the new velocity using the semi-Lagrangian method - ie. trace the particle back in time
	// and interpolate the velocity at that location, and then take that velocity to be the new velocity
	void solver::transport(float* new_vfield, float* old_vfield, float* xfield, float* yfield, float* zfield, velocity_component vc)
	{
		voxel_grid* vgrid = m_voxel_grid.value();
		float x0[vgrid->number_of_voxels()];
		k3d::point3 location;
		k3d::point3 result;

		// since velocities are located at each face, we have to trace each velocity component (ie. 3 velocity components per voxel)
		//

		if (vc == vel_x) {
			for (int i = 0; i < vgrid->num_rows(); ++i) {
				for (int j = 0; j < vgrid->num_cols(); ++j) {
					for (int k = 0; k < vgrid->num_slices(); ++k) {
						location = vgrid->get_velocity_i_pos(i,j,k);
						trace_particle(location, result, xfield, yfield, zfield);
						new_vfield[vgrid->indexof(i,j,k)] = m_voxel_grid.value()->velocity_x(result);
					}
				}
			}
		}
		else if (vc == vel_y) {
			for (int i = 0; i < vgrid->num_rows(); ++i) {
				for (int j = 0; j < vgrid->num_cols(); ++j) {
					for (int k = 0; k < vgrid->num_slices(); ++k) {
						location = vgrid->get_velocity_j_pos(i,j,k);
						trace_particle(location, result, xfield, yfield, zfield);
						new_vfield[vgrid->indexof(i,j,k)] = m_voxel_grid.value()->velocity_y(result);
					}
				}
			}

		}
		else if (vc == vel_z) {
			for (int i = 0; i < vgrid->num_rows(); ++i) {
				for (int j = 0; j < vgrid->num_cols(); ++j) {
					for (int k = 0; k < vgrid->num_slices(); ++k) {
						location = vgrid->get_velocity_k_pos(i,j,k);
						trace_particle(location, result, xfield, yfield, zfield);
						new_vfield[vgrid->indexof(i,j,k)] = m_voxel_grid.value()->velocity_z(result);
					}
				}
			}

		}
		
	}

	// back-trace the particle at position original and store the new position in result
	// we take this velocity to be the new velocity
	void solver::trace_particle(const k3d::point3& original, k3d::point3& result, float* xfield, float* yfield, float* zfield)
	{
		float dt = m_timestep.value();
		voxel_grid* grid = m_voxel_grid.value();

		int i = grid->get_i(original);
		int j = grid->get_j(original);
		int k = grid->get_k(original);

		result[0] = original[0] - dt*xfield[grid->indexof(i,j,k)];
		result[1] = original[1] - dt*yfield[grid->indexof(i,j,k)];
		result[2] = original[2] - dt*zfield[grid->indexof(i,j,k)];

	}


	// use the values from appendix B in SF for diffuse and project using pois3d form FISHPAK
	// LHS = U, RHS = F = old solution
	// the answer is stored in F
	void solver::diffuse(float* F)
	{
		// with non-periodic boundary conditions
		pois3d_(&m_lperod, &m_l, &m_K1, &m_mperod, &m_m, &m_K2, &m_nperod, &m_n, m_A, m_B, m_C, &m_ldimf, &m_mdimf, F, &m_ierror, m_W);
	}


	// projection step - u0* is the old velocity field
	// makes the field divergence free, tus enforcing equation 2 from Stable Fluids
	void solver::project(float* u1x, float* u1y, float* u1z, float* u0x, float* u0y, float* u0z)
	{
		voxel_grid* grid = m_voxel_grid.value();
		int rows = grid->num_rows();
		int cols = grid->num_cols();
		int slices = grid->num_slices();
		float vox_width = grid->voxel_width();
		
		float F[rows*cols*slices];
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				for (int k = 0; k < slices; ++k) {
					// RHS = divergence of the velocity field U
					F[grid->indexof(i,j,k)] = 0.5*((u0x[grid->indexof(i+1,j,k)]-u0x[grid->indexof(i-1,j,k)])/vox_width +
							         (u0y[grid->indexof(i,j+1,k)]-u0y[grid->indexof(i,j-1,k)])/vox_width +
								 (u0z[grid->indexof(i,j,k+1)]-u0z[grid->indexof(i,j,k-1)])/vox_width);


				}
			}
		}

		// with non-periodic boundary conditions - the outer voxels are taken to be solid (ie. walls)
		pois3d_(&m_lperod, &m_l, &m_K1, &m_mperod, &m_m, &m_K2, &m_nperod, &m_n, m_A, m_B, m_C, &m_ldimf, &m_mdimf, F, &m_ierror, m_W);

		// subtract the gradient of the solution from the pervious solution U
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				for (int k = 0; k < slices; ++k) {
					u1x[grid->indexof(i,j,k)] -= 0.5*(F[grid->indexof(i+1,j,k)] - F[grid->indexof(i-1,j,k)])/vox_width;
					u1y[grid->indexof(i,j,k)] -= 0.5*(F[grid->indexof(i,j+1,k)] - F[grid->indexof(i,j-1,k)])/vox_width;
					u1z[grid->indexof(i,j,k)] -= 0.5*(F[grid->indexof(i,j,k+1)] - F[grid->indexof(i,j,k-1)])/vox_width;
				}
			}
		}



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
