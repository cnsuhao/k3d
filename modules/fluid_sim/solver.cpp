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





}
