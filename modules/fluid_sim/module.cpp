#include <k3dsdk/module.h>


namespace fluid_sim
{
extern k3d::iplugin_factory& voxel_grid_visual_factory();
extern k3d::iplugin_factory& voxel_grid_factory();
extern k3d::iplugin_factory& solver_factory();
}

K3D_MODULE_START(Registry)
Registry.register_factory(fluid_sim::voxel_grid_factory());
Registry.register_factory(fluid_sim::voxel_grid_visual_factory());
Registry.register_factory(fluid_sim::solver_factory());
K3D_MODULE_END
