#include <k3dsdk/document_plugin_factory.h>
#include <k3dsdk/imaterial.h> 
#include <k3dsdk/node.h>
#include <k3dsdk/persistent.h> 
#include <k3dsdk/material.h>
#include <k3dsdk/material_client.h>
#include <k3dsdk/measurement.h>
#include <k3dsdk/mesh_source.h>
#include <k3dsdk/mesh_modifier.h>
#include <k3dsdk/mesh.h>
#include <k3dsdk/module.h>
#include <k3dsdk/property.h>

#include "voxel_grid.h"


// displays the voxel grid using the array based mesh design 
namespace fluid_sim
{
class voxel_grid_visual :
	public k3d::material_client<k3d::mesh_source<k3d::persistent<k3d::node> > >
{
	typedef k3d::material_client<k3d::mesh_source<k3d::persistent<k3d::node> > > base;

public:
	voxel_grid_visual(k3d::iplugin_factory& Factory, k3d::idocument& Document) : base(Factory, Document),
		 m_voxel_grid(init_owner(*this) + init_name("voxel_grid") + init_label(_("Voxel Grid")) + init_description(_("Voxel Grid")) + init_value<voxel_grid*>(0))

	{
		m_voxel_grid.changed_signal().connect(make_topology_changed_slot());
		m_voxel_grid.changed_signal().connect(make_reset_mesh_slot());
	}

	sigc::slot<void, iunknown*> make_reset_mesh_slot()		
	{
		return sigc::mem_fun(*this, &voxel_grid_visual::reset_mesh);
	}



	static k3d::iplugin_factory& get_factory()
 	{
		static k3d::document_plugin_factory<voxel_grid_visual> factory(
				k3d::uuid(0x6f951c6f, 0x3247f037, 0xd2df0fb6, 0x48a0cde6),
				"VoxelGridVisualPlugin",
				"Voxel Grid Visualization",
				"Fluid",
				k3d::iplugin_factory::EXPERIMENTAL);
		return factory;
	}

	void on_create_mesh_topology(k3d::mesh& Mesh)
	{
		typedef k3d::mesh mesh;
		
		if (m_voxel_grid.value() != 0) {
			
				
			m_rows = boost::any_cast<int>(get_property(*(m_voxel_grid.value()), "rows")->property_value());
			m_columns = boost::any_cast<int>(get_property(*(m_voxel_grid.value()), "cols")->property_value());
			m_slices = boost::any_cast<int>(get_property(*(m_voxel_grid.value()), "slices")->property_value());
			float voxel_width = boost::any_cast<float>(get_property(*(m_voxel_grid.value()), "voxel_width")->property_value());
			
			unsigned long point_rows = m_rows + 1;
			unsigned long point_columns = m_columns + 1;
			unsigned long point_slices = m_slices + 1;

			unsigned long rows = m_rows;
			unsigned long columns = m_columns;
			unsigned long slices = m_slices;
	
			unsigned long num_faces = rows*columns*point_slices + rows*point_columns*slices + point_rows*columns*slices;
		
			boost::shared_ptr<mesh::polyhedra_t> polyhedra(new mesh::polyhedra_t());
			boost::shared_ptr<mesh::indices_t> edge_points(new mesh::indices_t(4 * num_faces));
			boost::shared_ptr<mesh::points_t> points(new mesh::points_t(point_rows*point_columns*point_slices));
			boost::shared_ptr<mesh::indices_t> first_faces(new mesh::indices_t(1, 0));
			boost::shared_ptr<mesh::counts_t> face_counts(new mesh::counts_t(1, num_faces));
			boost::shared_ptr<mesh::polyhedra_t::types_t> types(new mesh::polyhedra_t::types_t(1, mesh::polyhedra_t::POLYGONS));
			boost::shared_ptr<mesh::indices_t> face_first_loops(new mesh::indices_t(num_faces));
			boost::shared_ptr<mesh::indices_t> loop_first_edges(new mesh::indices_t(num_faces));
			boost::shared_ptr<mesh::counts_t> face_loop_counts(new mesh::counts_t(num_faces, 1));
			boost::shared_ptr<mesh::selection_t> face_selection(new mesh::selection_t(num_faces, 0.0));
			boost::shared_ptr<mesh::materials_t> face_materials(new mesh::materials_t(num_faces, m_material.value()));
			boost::shared_ptr<mesh::indices_t> clockwise_edges(new mesh::indices_t(4 * num_faces));
			boost::shared_ptr<mesh::selection_t> edge_selection(new mesh::selection_t(4 * num_faces, 0.0));
			boost::shared_ptr<mesh::selection_t> point_selection(new mesh::selection_t(point_rows * point_columns * point_slices, 0.0));


			Mesh.points = points;

			mesh::indices_t::iterator face_first_loop = face_first_loops->begin();
			mesh::indices_t::iterator loop_first_edge = loop_first_edges->begin();
			mesh::indices_t::iterator edge_point = edge_points->begin();
			mesh::indices_t::iterator clockwise_edge = clockwise_edges->begin();
	
			size_t face_number = 0;

			
			/* phase 1 - front to back ---- from +y to -y */
			for (unsigned long slice = 0; slice != point_slices; ++slice) 
			{
				for (unsigned long row = 0; row != rows; ++row) 
				{
					for (unsigned long column = 0; column != columns; ++column)
					{
						*face_first_loop++ = face_number;

						*loop_first_edge++ = 4 * face_number;

						*edge_point++ = column + (row * point_columns) + slice * (point_rows * point_columns);
						*edge_point++ = column + (row * point_columns) + slice * (point_rows * point_columns) + 1;
						*edge_point++ = column + ((row + 1) * point_columns) + slice * (point_rows * point_columns) + 1;
						*edge_point++ = column + ((row + 1) * point_columns) + slice * (point_rows * point_columns);

						*clockwise_edge++ = (4 * face_number) + 1;
                        			*clockwise_edge++ = (4 * face_number) + 2;
			                        *clockwise_edge++ = (4 * face_number) + 3;
                        			*clockwise_edge++ = (4 * face_number);

						++face_number;
					}

				}
			}

			/* phase 2 - top to bottom ---- from +z to -z */
			for (unsigned long row = 0; row != point_rows; ++row) 
			{
				for (unsigned long slice = 0; slice != slices; ++slice) 
				{
					for (unsigned long column = 0; column != columns; ++column)
					{
						*face_first_loop++ = face_number;
						*loop_first_edge++ = 4 * face_number;

						*edge_point++ = column + (row * point_columns) + slice * (point_rows * point_columns);
						*edge_point++ = column + (row * point_columns) + slice * (point_rows * point_columns) + 1;
						*edge_point++ = column + (row * point_columns) + (slice + 1) * (point_rows * point_columns) + 1;
						*edge_point++ = column + (row * point_columns) + (slice + 1) * (point_rows * point_columns);
							
						*clockwise_edge++ = (4 * face_number) + 1;
                        			*clockwise_edge++ = (4 * face_number) + 2;
			                        *clockwise_edge++ = (4 * face_number) + 3;
                        			*clockwise_edge++ = (4 * face_number);

						++face_number;
					}

				}
			}

			/* phase 3 */
			for (unsigned long column = 0; column != point_columns; ++column) 
			{
				for (unsigned long row = 0; row != rows; ++row) 
				{
					for (unsigned long slice = 0; slice != slices; ++slice)
					{
						*face_first_loop++ = face_number;
						*loop_first_edge++ = 4 * face_number;

						*edge_point++ = column + (row * point_columns) + slice * (point_rows * point_columns);
						*edge_point++ = column + (row * point_columns) + (slice + 1) * (point_rows * point_columns);
						*edge_point++ = column + ((row + 1) * point_columns) + (slice + 1) * (point_rows * point_columns);
						*edge_point++ = column + ((row + 1) * point_columns) + slice * (point_rows * point_columns);
						

						*clockwise_edge++ = (4 * face_number) + 1;
                        			*clockwise_edge++ = (4 * face_number) + 2;
			                        *clockwise_edge++ = (4 * face_number) + 3;
                        			*clockwise_edge++ = (4 * face_number);

						++face_number;
					}

				}
			}

			polyhedra->first_faces = first_faces;
			polyhedra->face_counts = face_counts;
			polyhedra->types = types;
			polyhedra->face_first_loops = face_first_loops;
			polyhedra->face_loop_counts = face_loop_counts;
			polyhedra->face_selection = face_selection;
			polyhedra->face_materials = face_materials;
			polyhedra->loop_first_edges = loop_first_edges;
			polyhedra->edge_points = edge_points;
			polyhedra->clockwise_edges = clockwise_edges;
			polyhedra->edge_selection = edge_selection;

			Mesh.polyhedra = polyhedra;
			Mesh.points = points;
			Mesh.point_selection = point_selection;



		}
		
		

	}

	void on_update_mesh_geometry(k3d::mesh& Mesh)
	{	
		// create points here

		if (m_voxel_grid.value() != 0) {
			const float px = boost::any_cast<float>(get_property(*(m_voxel_grid.value()), "px")->property_value());
			const float nx = boost::any_cast<float>(get_property(*(m_voxel_grid.value()), "nx")->property_value());
			const float py = boost::any_cast<float>(get_property(*(m_voxel_grid.value()), "py")->property_value());
			const float ny = boost::any_cast<float>(get_property(*(m_voxel_grid.value()), "ny")->property_value());
			const float pz = boost::any_cast<float>(get_property(*(m_voxel_grid.value()), "pz")->property_value());
			const float nz = boost::any_cast<float>(get_property(*(m_voxel_grid.value()), "nz")->property_value());
			const float vw = boost::any_cast<float>(get_property(*(m_voxel_grid.value()), "voxel_width")->property_value());

			int rows = m_rows;
			int columns = m_columns;
			int slices = m_slices;

			const unsigned int point_rows = rows + 1;
			const unsigned int point_columns = columns + 1;
			const unsigned int point_slices = slices + 1;

			
			k3d::vector3 x, y, z, vox_width, vox_inc, vox_width_x, vox_width_y, vox_width_z;

			x = k3d::vector3(nx, 0, 0);
			y = k3d::vector3(0, py, 0);
			z = k3d::vector3(0, 0, pz);
			vox_width_x = k3d::vector3(vw, 0, 0);
			vox_width_y = k3d::vector3(0, -vw, 0);
			vox_width_z = k3d::vector3(0, 0, -vw);

			vox_inc = k3d::vector3(0,0,0);

		
			k3d::mesh::points_t::iterator point = const_cast<k3d::mesh::points_t&>(*Mesh.points).begin();
        	        for(unsigned long slice = 0; slice != point_slices; ++slice)
                	{
	                        for(unsigned long row = 0; row != point_rows; ++row)
        	                {
					for (unsigned long column = 0; column != point_columns; ++column)
					{
		                                *point = k3d::to_point(x + z + y + vox_inc);
						vox_inc = vox_inc + vox_width_x; 
						++point;
					}
					vox_inc[0] = 0;
					vox_inc = vox_inc + vox_width_z;
                	        }
				vox_inc[2] = 0;
				vox_inc = vox_inc + vox_width_y;
	                }

		}
	}

	void reset_mesh(iunknown* const Hint)
        {

		if (m_voxel_grid.value() != 0) {
			dynamic_cast<voxel_grid*>(m_voxel_grid.value())->voxel_grid_changed_signal.connect(make_topology_changed_slot());
		}
        }



protected:
	 k3d_data(voxel_grid*, immutable_name, change_signal, no_undo, node_storage, no_constraint, node_property, no_serialization) m_voxel_grid;
	
	 float m_width;
	 float m_height;
	 float m_length;
	 float m_voxel_length;

	 int m_rows;
	 int m_columns;
	 int m_slices;


};

k3d::iplugin_factory& voxel_grid_visual_factory()
{
	return voxel_grid_visual::get_factory();
}

}






