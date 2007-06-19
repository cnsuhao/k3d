// K-3D
// Copyright (c) 1995-2007, Timothy M. Shead
//
// Contact: tshead@k-3d.com
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

/** \file
		\author Timothy M. Shead (tshead@k-3d.com)
*/

#include <k3dsdk/document_plugin_factory.h>
#include <k3dsdk/i18n.h>
#include <k3dsdk/mesh_painter_gl.h>
#include <k3dsdk/mesh.h>
#include <k3dsdk/painter_render_state_gl.h>
#include <k3dsdk/persistent.h>
#include <k3dsdk/selection.h>

namespace libk3dquadremesh
{

/////////////////////////////////////////////////////////////////////////////
// face_normal_painter

class point_mean_curv_painter :
	public k3d::gl::mesh_painter
{
	typedef k3d::gl::mesh_painter base;

public:
	point_mean_curv_painter(k3d::iplugin_factory& Factory, k3d::idocument& Document) :
		base(Factory, Document),
		m_draw_selected(init_owner(*this) + init_name("draw_selected") + init_label(_("Draw Selected")) + init_description(_("Draw normals for selected polygons")) + init_value(true)),
		m_draw_unselected(init_owner(*this) + init_name("draw_unselected") + init_label(_("Draw Unselected")) + init_description(_("Draw normals for unselected polygons")) + init_value(false)),
		m_selected_color(init_owner(*this) + init_name("selected_color") + init_label(_("Selected Color")) + init_description(_("Normal color for selected polygons")) + init_value(k3d::color(0, 1, 1))),
		m_unselected_color(init_owner(*this) + init_name("unselected_color") + init_label(_("Unselected Color")) + init_description(_("Normal color for unselected polygons")) + init_value(k3d::color(0, 0.6, 0.6))),
		m_scale(init_owner(*this) + init_name("scale") + init_label(_("Scale")) + init_description(_("Scaling of vectors")) + init_value(1.0))
	{
		m_draw_selected.changed_signal().connect(make_async_redraw_slot());
		m_draw_unselected.changed_signal().connect(make_async_redraw_slot());
		m_selected_color.changed_signal().connect(make_async_redraw_slot());
		m_unselected_color.changed_signal().connect(make_async_redraw_slot());
		m_scale.changed_signal().connect(make_async_redraw_slot());
	}

	void on_paint_mesh(const k3d::mesh& Mesh, const k3d::gl::painter_render_state& RenderState)
	{
		if(!k3d::validate_polyhedra(Mesh))
			return;

		if(Mesh.vertex_data.find("PGPMeanCurv") == Mesh.vertex_data.end())
			return;
		if(Mesh.vertex_data.find("PGPMeanCurv") == Mesh.vertex_data.end())
			return;
		if(Mesh.vertex_data.find("PGPMeanCurv") == Mesh.vertex_data.end())
			return;
		
		const bool draw_selected = m_draw_selected.value();
		const bool draw_unselected = m_draw_unselected.value();
		k3d::log() << debug << "PGP MeanCurv" << std::endl;

		const k3d::mesh::points_t& points = *Mesh.points;

		const size_t vert_count = points.size();

		const k3d::typed_array < k3d::vector3 > & norm = dynamic_cast<const k3d::typed_array < k3d::vector3 > & >(*((*(Mesh.vertex_data.find("PGPMeanCurv"))).second)); 
		const k3d::typed_array < k3d::vector3 > & p1 = dynamic_cast<const k3d::typed_array < k3d::vector3 > & >(*((*(Mesh.vertex_data.find("PGPPrincCurv1"))).second)); 
		const k3d::typed_array < k3d::vector3 > & p2 = dynamic_cast<const k3d::typed_array < k3d::vector3 > & >(*((*(Mesh.vertex_data.find("PGPPrincCurv2"))).second)); 
		k3d::gl::store_attributes attributes;
		glDisable(GL_LIGHTING);
		double scale = m_scale.value();
		k3d::gl::color3d(m_selected_color.value());
		k3d::point3 x;
		glBegin(GL_LINES);
		for(size_t vert = 0; vert != vert_count; ++vert)
		{
			//k3d::gl::color3d(m_selected_color.value());
			//k3d::gl::vertex3d(points[vert]);
			//k3d::gl::vertex3d(points[vert] + k3d::to_point(norm[vert]));
			
			if(draw_selected) {
				k3d::gl::color3d(k3d::color(1,0,0));

				//k3d::gl::vertex3d(points[vert]);

				x =  k3d::to_point(p1[vert]);
				x *= scale;
				k3d::gl::vertex3d(points[vert] + x);
				k3d::gl::vertex3d(points[vert] - x);
			}

			if(draw_unselected) {
				k3d::gl::color3d(k3d::color(0,1,0));
				//k3d::gl::vertex3d(points[vert]);
				
				x = k3d::to_point(p2[vert]);
				x *= scale;
				k3d::gl::vertex3d(points[vert] + x);
				k3d::gl::vertex3d(points[vert] - x);
			}
		}
		glEnd();
	}
	
	static k3d::iplugin_factory& get_factory()
	{
		static k3d::document_plugin_factory<point_mean_curv_painter, k3d::interface_list<k3d::gl::imesh_painter > > factory(
			k3d::uuid(0x23b2a6a7, 0x4a4c9f08, 0x8bc5e6a2, 0xcc06655a),
			"OpenGLMeanCurvPainter",
			_("Renders mean curvature normal vector at each vertex"),
			"PGP",
			k3d::iplugin_factory::EXPERIMENTAL);

		return factory;
	}

private:
	k3d_data(bool, immutable_name, change_signal, with_undo, local_storage, no_constraint, writable_property, with_serialization) m_draw_selected;
	k3d_data(bool, immutable_name, change_signal, with_undo, local_storage, no_constraint, writable_property, with_serialization) m_draw_unselected;
	k3d_data(k3d::color, immutable_name, change_signal, with_undo, local_storage, no_constraint, writable_property, with_serialization) m_selected_color;
	k3d_data(k3d::color, immutable_name, change_signal, with_undo, local_storage, no_constraint, writable_property, with_serialization) m_unselected_color;
	k3d_data(double, immutable_name, change_signal, with_undo, local_storage, no_constraint, writable_property, with_serialization) m_scale;
};

/////////////////////////////////////////////////////////////////////////////
// face_normal_painter_factory

k3d::iplugin_factory& point_mean_curv_painter_factory()
{
	return point_mean_curv_painter::get_factory();
}

} // namespace libk3dquadremesh

