// K-3D
// Copyright (c) 1995-2006, Timothy M. Shead
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

#include "angle_axis_control.h"
#include "aqsis_layer_chooser.h"
#include "asynchronous_update.h"
#include "bitmap_preview.h"
#include "bounding_box.h"
#include "button.h"
#include "check_button.h"
#include "collapsible_frame.h"
#include "color_chooser.h"
#include "combo_box.h"
#include "document_state.h"
#include "entry.h"
#include "enumeration_chooser.h"
#include "file_chooser_dialog.h"
#include "icons.h"
#include "messages.h"
#include "node_chooser.h"
#include "node_properties.h"
#include "open_uri.h"
#include "path_chooser.h"
#include "point_control.h"
#include "property_button.h"
#include "property_label.h"
#include "render.h"
#include "scale.h"
#include "script_button.h"
#include "selection_button.h"
#include "spin_button.h"
#include "toolbar.h"
#include "ui_component.h"
#include "user_property.h"
#include "utility.h"
#include "widget_manip.h"

#include <k3dsdk/i18n.h>
#include <k3dsdk/ianimation_render_engine.h>
#include <k3dsdk/iaqsis.h>
#include <k3dsdk/icamera.h>
#include <k3dsdk/icamera_animation_render_engine.h>
#include <k3dsdk/icamera_preview_render_engine.h>
#include <k3dsdk/icamera_still_render_engine.h>
#include <k3dsdk/ideletable.h>
#include <k3dsdk/idocument.h>
#include <k3dsdk/ienumeration_property.h>
#include <k3dsdk/ilist_property.h>
#include <k3dsdk/imeasurement_property.h>
#include <k3dsdk/imesh_storage.h>
#include <k3dsdk/inode.h>
#include <k3dsdk/iplugin_factory.h>
#include <k3dsdk/ipreview_render_engine.h>
#include <k3dsdk/iproperty_group_collection.h>
#include <k3dsdk/iscript_property.h>
#include <k3dsdk/iselectable.h>
#include <k3dsdk/istill_render_engine.h>
#include <k3dsdk/iuser_property.h>
#include <k3dsdk/mesh_selection.h>
#include <k3dsdk/mesh.h>
#include <k3dsdk/options.h>
#include <k3dsdk/property.h>
#include <k3dsdk/types_ri.h>
#include <k3dsdk/state_change_set.h>
#include <k3dsdk/string_cast.h>
#include <k3dsdk/string_modifiers.h>
#include <k3dsdk/system.h>
#include <k3dsdk/types.h>
#include <k3dsdk/user_properties.h>
#include <k3dsdk/utility.h>

// Not strictly required to compile, but this #include ensures that we have a std::typeinfo for k3d::legacy::mesh that matches the SDK (i.e. we don't break the ODR)
#include <k3dsdk/legacy_mesh.h>

#include <gtkmm/arrow.h>
#include <gtkmm/box.h>
#include <gtkmm/frame.h>
#include <gtkmm/image.h>
#include <gtkmm/label.h>
#include <gtkmm/scrolledwindow.h>
#include <gtkmm/stock.h>
#include <gtkmm/table.h>

#include <k3dsdk/fstream.h>
#include <boost/format.hpp>

namespace libk3dngui
{

namespace node_properties
{

/////////////////////////////////////////////////////////////////////////////
// control::implementation

class control::implementation :
	public asynchronous_update
{
public:
	implementation(document_state& DocumentState, k3d::icommand_node& Parent) :
		m_document_state(DocumentState),
		m_node(0),
		m_parent(Parent),
		m_help_button(m_parent, "onlin_help", Gtk::Stock::HELP)
	{
		m_label.set_alignment(Gtk::ALIGN_LEFT);
		m_label.set_padding(5, 5);

		&m_help_button << connect_button(sigc::mem_fun(*this, &implementation::on_online_help));

		m_scrolled_window.set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);
		m_scrolled_window.add(m_vbox);

		// If only one node is selected, show its properties
		const k3d::nodes_t nodes = m_document_state.selected_nodes();
		if(1 == nodes.size())
			on_view_node_properties(nodes.front());

		m_document_state.document().close_signal().connect(sigc::mem_fun(*this, &implementation::on_document_closed));
		m_document_state.view_node_properties_signal().connect(sigc::mem_fun(*this, &implementation::on_view_node_properties));

		schedule_update();
	}

	void on_document_closed()
	{
		block_updates();
	}

	bool on_view_node_properties(k3d::inode* const Node)
	{
		if(Node != m_node)
		{
			m_node = Node;

			m_node_deleted_connection.disconnect();
			m_node_name_change_connection.disconnect();
			m_node_properties_changed_connection.disconnect();

			schedule_update();

			if(m_node)
			{
				m_node_deleted_connection = m_node->deleted_signal().connect(sigc::mem_fun(*this, &implementation::on_node_deleted));
				m_node_name_change_connection = m_node->name_changed_signal().connect(sigc::mem_fun(*this, &implementation::update_label));
				k3d::iproperty_collection* const property_collection = dynamic_cast<k3d::iproperty_collection*>(m_node);
				if(property_collection)
					m_node_properties_changed_connection = property_collection->connect_properties_changed_signal(sigc::hide(sigc::mem_fun(*this, &implementation::on_node_properties_changed)));

			}

			return true;
		}

		return false;
	}

	void on_node_properties_changed()
	{
		m_vbox.hide();
		schedule_update();
	}

	void on_node_deleted()
	{
		on_view_node_properties(0);
	}

	void update_label()
	{
		if(m_node)
		{
			m_label.set_text(m_node->name());
			m_help_button.set_sensitive(true);
		}
		else
		{
			m_label.set_text("");
			m_help_button.set_sensitive(false);
		}
	}

	void reset()
	{
		Glib::ListHandle<Widget*> children = m_vbox.get_children();
		std::for_each(children.begin(), children.end(), k3d::delete_object());
	}

	void on_update()
	{
		update_label();

		reset();

		// Create a toolbar ...
		toolbar::control* const toolbar_control = new toolbar::control(m_parent, "toolbar");
		m_vbox.pack_start(*manage(toolbar_control), Gtk::PACK_SHRINK);

		k3d::istate_recorder* const state_recorder = &m_document_state.document().state_recorder();

		// Add controls for cameras and camera render engines ...
		if(dynamic_cast<k3d::icamera*>(m_node) || dynamic_cast<k3d::icamera_preview_render_engine*>(m_node))
		{
			button::control* const control =
				new button::control(m_parent, "render_camera_preview", *Gtk::manage(new Gtk::Image(load_icon("render_preview", Gtk::ICON_SIZE_BUTTON))))
					<< connect_button(sigc::mem_fun(*this, &implementation::on_render_camera_preview))
					<< set_tooltip(_("Render Preview"));

			toolbar_control->row(0).pack_start(*Gtk::manage(control), Gtk::PACK_SHRINK);
		}

		if(dynamic_cast<k3d::icamera*>(m_node) || dynamic_cast<k3d::icamera_still_render_engine*>(m_node))
		{
			button::control* const control =
				new button::control(m_parent, "render_camera_frame", *Gtk::manage(new Gtk::Image(load_icon("render_frame", Gtk::ICON_SIZE_BUTTON))))
					<< connect_button(sigc::mem_fun(*this, &implementation::on_render_camera_frame))
					<< set_tooltip(_("Render Frame"));

			toolbar_control->row(0).pack_start(*Gtk::manage(control), Gtk::PACK_SHRINK);
		}

		if(dynamic_cast<k3d::icamera*>(m_node) || dynamic_cast<k3d::icamera_animation_render_engine*>(m_node))
		{
			button::control* const control =
				new button::control(m_parent, "render_camera_animation", *Gtk::manage(new Gtk::Image(load_icon("render_animation", Gtk::ICON_SIZE_BUTTON))))
					<< connect_button(sigc::mem_fun(*this, &implementation::on_render_camera_animation))
					<< set_tooltip(_("Render Animation"));

			toolbar_control->row(0).pack_start(*Gtk::manage(control), Gtk::PACK_SHRINK);
		}

		// Add controls for render engines
		if(dynamic_cast<k3d::ipreview_render_engine*>(m_node))
		{
			button::control* const control =
				new button::control(m_parent, "render_preview", *Gtk::manage(new Gtk::Image(load_icon("render_preview", Gtk::ICON_SIZE_BUTTON))))
					<< connect_button(sigc::mem_fun(*this, &implementation::on_render_preview))
					<< set_tooltip(_("Render Preview"));

			toolbar_control->row(0).pack_start(*Gtk::manage(control), Gtk::PACK_SHRINK);
		}

		if(dynamic_cast<k3d::istill_render_engine*>(m_node))
		{
			button::control* const control =
				new button::control(m_parent, "render_frame", *Gtk::manage(new Gtk::Image(load_icon("render_frame", Gtk::ICON_SIZE_BUTTON))))
					<< connect_button(sigc::mem_fun(*this, &implementation::on_render_frame))
					<< set_tooltip(_("Render Frame"));

			toolbar_control->row(0).pack_start(*Gtk::manage(control), Gtk::PACK_SHRINK);
		}

		if(dynamic_cast<k3d::ianimation_render_engine*>(m_node))
		{
			button::control* const control =
				new button::control(m_parent, "render_animation", *Gtk::manage(new Gtk::Image(load_icon("render_animation", Gtk::ICON_SIZE_BUTTON))))
					<< connect_button(sigc::mem_fun(*this, &implementation::on_render_animation))
					<< set_tooltip(_("Render Animation"));

			toolbar_control->row(0).pack_start(*Gtk::manage(control), Gtk::PACK_SHRINK);
		}

		// Add a "reset" button for nodes that implement k3d::imesh_storage (FrozenMesh, external file readers, etc) ...
		if(dynamic_cast<k3d::imesh_storage*>(m_node))
		{
			button::control* const control =
				new button::control(m_parent, "reset_mesh", _("Reset Mesh"))
					<< connect_button(sigc::mem_fun(*this, &implementation::on_reset_mesh))
					<< set_tooltip(_("Reset / Reload Mesh"));

			toolbar_control->row(0).pack_start(*Gtk::manage(control), Gtk::PACK_SHRINK);
		}

		// Get the node properties, grouped together ...
		k3d::iproperty_collection* const property_collection = dynamic_cast<k3d::iproperty_collection*>(m_node);
		k3d::iproperty_group_collection::groups_t property_groups;
		if(property_collection)
		{
			k3d::iproperty_collection::properties_t all_properties = property_collection->properties();
			k3d::iproperty_group_collection::groups_t groups;

			k3d::iproperty_group_collection* const property_group_collection = dynamic_cast<k3d::iproperty_group_collection*>(m_node);
			if(property_group_collection)
			{
				groups = property_group_collection->property_groups();
				for(k3d::iproperty_group_collection::groups_t::const_iterator group = groups.begin(); group != groups.end(); ++group)
				{
					for(k3d::iproperty_collection::properties_t::const_iterator property = group->properties.begin(); property != group->properties.end(); ++property)
						all_properties.erase(std::remove(all_properties.begin(), all_properties.end(), *property), all_properties.end());
				}
			}

			property_groups.insert(property_groups.end(), k3d::iproperty_group_collection::group(m_node->factory().name(), all_properties));
			property_groups.insert(property_groups.end(), groups.begin(), groups.end());
		}

		// For each property group ...
		for(k3d::iproperty_group_collection::groups_t::const_iterator property_group = property_groups.begin(); property_group != property_groups.end(); ++property_group)
		{
			if(property_group->properties.empty())
				continue;

			collapsible_frame::control* const frame = new collapsible_frame::control(property_group->name, m_collapsible_frame_group);
			m_vbox.pack_start(*manage(frame), Gtk::PACK_SHRINK);

			Gtk::Table* const table = new Gtk::Table(property_group->properties.size(), 5, false);
			frame->add(*manage(table));

			// Store entries for focus chain within table
			std::list<Gtk::Widget*> entry_list;

			const unsigned long prop_delete_begin = 0;
			const unsigned long prop_delete_end = 1;
			const unsigned long prop_button_begin = 1;
			const unsigned long prop_button_end = 2;
			const unsigned long prop_label_begin = 2;
			const unsigned long prop_label_end = 3;
			const unsigned long prop_control_begin = 3;
			const unsigned long prop_control_end = 4;

			// For each property within the group ...
			unsigned int row = 0;
			for(unsigned int i = 0; i != property_group->properties.size(); ++i, ++row)
			{
				k3d::iproperty& property = *property_group->properties[i];

				const std::string property_name = property.property_name();
				const std::type_info& property_type = property.property_type();

				// Provide a property button for the property ...
				table->attach(*manage(
					new property_button::control(m_parent, property_name + "_property", property_widget::proxy(m_document_state,property))),
					prop_button_begin, prop_button_end, row, row + 1, Gtk::SHRINK, Gtk::SHRINK);

				// Provide a label for the property ...
				table->attach(*manage(
					new property_label::control(m_parent, property_name + "_label", property_widget::proxy(m_document_state, property))),
					prop_label_begin, prop_label_end, row, row + 1, Gtk::FILL | Gtk::SHRINK, Gtk::FILL | Gtk::SHRINK);

				// Boolean properties ...
				if(property_type == typeid(bool))
				{
					check_button::control* const control = new check_button::control(m_parent, property_name, check_button::proxy(property, state_recorder, property_name));
					table->attach(*manage(control), prop_control_begin, prop_control_end, row, row + 1, Gtk::FILL | Gtk::SHRINK, Gtk::FILL | Gtk::SHRINK);

					entry_list.push_back(control);
				}
				// Scalar properties ...
				else if(property_type == typeid(double) || property_type == typeid(float) || property_type == typeid(long) || property_type == typeid(unsigned long) || property_type == typeid(int) || property_type == typeid(unsigned int))
				{
					spin_button::control* const control = new spin_button::control(m_parent, property_name, spin_button::proxy(property, state_recorder, property_name));
					if(k3d::imeasurement_property* const measurement_property = dynamic_cast<k3d::imeasurement_property*>(&property))
					{
						control->set_step_increment(measurement_property->property_step_increment());
						control->set_units(measurement_property->property_units());
					}
					table->attach(*manage(control), prop_control_begin, prop_control_end, row, row + 1, Gtk::FILL | Gtk::SHRINK, Gtk::FILL | Gtk::SHRINK);

					entry_list.push_back(control);
				}
				// Color properties ...
				else if(property_type == typeid(k3d::color))
				{
					color_chooser::control* const control = new color_chooser::control(m_parent, property_name, color_chooser::proxy(property, state_recorder, property_name));
					table->attach(*manage(control), prop_control_begin, prop_control_end, row, row + 1, Gtk::FILL | Gtk::SHRINK, Gtk::FILL | Gtk::SHRINK);

					entry_list.push_back(control);
				}
				// String properties ...
				else if(property_type == typeid(std::string))
				{
					if(dynamic_cast<k3d::ienumeration_property*>(&property))
					{
						enumeration_chooser::control* const control = new enumeration_chooser::control(m_parent, property_name, enumeration_chooser::proxy(property, state_recorder, property_name));
						table->attach(*manage(control), prop_control_begin, prop_control_end, row, row + 1, Gtk::FILL | Gtk::SHRINK, Gtk::FILL | Gtk::SHRINK);

						entry_list.push_back(control);
					}
					else if(dynamic_cast<k3d::iscript_property*>(&property))
					{
						script_button::control* const control = new script_button::control(m_parent, property_name, script_button::proxy(property, state_recorder, property_name));
						table->attach(*Gtk::manage(control), prop_control_begin, prop_control_end, row, row + 1, Gtk::FILL | Gtk::SHRINK, Gtk::FILL | Gtk::SHRINK);

						entry_list.push_back(control);
					}
					else if(k3d::ilist_property<std::string>* const list_property = dynamic_cast<k3d::ilist_property<std::string>*>(&property))
					{
						combo_box::control* const control = new combo_box::control(m_parent, property_name, combo_box::proxy(property, state_recorder, property_name));
						control->set_values(list_property->property_values());
						table->attach(*manage(control), prop_control_begin, prop_control_end, row, row + 1, Gtk::FILL | Gtk::SHRINK, Gtk::FILL | Gtk::SHRINK);

						entry_list.push_back(control);
					}
					else
					{
						entry::control* const control = new entry::control(m_parent, property_name, entry::proxy(property, state_recorder, property_name));
						table->attach(*manage(control), prop_control_begin, prop_control_end, row, row + 1, Gtk::FILL | Gtk::SHRINK, Gtk::FILL | Gtk::SHRINK);

						entry_list.push_back(control);
					}
				}
				// k3d::aqsis::ilayer_connection* properties
				else if(k3d::aqsis::ilayer_connection_property* const layer_connection_property = dynamic_cast<k3d::aqsis::ilayer_connection_property*>(&property))
				{
					aqsis_layer_chooser::control* const control = new aqsis_layer_chooser::control(m_document_state, *layer_connection_property, m_parent, property_name, state_recorder);
					table->attach(*manage(control), prop_control_begin, prop_control_end, row, row + 1, Gtk::FILL | Gtk::SHRINK, Gtk::FILL | Gtk::SHRINK);

					entry_list.push_back(control);
				}
				// k3d::inode* properties ...
				else if(property_type == typeid(k3d::inode*))
				{
					node_chooser::control* const control = new node_chooser::control(m_parent, property_name, node_chooser::proxy(m_document_state, property, state_recorder, property_name), node_chooser::filter(property));
					table->attach(*manage(control), prop_control_begin, prop_control_end, row, row + 1, Gtk::FILL | Gtk::SHRINK, Gtk::FILL | Gtk::SHRINK);

					entry_list.push_back(control);
				}
				// Bitmap properties ...
				else if(property_type == k3d::type_id_k3d_bitmap_ptr())
				{
					bitmap_preview::control* const control = new bitmap_preview::control(m_parent, property_name, bitmap_preview::proxy(m_document_state.document().dag(), property));
					table->attach(*manage(control), prop_control_begin, prop_control_end, row, row + 1, Gtk::FILL | Gtk::SHRINK, Gtk::FILL | Gtk::SHRINK);

					entry_list.push_back(control);
				}
				// Filesystem-path properties ...
				else if(property_type == typeid(k3d::filesystem::path))
				{
					path_chooser::control* const control = new path_chooser::control(m_parent, property_name, path_chooser::proxy(property, state_recorder, property_name));
					table->attach(*manage(control), prop_control_begin, prop_control_end, row, row + 1, Gtk::FILL | Gtk::SHRINK, Gtk::FILL | Gtk::SHRINK);

					entry_list.push_back(control);
				}
				// k3d::bounding_box3 properties ...
				else if(property_type == typeid(k3d::bounding_box3))
				{
					bounding_box::control* const control = new bounding_box::control(m_parent, property_name, bounding_box::proxy(property));
					table->attach(*manage(control), prop_control_begin, prop_control_end, row, row + 1, Gtk::FILL | Gtk::SHRINK, Gtk::FILL | Gtk::SHRINK);

					entry_list.push_back(control);
				}
				// k3d::point3 properties ...
				else if(property_type == typeid(k3d::point3) || property_type == typeid(k3d::vector3) || property_type == typeid(k3d::normal3))
				{
					point::control* const control = new point::control(m_parent, property_name, point::proxy(property));
					table->attach(*manage(control), prop_control_begin, prop_control_end, row, row + 1, Gtk::FILL | Gtk::SHRINK, Gtk::FILL | Gtk::SHRINK);

					entry_list.push_back(control);
				}
				// k3d::angle_axis properties ...
				else if(property_type == typeid(k3d::angle_axis))
				{
					angle_axis::control* const control = new angle_axis::control(m_parent, property_name, angle_axis::proxy(property, state_recorder, property_name));
					table->attach(*manage(control), prop_control_begin, prop_control_end, row, row + 1, Gtk::FILL | Gtk::SHRINK, Gtk::FILL | Gtk::SHRINK);

					entry_list.push_back(control);
				}
				// Transformation properties ...
				else if(property_type == typeid(k3d::matrix4))
				{
				}
				// Mesh properties ...
				else if(property_type == typeid(k3d::legacy::mesh*))
				{
				}
				else if(property_type == typeid(k3d::mesh*))
				{
				}
				// HPoint properties ...
				else if(property_type == typeid(k3d::point4))
				{
				}
				// Mesh Selection properties ...
				else if(property_type == typeid(k3d::mesh_selection))
				{
					selection_button::control* const control = new selection_button::control(m_parent, property_name, selection_button::proxy(property, state_recorder, property_name));
					table->attach(*manage(control), prop_control_begin, prop_control_end, row, row + 1, Gtk::FILL | Gtk::SHRINK, Gtk::FILL | Gtk::SHRINK);
				}
				else
				{
					k3d::log() << warning << k3d_file_reference << "unknown property type: " << property_type.name() << " name: " << property_name << std::endl;
				}

				// Provide a "delete" button for user properties ...
				if(dynamic_cast<k3d::iuser_property*>(&property))
				{
					button::control* const control =
						new button::control(m_parent, property_name + "_delete", *Gtk::manage(new Gtk::Image(Gtk::Stock::DELETE, Gtk::ICON_SIZE_BUTTON)))
						<< connect_button(sigc::bind(sigc::bind(sigc::mem_fun(*this, &implementation::on_delete_user_property), &property), property_collection))
						<< set_tooltip(_("Delete user property (no undo)"));

					table->attach(*manage(control), prop_delete_begin, prop_delete_end, row, row + 1, Gtk::SHRINK, Gtk::SHRINK);
				}
			}

			// Add controls for managing user properties ...
			if(property_collection)
			{
				button::control* const control =
					new button::control(m_parent, "add_user_property", *Gtk::manage(new Gtk::Image(Gtk::Stock::ADD, Gtk::ICON_SIZE_BUTTON)))
						<< connect_button(sigc::mem_fun(*this, &implementation::on_add_user_property))
						<< set_tooltip(_("Add a user property to this node"));

				table->attach(*Gtk::manage(control), prop_delete_begin, prop_delete_end, row, row + 1, Gtk::SHRINK, Gtk::SHRINK);
			}

			// Set focus chain
			table->set_focus_chain(entry_list);
		}

		m_vbox.show_all();
	}

	void on_online_help()
	{
		if(m_node)
		{
			open_uri("http://www.k-3d.org/wiki/" + m_node->factory().name());
		}
	}

	void on_render_camera_preview()
	{
		k3d::icamera* camera = dynamic_cast<k3d::icamera*>(m_node);
		if(!camera)
			camera = pick_camera(m_document_state);
		if(!camera)
			return;

		k3d::icamera_preview_render_engine* render_engine = dynamic_cast<k3d::icamera_preview_render_engine*>(m_node);
		if(!render_engine)
			render_engine = pick_camera_preview_render_engine(m_document_state);
		if(!render_engine)
			return;

		render_camera_preview(*camera, *render_engine);
	}

	void on_render_camera_frame()
	{
		k3d::icamera* camera = dynamic_cast<k3d::icamera*>(m_node);
		if(!camera)
			camera = pick_camera(m_document_state);
		if(!camera)
			return;

		k3d::icamera_still_render_engine* render_engine = dynamic_cast<k3d::icamera_still_render_engine*>(m_node);
		if(!render_engine)
			render_engine = pick_camera_still_render_engine(m_document_state);
		if(!render_engine)
			return;

		render_camera_frame(*camera, *render_engine);
	}

	void on_render_camera_animation()
	{
		k3d::icamera* camera = dynamic_cast<k3d::icamera*>(m_node);
		if(!camera)
			camera = pick_camera(m_document_state);
		if(!camera)
			return;

		k3d::icamera_animation_render_engine* render_engine = dynamic_cast<k3d::icamera_animation_render_engine*>(m_node);
		if(!render_engine)
			render_engine = pick_camera_animation_render_engine(m_document_state);
		if(!render_engine)
			return;

		render_camera_animation(m_document_state, *camera, *render_engine);
	}

	void on_render_preview()
	{
		k3d::ipreview_render_engine* render_engine = dynamic_cast<k3d::ipreview_render_engine*>(m_node);
		return_if_fail(render_engine);

		render_preview(*render_engine);
	}

	void on_render_frame()
	{
		k3d::istill_render_engine* render_engine = dynamic_cast<k3d::istill_render_engine*>(m_node);
		return_if_fail(render_engine);

		render_frame(*render_engine);
	}

	void on_render_animation()
	{
		k3d::ianimation_render_engine* render_engine = dynamic_cast<k3d::ianimation_render_engine*>(m_node);
		return_if_fail(render_engine);

		render_animation(m_document_state, *render_engine);
	}

    void on_reset_mesh()
    {
            k3d::imesh_storage* const mesh_storage = dynamic_cast<k3d::imesh_storage*>(m_node);
            return_if_fail(mesh_storage);

            mesh_storage->reset_mesh(0);
    }

	void on_delete_user_property(k3d::iproperty_collection* Collection, k3d::iproperty* Property)
	{
		return_if_fail(Collection);
		return_if_fail(Property);
		return_if_fail(dynamic_cast<k3d::iuser_property*>(Property));

		k3d::record_state_change_set change_set(m_document_state.document(), "Delete user property", __PRETTY_FUNCTION__);

		if(m_document_state.document().state_recorder().current_change_set())
			m_document_state.document().state_recorder().current_change_set()->record_old_state(new k3d::user::property_container(*Collection));

		Collection->unregister_property(*Property);
		if(k3d::ipersistent* const persistent = dynamic_cast<k3d::ipersistent*>(Property))
		{
			if(k3d::ipersistent_container* const persistent_container = dynamic_cast<k3d::ipersistent_container*>(Collection))
				persistent_container->disable_serialization(*persistent);
		}

		if(k3d::ideletable* const deletable = dynamic_cast<k3d::ideletable*>(Property))
			undoable_delete(deletable, m_document_state.document());

		if(m_document_state.document().state_recorder().current_change_set())
			m_document_state.document().state_recorder().current_change_set()->record_new_state(new k3d::user::property_container(*Collection));
	}

	void on_add_user_property()
	{
		return_if_fail(m_node);
		add_user_property* const window = new add_user_property(*m_node, m_parent);
		window->show();
	}

	/// Stores a reference to the owning document
	document_state& m_document_state;
	/// Stores a reference to the currently-selected node (if any)
	k3d::inode* m_node;
	k3d::icommand_node& m_parent;
	/// Tracks whether the currently-visible node is deleted
	sigc::connection m_node_deleted_connection;
	/// Keeps track of node name changes
	sigc::connection m_node_name_change_connection;
	/// Keeps track of changes to the set of node properties
	sigc::connection m_node_properties_changed_connection;
	/// Displays the current node name
	Gtk::Label m_label;
	/// Online help button
	button::control m_help_button;
	/// Contains the set of node properties
	Gtk::ScrolledWindow m_scrolled_window;
	/// Parent widget for the rest of the implementation
	Gtk::VBox m_vbox;
	/// Groups collapsible frames together
	collapsible_frame::group m_collapsible_frame_group;

	sigc::signal<void, const std::string&, const std::string&> m_command_signal;

	/// Signal that will be emitted whenever this control should grab the panel focus
	sigc::signal<void> m_panel_grab_signal;
};

/////////////////////////////////////////////////////////////////////////////
// control

control::control(document_state& DocumentState, k3d::icommand_node& Parent) :
	ui_component("node_properties", &Parent),
	m_implementation(new implementation(DocumentState, *this))
{
	m_implementation->m_command_signal.connect(sigc::mem_fun(*this, &control::record_command));

	m_implementation->m_scrolled_window.signal_button_press_event().connect(sigc::bind_return(sigc::hide(m_implementation->m_panel_grab_signal.make_slot()), false), false);

	Gtk::HBox* const hbox = new Gtk::HBox();
	hbox->pack_start(m_implementation->m_label, Gtk::PACK_EXPAND_WIDGET);
	hbox->pack_start(m_implementation->m_help_button, Gtk::PACK_SHRINK);

	pack_start(*Gtk::manage(hbox), Gtk::PACK_SHRINK);
	pack_start(m_implementation->m_scrolled_window, Gtk::PACK_EXPAND_WIDGET);
	show_all();
}

control::~control()
{
	delete m_implementation;
}

sigc::connection control::connect_focus_signal(const sigc::slot<void>& Slot)
{
	return m_implementation->m_panel_grab_signal.connect(Slot);
}

} // namespace node_properties

} // namespace libk3dngui
