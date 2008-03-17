// K-3D
// Copyright (c) 1995-2008, Timothy M. Shead
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

#include "idocument_python.h"
#include "node_python.h"

#include <k3dsdk/classes.h>
#include <k3dsdk/command_node.h>
#include <k3dsdk/plugins.h>
#include <k3dsdk/plugins.h>
#include <k3dsdk/idocument.h>
#include <k3dsdk/idocument_exporter.h>
#include <k3dsdk/ipipeline.h>
#include <k3dsdk/iplugin_factory_collection.h>
#include <k3dsdk/nodes.h>
#include <k3dsdk/state_change_set.h>
#include <k3dsdk/utility_gl.h>

#include <boost/python.hpp>
using namespace boost::python;

#include <boost/scoped_ptr.hpp>

namespace k3d
{

namespace python
{

idocument::idocument() :
	base()
{
}

idocument::idocument(k3d::idocument* Document) :
	base(Document)
{
}

idocument::idocument(k3d::idocument& Document) :
	base(Document)
{
}

const bool idocument::save(const std::string& Path)
{
	boost::scoped_ptr<k3d::idocument_exporter> exporter(k3d::plugin::create<k3d::idocument_exporter>(k3d::classes::DocumentExporter()));
	if(!exporter)
		throw std::runtime_error("no exporter plugin available");

	return exporter->write_file(wrapped(), filesystem::native_path(ustring::from_utf8(Path)));
}

void idocument::start_change_set()
{
	k3d::start_state_change_set(wrapped(), K3D_CHANGE_SET_CONTEXT);
}

void idocument::cancel_change_set()
{
	k3d::cancel_state_change_set(wrapped(), K3D_CHANGE_SET_CONTEXT);
}

void idocument::finish_change_set(const std::string& Label)
{
	k3d::finish_state_change_set(wrapped(), Label, K3D_CHANGE_SET_CONTEXT);
}

void idocument::redraw_all()
{
	k3d::gl::redraw_all(wrapped(), k3d::gl::irender_viewport::ASYNCHRONOUS);
}

const list idocument::nodes()
{
	list results;

	const k3d::inode_collection::nodes_t nodes = wrapped().nodes().collection();
	for(k3d::inode_collection::nodes_t::const_iterator n = nodes.begin(); n != nodes.end(); ++n)
		results.append(node(*n));

	return results;
}

const object idocument::new_node(const object& Type)
{
	extract<std::string> plugin_name(Type);
	if(plugin_name.check())
	{
		k3d::iplugin_factory* const plugin_factory = k3d::plugin::factory::lookup(plugin_name());
		if(!plugin_factory)
			throw std::runtime_error("no factory for plugin type " + plugin_name());

		return object(node(k3d::plugin::create<k3d::iunknown>(*plugin_factory, wrapped(), k3d::unique_name(wrapped().nodes(), plugin_name()))));
	}

	extract<iplugin_factory> plugin_factory(Type);
	if(plugin_factory.check())
	{
		return object(node(k3d::plugin::create<k3d::iunknown>(plugin_factory().wrapped(), wrapped())));
	}

	throw std::invalid_argument("can't create new node from given argument");
}

const object idocument::get_node(const std::string& Name)
{
	return object(node(k3d::find_node(wrapped().nodes(), Name)));
}

void idocument::delete_node(object& Node)
{
	extract<node> node(Node);
	if(!node.check())
		throw std::invalid_argument("argument isn't a node");

	k3d::delete_nodes(wrapped(), k3d::make_collection<k3d::nodes_t>(node().interface_wrapper<k3d::inode>::wrapped_ptr()));
}

object idocument::get_dependency(iproperty& Property)
{
	k3d::iproperty* const property = Property.wrapped_ptr();
	if(!property)
		throw std::invalid_argument("property cannot be null");

	k3d::iproperty* const dependency = wrapped().pipeline().dependency(*property);
	return dependency ? object(iproperty(dependency)) : object();
}

void idocument::set_dependency(iproperty& From, boost::python::object& To)
{
	k3d::iproperty* to = 0;

	extract<iproperty> iproperty_value(To);
	if(iproperty_value.check())
	{
		to = iproperty_value().wrapped_ptr();
	}
	else if(To.ptr() == boost::python::object().ptr())
	{
		to = 0;
	}
	else
	{
		throw std::invalid_argument("to property must be an iproperty instance or None");
	}

	k3d::iproperty* const from = From.wrapped_ptr();
	if(!from)
		throw std::invalid_argument("from property cannot be null");

	if(from && to && from->property_type() != to->property_type())
		throw std::invalid_argument("property types do not match");

	k3d::ipipeline::dependencies_t dependencies;
	dependencies[from] = to;
	wrapped().pipeline().set_dependencies(dependencies);
}

void idocument::define_class()
{
	class_<idocument>("idocument")
		.def("save", &idocument::save)
		.def("start_change_set", &idocument::start_change_set)
		.def("cancel_change_set", &idocument::cancel_change_set)
		.def("finish_change_set", &idocument::finish_change_set)
		.def("redraw_all", &idocument::redraw_all)
		.def("nodes", &idocument::nodes)
		.def("new_node", &idocument::new_node)
		.def("get_node", &idocument::get_node)
		.def("delete_node", &idocument::delete_node)
		.def("get_dependency", &idocument::get_dependency)
		.def("set_dependency", &idocument::set_dependency);
}

} // namespace python

} // namespace k3d

