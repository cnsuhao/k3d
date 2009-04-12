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
		\author Tim Shead <tshead@k-3d.com>
*/

#include <k3d-i18n-config.h>
#include <k3dsdk/classes.h>
#include <k3dsdk/document_plugin_factory.h>
#include <k3dsdk/imaterial.h>
#include <k3dsdk/node.h>

namespace module
{

namespace core
{

/////////////////////////////////////////////////////////////////////////////
// material

class multi_material :
	public k3d::node,
	public k3d::imaterial
{
	typedef k3d::node base;

public:
	multi_material(k3d::iplugin_factory& Factory, k3d::idocument& Document) :
		base(Factory, Document)
	{
	}

	static k3d::iplugin_factory& get_factory()
	{
		static k3d::document_plugin_factory<multi_material,
			k3d::interface_list<k3d::imaterial> > factory(
				k3d::classes::MultiMaterial(),
				"MultiMaterial",
				_("Material"),
				"Material",
				k3d::iplugin_factory::STABLE);

		return factory;
	}
};

k3d::iplugin_factory& multi_material_factory()
{
	return multi_material::get_factory();
}

} // namespace core

} // namespace module

