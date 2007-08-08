// K-3D
// Copyright (c) 1995-2005, Timothy M. Shead
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
	\author Anders Dahnielson (anders@dahnielson.com)
*/

#include "simple_bitmap_modifier.h"

#include <k3dsdk/document_plugin_factory.h>
#include <k3d-i18n-config.h>

namespace libk3dbitmap
{

/////////////////////////////////////////////////////////////////////////////
// bitmap_matte_colordiff

class bitmap_matte_colordiff :
	public simple_bitmap_modifier
{
	typedef simple_bitmap_modifier base;

public:
	bitmap_matte_colordiff(k3d::iplugin_factory& Factory, k3d::idocument& Document) :
		base(Factory, Document),
		m_value(init_owner(*this) + init_name("value") + init_label(_("Threshold Value")) + init_description(_("Set alpha channel using color difference")) + init_value(0.0))
	{
		m_value.changed_signal().connect(make_update_bitmap_slot());
	}

	struct functor
	{
		functor(const double Threshold) :
			threshold(Threshold)
		{
		}

		k3d::pixel operator()(const k3d::pixel& Input)
		{
			const half red = boost::gil::get_color(Input, boost::gil::red_t());
			const half green = boost::gil::get_color(Input, boost::gil::green_t());
			const half blue = boost::gil::get_color(Input, boost::gil::blue_t());
			const half alpha = boost::gil::get_color(Input, boost::gil::alpha_t());

			return k3d::pixel(
				red,
				green,
				std::min(green, blue),
				blue > threshold ? 1 - (blue - std::max(red, green)) : 1);
		}

		const double threshold;
	};

	void on_update_bitmap(const k3d::bitmap& Input, k3d::bitmap& Output)
	{
		boost::gil::transform_pixels(const_view(Input), view(Output), functor(m_value.pipeline_value()));
	}

	static k3d::iplugin_factory& get_factory()
	{
		static k3d::document_plugin_factory<bitmap_matte_colordiff,
			k3d::interface_list<k3d::ibitmap_source,
			k3d::interface_list<k3d::ibitmap_sink> > > factory(
				k3d::uuid(0xeefaccf2, 0x65bc4c78, 0xbd46cbdb, 0x5ca5d3e2),
				"BitmapMatteColorDiff",
				_("Create alpha channel using color difference"),
				"Bitmap",
				k3d::iplugin_factory::STABLE);

		return factory;
	}

private:
	k3d_data(double, immutable_name, change_signal, with_undo, local_storage, no_constraint, writable_property, with_serialization) m_value;
};

/////////////////////////////////////////////////////////////////////////////
// bitmap_matte_colordiff_factory

k3d::iplugin_factory& bitmap_matte_colordiff_factory()
{
	return bitmap_matte_colordiff::get_factory();
}

} // namespace libk3dbitmap


