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

/** \file
		\author Tim Shead (tshead@k-3d.com)
*/

#include <k3d-i18n-config.h>
#include <k3d-platform-config.h>
#include <k3dsdk/document_plugin_factory.h>
#include <k3dsdk/iuser_interface.h>
#include <k3dsdk/measurement.h>
#include <k3dsdk/node.h>
#include <k3dsdk/persistent.h>
#include <k3dsdk/user_interface.h>

#include <limits>

#ifdef K3D_API_WIN32
#include <sys/types.h>
#include <sys/timeb.h>
#else // K3D_API_WIN32
#include <sys/time.h>
#endif // !K3D_API_WIN32

namespace module
{

namespace time
{

/////////////////////////////////////////////////////////////////////////////
// realtime_source

class realtime_source :
	public k3d::persistent<k3d::node>
{
	typedef k3d::persistent<k3d::node>  base;

public:
	realtime_source(k3d::iplugin_factory& Factory, k3d::idocument& Document) :
		base(Factory, Document),
		m_frame_rate(init_owner(*this) + init_name("frame_rate") + init_label(_("Frame rate")) + init_description(_("Frame rate")) + init_value(1.0) + init_step_increment(1) + init_units(typeid(k3d::measurement::scalar)) + init_constraint(constraint::minimum(std::numeric_limits<double>::epsilon()))),
		m_time(init_owner(*this) + init_name("time") + init_label(_("Time")) + init_description(_("Time")) + init_slot(sigc::mem_fun(*this, &realtime_source::get_time)))
	{
		m_frame_rate.changed_signal().connect(sigc::mem_fun(*this, &realtime_source::on_reset_source));
		on_reset_source(0);
	}

	~realtime_source()
	{
		m_timeout_connection.disconnect();
	}

	void on_reset_source(k3d::iunknown*)
	{
		m_timeout_connection.disconnect();
		m_timeout_connection = k3d::user_interface().get_timer(m_frame_rate.pipeline_value(), sigc::bind(m_time.make_reset_slot(), static_cast<k3d::iunknown*>(0)));
	}

	double get_time()
	{
#ifdef K3D_API_WIN32
		timeb tv;
		ftime(&tv);
		return tv.time + (static_cast<double>(tv.millitm) / 1000);
#else // K3D_API_WIN32
		timeval tv;
		gettimeofday(&tv, 0);
		return tv.tv_sec + static_cast<double>(tv.tv_usec) / 1000000;
#endif // !K3D_API_WIN32
	}

	static k3d::iplugin_factory& get_factory()
	{
		static k3d::document_plugin_factory<realtime_source > factory(
			k3d::uuid(0xf7c119e7, 0xd87b465a, 0x9b7e0289, 0x82d91702),
			"RealtimeSource",
			_("Provides a realtime time source for a K-3D document"),
			"Animation",
			k3d::iplugin_factory::STABLE);

		return factory;
	}

private:
	k3d_data(double, immutable_name, change_signal, with_undo, local_storage, with_constraint, measurement_property, with_serialization) m_frame_rate;
	k3d_data(double, immutable_name, change_signal, no_undo, computed_storage, no_constraint, read_only_property, no_serialization) m_time;

	sigc::connection m_timeout_connection;
};

/////////////////////////////////////////////////////////////////////////////
// realtime_source_factory

k3d::iplugin_factory& realtime_source_factory()
{
	return realtime_source::get_factory();
}

} // namespace time

} // namespace module

