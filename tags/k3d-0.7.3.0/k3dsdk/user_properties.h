#ifndef K3DSDK_USER_PROPERTIES_H
#define K3DSDK_USER_PROPERTIES_H

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
	\author Tim Shead (tshead@k-3d.com)
*/

#include "ipersistent_container.h"
#include "iproperty_collection.h"
#include "istate_container.h"
#include "types.h"

namespace k3d
{

namespace user
{

////////////////////////////////////////////////////////////////////////////////////////////////////
// property_container

/// istate_container implementation that handles undo/redo for the set of user properties within a property collection
class property_container :
	public istate_container
{
public:
	property_container(iunknown& Owner);
	void restore_state();

private:
	iproperty_collection* const m_property_collection;
	ipersistent_container* const m_persistent_container;
	iproperty_collection::properties_t m_user_properties;
	ipersistent_container::named_children_t m_persistent_properties;
};

} // namespace user

} // namespace k3d

#endif // K3DSDK_USER_PROPERTIES_H

