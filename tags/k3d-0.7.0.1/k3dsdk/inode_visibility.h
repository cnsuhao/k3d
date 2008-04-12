#ifndef K3DSDK_INODE_VISIBILITY_H
#define K3DSDK_INODE_VISIBILITY_H

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
	\author Tim Shead (tshead@k-3d.com)
*/

#include "iunknown.h"

namespace k3d
{

class iproperty;

/// Abstract interface for "render engine" objects that maintain an explicit list of visible nodes
class inode_visibility :
	public virtual iunknown
{
public:
	/// Returns the property that will be used to store the list of visible nodes
	virtual iproperty& visible_nodes() = 0;

protected:
	inode_visibility() {}
	inode_visibility(const inode_visibility&) {}
	inode_visibility& operator=(const inode_visibility&) { return *this; }
	virtual ~inode_visibility() {}
};

} // namespace k3d

#endif // K3DSDK_INODE_VISIBILITY_H
