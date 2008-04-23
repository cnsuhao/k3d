#ifndef K3DSDK_HINTS_H
#define K3DSDK_HINTS_H

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

/** \file
	\author Timothy M. Shead (tshead@k-3d.com)
*/

#include "algebra.h"
#include "iunknown.h"
#include "mesh_selection.h"
#include "mesh.h"

#include <boost/any.hpp>

#include <iosfwd>

namespace k3d
{

class mesh_selection;
	
/// Namespace reserved for "hint" object that pass metadata about an upstream change to downstream nodes
namespace hint
{

/// Hint object that indicates that an object's selection state has changed
class selection_changed_t :
	public iunknown
{
};

/// Convenience function that returns a reference to a static instance of selection_changed_t
selection_changed_t* selection_changed();

/// Hint object that indicates that a mesh's geometry (the locations of its points) has changed
class mesh_geometry_changed_t :
	public iunknown
{
public:
	mesh_geometry_changed_t() {}
	/// Indices of the points affected by the change
	k3d::mesh::indices_t changed_points;
	/// Transformation matrix used for the change
	k3d::matrix4 transformation_matrix;
};

/// Convenience function that returns a reference to a static instance of mesh_geometry_changed_t
mesh_geometry_changed_t* mesh_geometry_changed();

/// Hint object that indicates that a mesh's topology has changed
class mesh_topology_changed_t :
	public iunknown
{
};

/// Convenience function that returns a reference to a static instance of mesh_topology_changed_t
mesh_topology_changed_t* mesh_topology_changed();

/// Hint object that indicates that only the mesh's address has changed and provides access to the old addresses
class mesh_address_changed_t :
	public iunknown
{
public:
	boost::shared_ptr<const k3d::mesh::points_t> old_points_address;
	boost::shared_ptr<const k3d::mesh::indices_t> old_edge_points_address;
	boost::shared_ptr<const k3d::mesh::indices_t> old_face_first_loops_address;
};

/// Convenience function that returns a reference to a static instance of mesh_address_changed_t
mesh_address_changed_t* mesh_address_changed();

/// Hint object that indicates a mesh was deleted
class mesh_deleted_t :
	public iunknown
{
};

/// Convenience function that returns a reference to a static instance of mesh_deleted_t
mesh_deleted_t* mesh_deleted();

/// iostream-compatible manipulator object that serializes information about a hint object
class print
{
public:
	print(iunknown* Hint);

	iunknown* const hint;
};

/// Stream serialization
std::ostream& operator<<(std::ostream& Stream, const print& RHS);

/// Common interface for classes that need to process incoming hints
class hint_processor
{
public:
	/// Process the given hint, calling the required on_... method
	virtual void process(const k3d::mesh& Mesh, k3d::iunknown* Hint);
	
	/// Process the given hint (stored as boost::any), calling the required on_... method
	virtual void process(const k3d::mesh& Mesh, boost::any& Hint);

protected:
	/// Called when only the geometry changed
	virtual void on_geometry_changed(const k3d::mesh& Mesh, k3d::iunknown* Hint) {}
	
	/// Called when only the selection changed
	virtual void on_selection_changed(const k3d::mesh& Mesh, k3d::iunknown* Hint) {}
	
	/// Called when the mesh topology changed
	virtual void on_topology_changed(const k3d::mesh& Mesh, k3d::iunknown* Hint) {}
	
	/// Called when the mesh address changed
	virtual void on_address_changed(const k3d::mesh& Mesh, k3d::iunknown* Hint) {}
	
	/// Called when the mesh is deleted
	virtual void on_mesh_deleted(const k3d::mesh& Mesh, k3d::iunknown* Hint) {}
	
	/// Called when an unknown hint is encountered
	virtual void on_unknown_change(const k3d::mesh& Mesh, k3d::iunknown* Hint)
	{
		on_topology_changed(Mesh, Hint);
	}
};

} // namespace hint

} // namespace k3d

#endif // K3DSDK_HINTS_H
