#ifndef K3DSDK_CONST_NAMED_ARRAYS_PYTHON_H
#define K3DSDK_CONST_NAMED_ARRAYS_PYTHON_H

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

#include <k3dsdk/mesh.h>
#include <boost/python.hpp>

namespace k3d
{

namespace python
{

//////////////////////////////////////////////////////////////////////////////////////////////////
// const_named_arrays

class const_named_arrays
{
public:
	const_named_arrays(const k3d::mesh::named_arrays& NamedArrays);

	boost::python::list array_names();
	boost::python::object array(const std::string& Name);

	int len();
	boost::python::object get_item(int item);

	static void define_class();

private:
	boost::python::object wrap_array(const k3d::array* const Array);
	const k3d::mesh::named_arrays& wrapped;
};

} // namespace python

} // namespace k3d

#endif // !K3DSDK_CONST_NAMED_ARRAYS_PYTHON_H

