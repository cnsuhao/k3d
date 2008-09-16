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
	\author Carlos Andres Dominguez Caballero (carlosadc@gmail.com)
*/

#ifndef __integration_h__
#define __integration_h__

#include <k3dsdk/idocument.h>
#include <k3dsdk/mesh.h>

namespace module
{

namespace f3ds
{

namespace io
{

class f3dsParser
{
public:
	f3dsParser(const char *filename, k3d::mesh &Mesh);
};

} // namespace io

} // namespace f3ds

} // namespace module

#endif

