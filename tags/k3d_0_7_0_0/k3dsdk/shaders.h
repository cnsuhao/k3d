#ifndef K3DSDK_SHADERS_H
#define K3DSDK_SHADERS_H

// K-3D
// Copyright (c) 1995-2004, Timothy M. Shead
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
		\brief defines shader helpers
		\author Tim Shead (tshead@k-3d.com)
*/

#include "sl.h"
#include <memory>

namespace k3d
{

/// Recompiles a shader from source iff it's out-of-date
bool compile_shader(const filesystem::path& SourceFile, const std::string& RenderEngineType, const std::string& RenderEngine);

} // namespace k3d

#endif // !K3DSDK_SHADERS_H

