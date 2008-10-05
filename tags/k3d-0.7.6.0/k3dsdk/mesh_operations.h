#ifndef K3DSDK_MESH_OPERATIONS_H
#define K3DSDK_MESH_OPERATIONS_H

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

#include "bounding_box3.h"
#include "mesh.h"
#include "mesh_selection.h"
#include "normal3.h"
#include "point3.h"

namespace k3d
{

class imaterial;

/// Returns a mesh containing a topological "grid" of polygons with the given number of polys in each parametric direction
const mesh create_grid(const uint_t Rows, const uint_t Columns, imaterial* const Material = 0);
/// Returns a mesh containing a topological "cylinder" of polygons with the given number of polys in each parametric direction
const mesh create_cylinder(const uint_t Rows, const uint_t Columns, imaterial* const Material = 0);

/// Returns true iff every polyhedron in the given mesh is solid volume
const bool_t is_solid(const mesh& Mesh);
/// Returns true iff every face in the given mesh is a triangle
const bool_t is_triangles(const mesh& Mesh);
/// Returns true iff every array and primitive shared pointer in the mesh is unitinialized
const bool_t is_uninitialized(const mesh& Mesh);

/// Copies the selection state of a mesh into a mesh_selection
void store_selection(const mesh& Mesh, mesh_selection& Selection);
/// Merges a mesh_selection with the current selection state in the given mesh
void merge_selection(const mesh_selection& MeshSelection, mesh& Mesh);
/// Merges a set of mesh_selection records with the current selection state in the given array
void merge_selection(const mesh_selection::records_t& Records, mesh::selection_t& Selection);
/// Clears the selection of mesh components (all array elements to 0.0)
void clear_component_selection(mesh& Mesh);

/// Returns a bounding-box containing every point in the given mesh
const bounding_box3 bounds(const mesh& Mesh);
/// Returns a bounding-box containing every point in the given collection
const bounding_box3 bounds(const mesh::points_t& Points);
/// Calculates the center (average) for an edge loop (returns the origin for degenerate cases)
const point3 center(const mesh::indices_t& EdgePoints, const mesh::indices_t& ClockwiseEdges, const mesh::points_t& Points, const uint_t EdgeIndex);
/// Calculates the normal for an edge loop (returns a zero-length normal for degenerate cases)
const normal3 normal(const mesh::indices_t& EdgePoints, const mesh::indices_t& ClockwiseEdges, const mesh::points_t& Points, const uint_t EdgeIndex);

/// Performs a deep-copy from one mesh to another (the new mesh doesn't share any memory with the old)
void deep_copy(const mesh& From, mesh& To);

/// Performs sanity-checking on a mesh, validating all constraints - returns true iff the mesh is valid
const bool_t validate(mesh& Mesh);

/// Returns true iff the given mesh contains valid point data (i.e. both point and point_selection arrays are defined)
const bool_t validate_points(const mesh& Mesh);
/// Returns true iff the given mesh contains valid point group data (i.e. every array is defined)
const bool_t validate_point_groups(const mesh& Mesh);
/// Returns true iff the given mesh contains valid linear curve group data (i.e. every array is defined)
const bool_t validate_linear_curve_groups(const mesh& Mesh);
/// Returns true iff the given mesh contains valid cubic curve group data (i.e. every array is defined)
const bool_t validate_cubic_curve_groups(const mesh& Mesh);
/// Returns true iff the given mesh contains valid nurbs curve group data (i.e. every array is defined)
const bool_t validate_nurbs_curve_groups(const mesh& Mesh);
/// Returns true iff the given mesh contains valid bilinear patch data (i.e. every array is defined)
const bool_t validate_bilinear_patches(const mesh& Mesh);
/// Returns true iff the given mesh contains valid bicubic patch data (i.e. every array is defined)
const bool_t validate_bicubic_patches(const mesh& Mesh);
/// Returns true iff the given mesh contains valid nurbs patch data (i.e. every array is defined)
const bool_t validate_nurbs_patches(const mesh& Mesh);
/// Returns true iff the given mesh contains valid polyhedron data (i.e. every array is defined)
const bool_t validate_polyhedra(const mesh& Mesh);
/// Returns true iff the given mesh contains valid blobby data (i.e. every array is defined)
const bool_t validate_blobbies(const mesh& Mesh);

/// Returns true iff the given mesh should be rendered as SDS
const bool_t is_sds(const mesh& Mesh);

/// Traverse polygonal mesh, visiting faces, loops, and points
template<typename visitor_t>
void traverse_polyhedra(const mesh& Mesh, visitor_t& Visitor)
{
	return_if_fail(validate_polyhedra(Mesh));
	const mesh::points_t& points = *Mesh.points;
	const mesh::indices_t& face_first_loops = *Mesh.polyhedra->face_first_loops;
	const mesh::counts_t& face_loop_counts = *Mesh.polyhedra->face_loop_counts;
	const mesh::indices_t& loop_first_edges = *Mesh.polyhedra->loop_first_edges;
	const mesh::indices_t& edge_points = *Mesh.polyhedra->edge_points;
	const mesh::indices_t& clockwise_edges = *Mesh.polyhedra->clockwise_edges;
	
	const uint_t face_count = face_first_loops.size();
	for(uint_t face = 0; face != face_count; ++face)
	{
		Visitor.on_face_start(face);
		const uint_t loop_begin = face_first_loops[face];
		const uint_t loop_end = loop_begin + face_loop_counts[face];

		for(uint_t loop = loop_begin; loop != loop_end; ++loop)
		{
			Visitor.on_loop(loop);
			const uint_t first_edge = loop_first_edges[loop];

			for(uint_t edge = first_edge; ; )
			{
				Visitor.on_edge(edge, edge_points[edge], edge_points[clockwise_edges[edge]], points[edge_points[edge]], points[edge_points[clockwise_edges[edge]]]);
				edge = clockwise_edges[edge];
				if(edge == first_edge)
					break;
			}
		}
		Visitor.on_face_end(face);
	}
}

/// Traverse selected mesh points
template<typename visitor_t>
void traverse_selected_points(const mesh& Mesh, visitor_t& Visitor)
{
	return_if_fail(validate_points(Mesh));
	for (uint_t point = 0; point != Mesh.points->size(); ++point)
	{
		if (Mesh.point_selection->at(point))
		{
			Visitor(point, Mesh.points->at(point));
		}
	}
}

/// For each selected edge, visit the start and end point (multiple visits per point possible!)
template<typename visitor_t>
void traverse_selected_edge_points(const mesh& Mesh, visitor_t& Visitor)
{
	return_if_fail(validate_polyhedra(Mesh));
	const mesh::points_t& points = *Mesh.points;
	const mesh::indices_t& edge_points = *Mesh.polyhedra->edge_points;
	const mesh::indices_t& clockwise_edges = *Mesh.polyhedra->clockwise_edges;
	const mesh::selection_t& edge_selection = *Mesh.polyhedra->edge_selection;
	for (uint_t edge = 0; edge != edge_points.size(); ++edge)
	{
		if (edge_selection[edge])
		{
			Visitor(edge_points[edge], points[edge_points[edge]]);
			Visitor(edge_points[clockwise_edges[edge]], points[edge_points[clockwise_edges[edge]]]);
		}
	}
}

// For each selected face, visit all of its points (multiple visits per point possible!)
template<typename visitor_t>
void traverse_selected_face_points(const mesh& Mesh, visitor_t& Visitor)
{
	return_if_fail(validate_polyhedra(Mesh));
	const mesh::points_t& points = *Mesh.points;
	const mesh::indices_t& face_first_loops = *Mesh.polyhedra->face_first_loops;
	const mesh::counts_t& face_loop_counts = *Mesh.polyhedra->face_loop_counts;
	const mesh::indices_t& loop_first_edges = *Mesh.polyhedra->loop_first_edges;
	const mesh::indices_t& edge_points = *Mesh.polyhedra->edge_points;
	const mesh::indices_t& clockwise_edges = *Mesh.polyhedra->clockwise_edges;
	const mesh::selection_t& face_selection = *Mesh.polyhedra->face_selection;
	for(uint_t face = 0; face != face_first_loops.size(); ++face)
	{
		if (!face_selection[face])
			continue;
		
		const uint_t loop_begin = face_first_loops[face];
		const uint_t loop_end = loop_begin + face_loop_counts[face];
		for(uint_t loop = loop_begin; loop != loop_end; ++loop)
		{
			const uint_t first_edge = loop_first_edges[loop];
			for(uint_t edge = first_edge; ; )
			{
				Visitor(edge_points[edge], points[edge_points[edge]]);

				edge = clockwise_edges[edge];
				if(edge == first_edge)
					break;
			}
		}
	}
}

} // namespace k3d

#endif // K3DSDK_MESH_OPERATIONS_H
