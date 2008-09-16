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

#include "imaterial.h"
#include "mesh_operations.h"
#include "mesh_selection.h"
#include "polyhedra.h"
#include "shared_pointer.h"
#include "vectors.h"

namespace k3d
{

namespace detail
{

void store_selection(const boost::shared_ptr<const mesh::selection_t>& MeshSelection, mesh_selection::records_t& Records)
{
	if(!MeshSelection)
		return;

	const mesh::selection_t& selection_weight = *MeshSelection;

	const uint_t selection_begin = 0;
	const uint_t selection_end = selection_begin + selection_weight.size();
	for(uint_t selection = selection_begin; selection != selection_end; ++selection)
		Records.push_back(mesh_selection::record(selection, selection+1, selection_weight[selection]));
}

void merge_selection(const mesh_selection::records_t& Records, mesh::selection_t& Selection)
{
	const uint_t selection_count = Selection.size();

	for(mesh_selection::records_t::const_iterator record = Records.begin(); record != Records.end(); ++record)
	{
		if(record->begin >= selection_count)
			break;

		mesh::selection_t::iterator begin = Selection.begin() + record->begin;
		mesh::selection_t::iterator end = Selection.begin() + std::min(selection_count, record->end);
		std::fill(begin, end, record->weight);
	}
}

template<typename gprims_type>
void merge_selection(const mesh_selection::records_t& Records, const gprims_type& GPrims, boost::shared_ptr<const mesh::selection_t>& Selection)
{
	return_if_fail(GPrims);

	const uint_t gprim_count = GPrims->size();

	if(!Selection || Selection->size() != gprim_count)
		Selection.reset(new mesh::selection_t(gprim_count));

	mesh::selection_t& selection = *make_unique(Selection);

	detail::merge_selection(Records, selection);
}

template<typename gprims_type>
void clear_selection(const gprims_type& GPrims, boost::shared_ptr<const mesh::selection_t>& Selection)
{
	return_if_fail(GPrims);

	const uint_t gprim_count = GPrims->size();

	Selection.reset(new mesh::selection_t(gprim_count, 0.0));
}

} // namespace detail

const mesh create_grid(const uint_t Rows, const uint_t Columns, imaterial* const Material)
{
	return_val_if_fail(Rows, mesh());
	return_val_if_fail(Columns, mesh());

	mesh result;

	const unsigned long rows = Rows;
	const unsigned long columns = Columns;
	const unsigned long point_rows = rows + 1;
	const unsigned long point_columns = columns + 1;

	boost::shared_ptr<mesh::polyhedra_t> polyhedra(new mesh::polyhedra_t());
	boost::shared_ptr<mesh::indices_t> first_faces(new mesh::indices_t(1, 0));
	boost::shared_ptr<mesh::counts_t> face_counts(new mesh::counts_t(1, rows * columns));
	boost::shared_ptr<mesh::polyhedra_t::types_t> types(new mesh::polyhedra_t::types_t(1, mesh::polyhedra_t::POLYGONS));
	boost::shared_ptr<mesh::indices_t> face_first_loops(new mesh::indices_t(rows * columns));
	boost::shared_ptr<mesh::counts_t> face_loop_counts(new mesh::counts_t(rows * columns, 1));
	boost::shared_ptr<mesh::selection_t> face_selection(new mesh::selection_t(rows * columns, 0.0));
	boost::shared_ptr<mesh::materials_t> face_materials(new mesh::materials_t(rows * columns, Material));
	boost::shared_ptr<mesh::indices_t> loop_first_edges(new mesh::indices_t(rows * columns));
	boost::shared_ptr<mesh::indices_t> edge_points(new mesh::indices_t(4 * rows * columns));
	boost::shared_ptr<mesh::indices_t> clockwise_edges(new mesh::indices_t(4 * rows * columns));
	boost::shared_ptr<mesh::selection_t> edge_selection(new mesh::selection_t(4 * rows * columns, 0.0));
	boost::shared_ptr<mesh::points_t> points(new mesh::points_t(point_rows * point_columns));
	boost::shared_ptr<mesh::selection_t> point_selection(new mesh::selection_t(point_rows * point_columns, 0.0));

	mesh::indices_t::iterator face_first_loop = face_first_loops->begin();
	mesh::indices_t::iterator loop_first_edge = loop_first_edges->begin();
	mesh::indices_t::iterator edge_point = edge_points->begin();
	mesh::indices_t::iterator clockwise_edge = clockwise_edges->begin();

	uint_t face_number = 0;

	for(unsigned long row = 0; row != rows; ++row)
	{
		for(unsigned long column = 0; column != columns; ++column)
		{
			*face_first_loop++ = face_number;

			*loop_first_edge++ = 4 * face_number;

			*edge_point++ = column + (row * point_columns);
			*edge_point++ = column + (row * point_columns) + 1;
			*edge_point++ = column + ((row + 1) * point_columns) + 1;
			*edge_point++ = column + ((row + 1) * point_columns);

			*clockwise_edge++ = (4 * face_number) + 1;
			*clockwise_edge++ = (4 * face_number) + 2;
			*clockwise_edge++ = (4 * face_number) + 3;
			*clockwise_edge++ = (4 * face_number);

			++face_number;
		}
	}

	polyhedra->first_faces = first_faces;
	polyhedra->face_counts = face_counts;
	polyhedra->types = types;
	polyhedra->face_first_loops = face_first_loops;
	polyhedra->face_loop_counts = face_loop_counts;
	polyhedra->face_selection = face_selection;
	polyhedra->face_materials = face_materials;
	polyhedra->loop_first_edges = loop_first_edges;
	polyhedra->edge_points = edge_points;
	polyhedra->clockwise_edges = clockwise_edges;
	polyhedra->edge_selection = edge_selection;

	result.polyhedra = polyhedra;
	result.points = points;
	result.point_selection = point_selection;

	return result;
}

const mesh create_cylinder(const uint_t Rows, const uint_t Columns, imaterial* const Material)
{
	return_val_if_fail(Rows, mesh());
	return_val_if_fail(Columns > 1, mesh());

	mesh result;

	const unsigned long rows = Rows;
	const unsigned long columns = Columns;
	const unsigned long point_rows = rows + 1;
	const unsigned long point_columns = columns;

	boost::shared_ptr<mesh::polyhedra_t> polyhedra(new mesh::polyhedra_t());
	boost::shared_ptr<mesh::indices_t> first_faces(new mesh::indices_t(1, 0));
	boost::shared_ptr<mesh::counts_t> face_counts(new mesh::counts_t(1, rows * columns));
	boost::shared_ptr<mesh::polyhedra_t::types_t> types(new mesh::polyhedra_t::types_t(1, mesh::polyhedra_t::POLYGONS));
	boost::shared_ptr<mesh::indices_t> face_first_loops(new mesh::indices_t(rows * columns));
	boost::shared_ptr<mesh::counts_t> face_loop_counts(new mesh::counts_t(rows * columns, 1));
	boost::shared_ptr<mesh::selection_t> face_selection(new mesh::selection_t(rows * columns, 0.0));
	boost::shared_ptr<mesh::materials_t> face_materials(new mesh::materials_t(rows * columns, Material));
	boost::shared_ptr<mesh::indices_t> loop_first_edges(new mesh::indices_t(rows * columns));
	boost::shared_ptr<mesh::indices_t> edge_points(new mesh::indices_t(4 * rows * columns));
	boost::shared_ptr<mesh::indices_t> clockwise_edges(new mesh::indices_t(4 * rows * columns));
	boost::shared_ptr<mesh::selection_t> edge_selection(new mesh::selection_t(4 * rows * columns, 0.0));
	boost::shared_ptr<mesh::points_t> points(new mesh::points_t(point_rows * point_columns));
	boost::shared_ptr<mesh::selection_t> point_selection(new mesh::selection_t(point_rows * point_columns, 0.0));

	mesh::indices_t::iterator face_first_loop = face_first_loops->begin();
	mesh::indices_t::iterator loop_first_edge = loop_first_edges->begin();
	mesh::indices_t::iterator edge_point = edge_points->begin();
	mesh::indices_t::iterator clockwise_edge = clockwise_edges->begin();

	uint_t face_number = 0;

	for(unsigned long row = 0; row != rows; ++row)
	{
		for(unsigned long column = 0; column != columns; ++column)
		{
			*face_first_loop++ = face_number;

			*loop_first_edge++ = 4 * face_number;

			*edge_point++ = column + (row * point_columns);
			*edge_point++ = (column + (row * point_columns) + 1) % point_columns;
			*edge_point++ = (column + ((row + 1) * point_columns) + 1) % point_columns;
			*edge_point++ = column + ((row + 1) * point_columns);

			*clockwise_edge++ = (4 * face_number) + 1;
			*clockwise_edge++ = (4 * face_number) + 2;
			*clockwise_edge++ = (4 * face_number) + 3;
			*clockwise_edge++ = (4 * face_number);

			++face_number;
		}
	}

	polyhedra->first_faces = first_faces;
	polyhedra->face_counts = face_counts;
	polyhedra->types = types;
	polyhedra->face_first_loops = face_first_loops;
	polyhedra->face_loop_counts = face_loop_counts;
	polyhedra->face_selection = face_selection;
	polyhedra->face_materials = face_materials;
	polyhedra->loop_first_edges = loop_first_edges;
	polyhedra->edge_points = edge_points;
	polyhedra->clockwise_edges = clockwise_edges;
	polyhedra->edge_selection = edge_selection;

	result.polyhedra = polyhedra;
	result.points = points;
	result.point_selection = point_selection;

	return result;
}

const bool_t is_solid(const mesh& Mesh)
{
	const uint_t polyhedron_begin = 0;
	const uint_t polyhedron_end = polyhedron_begin + polyhedron::count(Mesh);
	for(uint_t polyhedron = polyhedron_begin; polyhedron != polyhedron_end; ++polyhedron)
	{
		if(!polyhedron::is_solid(Mesh, polyhedron))
			return false;
	}

	return true;
}

const bool_t is_triangles(const mesh& Mesh)
{
	if(!validate_polyhedra(Mesh))
		return true;

	const k3d::mesh::indices_t& face_first_loops = *Mesh.polyhedra->face_first_loops;
	const k3d::mesh::indices_t& loop_first_edges = *Mesh.polyhedra->loop_first_edges;
	const k3d::mesh::indices_t& edge_points = *Mesh.polyhedra->edge_points;
	const k3d::mesh::indices_t& clockwise_edges = *Mesh.polyhedra->clockwise_edges;

	const uint_t face_begin = 0;
	const uint_t face_end = face_begin + face_first_loops.size();
	for(uint_t face = face_begin; face != face_end; ++face)
	{
		uint_t edge_count = 0;
		const uint_t first_edge = loop_first_edges[face_first_loops[face]];
		for(uint_t edge = first_edge; ; )
		{
			++edge_count;

			edge = clockwise_edges[edge];
			if(edge == first_edge)
				break;
		}

		if(edge_count != 3)
			return false;
	}

	return true;
}

const bool_t is_uninitialized(const mesh& Mesh)
{
	if (Mesh.points)
		return false;
	if (Mesh.bicubic_patches)
		return false;
	if (Mesh.bilinear_patches)
		return false;
	if (Mesh.blobbies)
		return false;
	if (Mesh.cubic_curve_groups)
		return false;
	if (Mesh.linear_curve_groups)
		return false;
	if (Mesh.nurbs_curve_groups)
		return false;
	if (Mesh.nurbs_patches)
		return false;
	if (Mesh.point_groups)
		return false;
	if (Mesh.polyhedra)
		return false;
	if (Mesh.point_selection)
		return false;
}

/////////////////////////////////////////////////////////////////////////////
// store_selection

void store_selection(const mesh& Mesh, mesh_selection& Selection)
{
	Selection.points.clear();
	Selection.edges.clear();
	Selection.faces.clear();
	Selection.linear_curves.clear();
	Selection.cubic_curves.clear();
	Selection.nurbs_curves.clear();
	Selection.bilinear_patches.clear();
	Selection.bicubic_patches.clear();
	Selection.nurbs_patches.clear();

	detail::store_selection(Mesh.point_selection, Selection.points);

	if(Mesh.polyhedra)
	{
		detail::store_selection(Mesh.polyhedra->edge_selection, Selection.edges);
		detail::store_selection(Mesh.polyhedra->face_selection, Selection.faces);
	}

	if(Mesh.linear_curve_groups)
	{
		detail::store_selection(Mesh.linear_curve_groups->curve_selection, Selection.linear_curves);
	}

	if(Mesh.cubic_curve_groups)
	{
		detail::store_selection(Mesh.cubic_curve_groups->curve_selection, Selection.cubic_curves);
	}

	if(Mesh.nurbs_curve_groups)
	{
		detail::store_selection(Mesh.nurbs_curve_groups->curve_selection, Selection.nurbs_curves);
	}

	if(Mesh.bilinear_patches)
	{
		detail::store_selection(Mesh.bilinear_patches->patch_selection, Selection.bilinear_patches);
	}

	if(Mesh.bicubic_patches)
	{
		detail::store_selection(Mesh.bicubic_patches->patch_selection, Selection.bicubic_patches);
	}

	if(Mesh.nurbs_patches)
	{
		detail::store_selection(Mesh.nurbs_patches->patch_selection, Selection.nurbs_patches);
	}

	/** \todo Store blobby and quadric selections */
}

/////////////////////////////////////////////////////////////////////////////
// merge_selection

void merge_selection(const mesh_selection& MeshSelection, mesh& Mesh)
{
	detail::merge_selection(MeshSelection.points, Mesh.points, Mesh.point_selection);

	if(Mesh.polyhedra && Mesh.polyhedra->edge_points)
	{
		k3d::mesh::polyhedra_t* const polyhedra = make_unique(Mesh.polyhedra);
		detail::merge_selection(MeshSelection.edges, polyhedra->edge_points, polyhedra->edge_selection);
	}

	if(Mesh.polyhedra && Mesh.polyhedra->face_first_loops)
	{
		k3d::mesh::polyhedra_t* const polyhedra = make_unique(Mesh.polyhedra);
		detail::merge_selection(MeshSelection.faces, polyhedra->face_first_loops, polyhedra->face_selection);
	}

	if(Mesh.linear_curve_groups)
	{
		k3d::mesh::linear_curve_groups_t* const linear_curve_groups = make_unique(Mesh.linear_curve_groups);
		detail::merge_selection(MeshSelection.linear_curves, linear_curve_groups->curve_first_points, linear_curve_groups->curve_selection);
	}

	if(Mesh.cubic_curve_groups)
	{
		k3d::mesh::cubic_curve_groups_t* const cubic_curve_groups = make_unique(Mesh.cubic_curve_groups);
		detail::merge_selection(MeshSelection.cubic_curves, cubic_curve_groups->curve_first_points, cubic_curve_groups->curve_selection);
	}

	if(Mesh.nurbs_curve_groups)
	{
		k3d::mesh::nurbs_curve_groups_t* const nurbs_curve_groups = make_unique(Mesh.nurbs_curve_groups);
		detail::merge_selection(MeshSelection.nurbs_curves, nurbs_curve_groups->curve_first_points, nurbs_curve_groups->curve_selection);
	}

	if(Mesh.bilinear_patches)
	{
		k3d::mesh::bilinear_patches_t* const bilinear_patches = make_unique(Mesh.bilinear_patches);
		detail::merge_selection(MeshSelection.bilinear_patches, bilinear_patches->patch_materials, bilinear_patches->patch_selection);
	}

	if(Mesh.bicubic_patches)
	{
		k3d::mesh::bicubic_patches_t* const bicubic_patches = make_unique(Mesh.bicubic_patches);
		detail::merge_selection(MeshSelection.bicubic_patches, bicubic_patches->patch_materials, bicubic_patches->patch_selection);
	}

	if(Mesh.nurbs_patches)
	{
		k3d::mesh::nurbs_patches_t* const nurbs_patches = make_unique(Mesh.nurbs_patches);
		detail::merge_selection(MeshSelection.nurbs_patches, nurbs_patches->patch_materials, nurbs_patches->patch_selection);
	}
}

void merge_selection(const mesh_selection::records_t& Records, mesh::selection_t& Selection)
{
	detail::merge_selection(Records, Selection);
}

void clear_component_selection(mesh& Mesh)
{
	detail::clear_selection(Mesh.points, Mesh.point_selection);

	if(Mesh.polyhedra && Mesh.polyhedra->edge_points)
	{
		k3d::mesh::polyhedra_t* const polyhedra = make_unique(Mesh.polyhedra);
		detail::clear_selection(polyhedra->edge_points, polyhedra->edge_selection);
	}

	if(Mesh.polyhedra && Mesh.polyhedra->face_first_loops)
	{
		k3d::mesh::polyhedra_t* const polyhedra = make_unique(Mesh.polyhedra);
		detail::clear_selection(polyhedra->face_first_loops, polyhedra->face_selection);
	}

	if(Mesh.linear_curve_groups)
	{
		k3d::mesh::linear_curve_groups_t* const linear_curve_groups = make_unique(Mesh.linear_curve_groups);
		detail::clear_selection(linear_curve_groups->curve_first_points, linear_curve_groups->curve_selection);
	}

	if(Mesh.cubic_curve_groups)
	{
		k3d::mesh::cubic_curve_groups_t* const cubic_curve_groups = make_unique(Mesh.cubic_curve_groups);
		detail::clear_selection(cubic_curve_groups->curve_first_points, cubic_curve_groups->curve_selection);
	}

	if(Mesh.nurbs_curve_groups)
	{
		k3d::mesh::nurbs_curve_groups_t* const nurbs_curve_groups = make_unique(Mesh.nurbs_curve_groups);
		detail::clear_selection(nurbs_curve_groups->curve_first_points, nurbs_curve_groups->curve_selection);
	}

	if(Mesh.bilinear_patches)
	{
		k3d::mesh::bilinear_patches_t* const bilinear_patches = make_unique(Mesh.bilinear_patches);
		detail::clear_selection(bilinear_patches->patch_materials, bilinear_patches->patch_selection);
	}

	if(Mesh.bicubic_patches)
	{
		k3d::mesh::bicubic_patches_t* const bicubic_patches = make_unique(Mesh.bicubic_patches);
		detail::clear_selection(bicubic_patches->patch_materials, bicubic_patches->patch_selection);
	}

	if(Mesh.nurbs_patches)
	{
		k3d::mesh::nurbs_patches_t* const nurbs_patches = make_unique(Mesh.nurbs_patches);
		detail::clear_selection(nurbs_patches->patch_materials, nurbs_patches->patch_selection);
	}
}

/////////////////////////////////////////////////////////////////////////////
// bounds

const bounding_box3 bounds(const k3d::mesh& Mesh)
{
	bounding_box3 results;

	if(Mesh.points)
        return bounds(*Mesh.points);

	return results;
}

/////////////////////////////////////////////////////////////////////////////
// bounds

const bounding_box3 bounds(const k3d::mesh::points_t& Points)
{
	bounding_box3 results;

    const uint_t point_begin = 0;
    const uint_t point_end = point_begin + Points.size();
    for(uint_t point = point_begin; point != point_end; ++point)
        results.insert(Points[point]);

	return results;
}

//////////////////////////////////////////////////////////////////////////////
// center

const point3 center(const mesh::indices_t& EdgePoints, const mesh::indices_t& ClockwiseEdges, const mesh::points_t& Points, const uint_t EdgeIndex)
{
	point3 result(0, 0, 0);

	uint_t count = 0;
	for(uint_t edge = EdgeIndex; ; )
	{
		result += to_vector(Points[EdgePoints[edge]]);
		++count;

		edge = ClockwiseEdges[edge];
		if(edge == EdgeIndex)
			break;
	}

	if(count)
		result /= static_cast<double>(count);

	return result;
}

//////////////////////////////////////////////////////////////////////////////
// normal

const normal3 normal(const mesh::indices_t& EdgePoints, const mesh::indices_t& ClockwiseEdges, const mesh::points_t& Points, const uint_t EdgeIndex)
{
	// Calculates the normal for an edge loop using the summation method, which is more robust than the three-point methods (handles zero-length edges)
	normal3 result(0, 0, 0);

	for(uint_t edge = EdgeIndex; ; )
	{
		const point3& i = Points[EdgePoints[edge]];
		const point3& j = Points[EdgePoints[ClockwiseEdges[edge]]];

		result[0] += (i[1] + j[1]) * (j[2] - i[2]);
		result[1] += (i[2] + j[2]) * (j[0] - i[0]);
		result[2] += (i[0] + j[0]) * (j[1] - i[1]);

		edge = ClockwiseEdges[edge];
		if(edge == EdgeIndex)
			break;
	}

	return 0.5 * result;
}

////////////////////////////////////////////////////////////////////////////////
// deep_copy

void deep_copy(const mesh& From, mesh& To)
{
	To = From;
	assert_not_implemented(); // Need to ensure that all storage is unique
}


////////////////////////////////////////////////////////////////////////////////
// validate

/** \todo Handle all gprim types */
const bool_t validate(mesh& Mesh)
{
	bool_t result = true;

	if(Mesh.points && !validate_points(Mesh))
		result = false;

	if(Mesh.point_groups && !validate_point_groups(Mesh))
		result = false;

	if(Mesh.linear_curve_groups && !validate_linear_curve_groups(Mesh))
		result = false;

	if(Mesh.cubic_curve_groups && !validate_cubic_curve_groups(Mesh))
		result = false;

	if(Mesh.nurbs_curve_groups && !validate_nurbs_curve_groups(Mesh))
		result = false;

	if(Mesh.bilinear_patches && !validate_bilinear_patches(Mesh))
		result = false;

	if(Mesh.bicubic_patches && !validate_bicubic_patches(Mesh))
		result = false;

	if(Mesh.nurbs_patches && !validate_nurbs_patches(Mesh))
		result = false;

	if(Mesh.polyhedra && !validate_polyhedra(Mesh))
		result = false;

	if(Mesh.blobbies && !validate_blobbies(Mesh))
		result = false;

	return result;
}

const bool_t validate_points(const mesh& Mesh)
{
	if(!Mesh.points)
		return false;

	// If the mesh has points, it must have storage for point selections
	return_val_if_fail(Mesh.point_selection, false);

	// The point and point selection arrays must be the same length
	return_val_if_fail(Mesh.points->size() == Mesh.point_selection->size(), false);

	return true;
}

const bool_t validate_point_groups(const mesh& Mesh)
{
	if(!Mesh.point_groups)
		return false;

	return_val_if_fail(validate_points(Mesh), false);
	return_val_if_fail(Mesh.point_groups->first_points, false);
	return_val_if_fail(Mesh.point_groups->point_counts, false);
	return_val_if_fail(Mesh.point_groups->points, false);
	return_val_if_fail(Mesh.point_groups->materials, false);

	return true;
}

const bool_t validate_linear_curve_groups(const mesh& Mesh)
{
	if(!Mesh.linear_curve_groups)
		return false;

	return_val_if_fail(validate_points(Mesh), false);
	return_val_if_fail(Mesh.linear_curve_groups->first_curves, false);
	return_val_if_fail(Mesh.linear_curve_groups->curve_counts, false);
	return_val_if_fail(Mesh.linear_curve_groups->periodic_curves, false);
	return_val_if_fail(Mesh.linear_curve_groups->curve_first_points, false);
	return_val_if_fail(Mesh.linear_curve_groups->curve_point_counts, false);
	return_val_if_fail(Mesh.linear_curve_groups->curve_points, false);
	return_val_if_fail(Mesh.linear_curve_groups->curve_selection, false);
	return_val_if_fail(Mesh.linear_curve_groups->materials, false);

	return true;
}

const bool_t validate_cubic_curve_groups(const mesh& Mesh)
{
	if(!Mesh.cubic_curve_groups)
		return false;

	return_val_if_fail(validate_points(Mesh), false);
	return_val_if_fail(Mesh.cubic_curve_groups->first_curves, false);
	return_val_if_fail(Mesh.cubic_curve_groups->curve_counts, false);
	return_val_if_fail(Mesh.cubic_curve_groups->periodic_curves, false);
	return_val_if_fail(Mesh.cubic_curve_groups->curve_first_points, false);
	return_val_if_fail(Mesh.cubic_curve_groups->curve_point_counts, false);
	return_val_if_fail(Mesh.cubic_curve_groups->curve_points, false);
	return_val_if_fail(Mesh.cubic_curve_groups->curve_selection, false);
	return_val_if_fail(Mesh.cubic_curve_groups->materials, false);

	return true;
}

const bool_t validate_nurbs_curve_groups(const mesh& Mesh)
{
	if(!Mesh.nurbs_curve_groups)
		return false;

	return_val_if_fail(validate_points(Mesh), false);
	return_val_if_fail(Mesh.nurbs_curve_groups->first_curves, false);
	return_val_if_fail(Mesh.nurbs_curve_groups->curve_counts, false);
	return_val_if_fail(Mesh.nurbs_curve_groups->curve_first_points, false);
	return_val_if_fail(Mesh.nurbs_curve_groups->curve_point_counts, false);
	return_val_if_fail(Mesh.nurbs_curve_groups->curve_orders, false);
	return_val_if_fail(Mesh.nurbs_curve_groups->curve_first_knots, false);
	return_val_if_fail(Mesh.nurbs_curve_groups->curve_points, false);
	return_val_if_fail(Mesh.nurbs_curve_groups->curve_point_weights, false);
	return_val_if_fail(Mesh.nurbs_curve_groups->curve_knots, false);
	return_val_if_fail(Mesh.nurbs_curve_groups->curve_selection, false);
	return_val_if_fail(Mesh.nurbs_curve_groups->materials, false);

	return true;
}

const bool_t validate_bilinear_patches(const mesh& Mesh)
{
	if(!Mesh.bilinear_patches)
		return false;

	return_val_if_fail(validate_points(Mesh), false);
	return_val_if_fail(Mesh.bilinear_patches->patch_points, false);
	return_val_if_fail(Mesh.bilinear_patches->patch_selection, false);
	return_val_if_fail(Mesh.bilinear_patches->patch_materials, false);

	return true;
}

const bool_t validate_bicubic_patches(const mesh& Mesh)
{
	if(!Mesh.bicubic_patches)
		return false;

	return_val_if_fail(validate_points(Mesh), false);
	return_val_if_fail(Mesh.bicubic_patches->patch_points, false);
	return_val_if_fail(Mesh.bicubic_patches->patch_selection, false);
	return_val_if_fail(Mesh.bicubic_patches->patch_materials, false);

	return true;
}

const bool_t validate_nurbs_patches(const mesh& Mesh)
{
	if(!Mesh.nurbs_patches)
		return false;

	return_val_if_fail(validate_points(Mesh), false);
	return_val_if_fail(Mesh.nurbs_patches->patch_first_points, false);
	return_val_if_fail(Mesh.nurbs_patches->patch_u_point_counts, false);
	return_val_if_fail(Mesh.nurbs_patches->patch_v_point_counts, false);
	return_val_if_fail(Mesh.nurbs_patches->patch_u_orders, false);
	return_val_if_fail(Mesh.nurbs_patches->patch_v_orders, false);
	return_val_if_fail(Mesh.nurbs_patches->patch_u_first_knots, false);
	return_val_if_fail(Mesh.nurbs_patches->patch_v_first_knots, false);
	return_val_if_fail(Mesh.nurbs_patches->patch_points, false);
	return_val_if_fail(Mesh.nurbs_patches->patch_point_weights, false);
	return_val_if_fail(Mesh.nurbs_patches->patch_u_knots, false);
	return_val_if_fail(Mesh.nurbs_patches->patch_v_knots, false);
	return_val_if_fail(Mesh.nurbs_patches->patch_selection, false);
	return_val_if_fail(Mesh.nurbs_patches->patch_materials, false);

	return true;
}

const bool_t validate_polyhedra(const mesh& Mesh)
{
	if(!Mesh.polyhedra)
		return false;

	return_val_if_fail(validate_points(Mesh), false);
	return_val_if_fail(Mesh.polyhedra->first_faces, false);
	return_val_if_fail(Mesh.polyhedra->face_counts, false);
	return_val_if_fail(Mesh.polyhedra->types, false);
	return_val_if_fail(Mesh.polyhedra->face_first_loops, false);
	return_val_if_fail(Mesh.polyhedra->face_loop_counts, false);
	return_val_if_fail(Mesh.polyhedra->face_selection, false);
	return_val_if_fail(Mesh.polyhedra->face_materials, false);
	return_val_if_fail(Mesh.polyhedra->loop_first_edges, false);
	return_val_if_fail(Mesh.polyhedra->edge_points, false);
	return_val_if_fail(Mesh.polyhedra->clockwise_edges, false);
	return_val_if_fail(Mesh.polyhedra->edge_selection, false);

	return_val_if_fail(Mesh.polyhedra->first_faces->size() == Mesh.polyhedra->face_counts->size(), false);
	return_val_if_fail(Mesh.polyhedra->first_faces->size() == Mesh.polyhedra->types->size(), false);

	return_val_if_fail(Mesh.polyhedra->face_first_loops->size() == Mesh.polyhedra->face_loop_counts->size(), false);
	return_val_if_fail(Mesh.polyhedra->face_first_loops->size() == Mesh.polyhedra->face_selection->size(), false);
	return_val_if_fail(Mesh.polyhedra->face_first_loops->size() == Mesh.polyhedra->face_materials->size(), false);

	return_val_if_fail(Mesh.polyhedra->edge_points->size() == Mesh.polyhedra->clockwise_edges->size(), false);
	return_val_if_fail(Mesh.polyhedra->edge_points->size() == Mesh.polyhedra->edge_selection->size(), false);

	// Check for infinite loops
	// TODO:  Catch segmentation faults
	const mesh::indices_t& loop_first_edges = *Mesh.polyhedra->loop_first_edges;
	const mesh::indices_t& edge_points = *Mesh.polyhedra->edge_points;
	const mesh::indices_t& clockwise_edges = *Mesh.polyhedra->clockwise_edges;

	const uint_t loop_begin = 0;
	const uint_t loop_end = loop_begin + loop_first_edges.size();
	for(uint_t loop = loop_begin; loop != loop_end; ++loop)
	{
		const uint_t first_edge = loop_first_edges[loop];
		uint_t edge_slow = first_edge;
		uint_t edge_fast = first_edge;
		uint_t cycle_count = 0;
		while(true)
		{
			edge_slow = clockwise_edges[edge_slow];
			edge_fast = clockwise_edges[clockwise_edges[edge_fast]];

			if(edge_slow == edge_fast)
				++cycle_count;

			if(cycle_count > 2)
			{
				log() << error << "infinite loop at loop index " << loop << std::endl;
				return false;
			}

			if(edge_slow == first_edge)
				break;
		}
	}
	return true;
}

const bool_t validate_blobbies(const mesh& Mesh)
{
	if(!Mesh.blobbies)
		return false;

	return_val_if_fail(Mesh.blobbies->first_primitives, false);
	return_val_if_fail(Mesh.blobbies->primitive_counts, false);
	return_val_if_fail(Mesh.blobbies->first_operators, false);
	return_val_if_fail(Mesh.blobbies->operator_counts, false);
	return_val_if_fail(Mesh.blobbies->primitives, false);
	return_val_if_fail(Mesh.blobbies->primitive_first_floats, false);
	return_val_if_fail(Mesh.blobbies->primitive_float_counts, false);
	return_val_if_fail(Mesh.blobbies->operators, false);
	return_val_if_fail(Mesh.blobbies->operator_first_operands, false);
	return_val_if_fail(Mesh.blobbies->floats, false);
	return_val_if_fail(Mesh.blobbies->operands, false);
	return_val_if_fail(Mesh.blobbies->materials, false);

	return true;
}

const bool_t is_sds(const mesh& Mesh)
{
	if(!validate_polyhedra(Mesh))
		return false;

	const mesh::polyhedra_t::types_t& types = *Mesh.polyhedra->types;
	for(uint_t type = 0; type != types.size(); ++type)
	{
		if(types[type] == mesh::polyhedra_t::CATMULL_CLARK)
			return true;
	}
	return false;
}

} // namespace k3d
