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
	\author Timothy M. Shead (tshead@k-3d.com)
	\author Romain Behar (romainbehar@yahoo.com)
	\author Bart Janssens (bart.janssens@lid.kviv.be)
*/

#include "detail.h"

#include <k3dsdk/imaterial.h>
#include <k3dsdk/attribute_array_copier.h>
#include <k3dsdk/shared_pointer.h>

namespace module
{

namespace blobby
{

namespace detail
{

void merge(const mesh_collection& Inputs, k3d::imaterial* const Material, const k3d::mesh::blobbies_t::operator_type Operator, const k3d::bool_t VariableArguments, k3d::mesh& Output)
{
	// Collect all of the varying and vertex arrays to be merged ...
	k3d::mesh::attribute_arrays_t::attribute_arrays_collection source_varying_data;
	k3d::mesh::attribute_arrays_t::attribute_arrays_collection source_vertex_data;
	for(mesh_collection::const_iterator mesh = Inputs.begin(); mesh != Inputs.end(); ++mesh)
	{
		source_varying_data.push_back(&(**mesh).blobbies->varying_data);
		source_vertex_data.push_back(&(**mesh).blobbies->vertex_data);
	}

	// Setup the initial state of the output mesh ...
	k3d::mesh::blobbies_t& target_blobbies = *k3d::make_unique(Output.blobbies);
	k3d::mesh::indices_t& target_first_primitives = *k3d::make_unique(target_blobbies.first_primitives);
	k3d::mesh::counts_t& target_primitive_counts = *k3d::make_unique(target_blobbies.primitive_counts);
	k3d::mesh::indices_t& target_first_operators = *k3d::make_unique(target_blobbies.first_operators);
	k3d::mesh::counts_t& target_operator_counts = *k3d::make_unique(target_blobbies.operator_counts);
	k3d::mesh::materials_t& target_materials = *k3d::make_unique(target_blobbies.materials);
	k3d::mesh::blobbies_t::primitives_t& target_primitives = *k3d::make_unique(target_blobbies.primitives);
	k3d::mesh::indices_t& target_primitive_first_floats = *k3d::make_unique(target_blobbies.primitive_first_floats);
	k3d::mesh::counts_t& target_primitive_float_counts = *k3d::make_unique(target_blobbies.primitive_float_counts);
	k3d::mesh::attribute_arrays_t& target_varying_data = target_blobbies.varying_data;
	k3d::mesh::attribute_arrays_t& target_vertex_data = target_blobbies.vertex_data;
	k3d::mesh::blobbies_t::operators_t& target_operators = *k3d::make_unique(target_blobbies.operators);
	k3d::mesh::indices_t& target_operator_first_operands = *k3d::make_unique(target_blobbies.operator_first_operands);
	k3d::mesh::counts_t& target_operator_operand_counts = *k3d::make_unique(target_blobbies.operator_operand_counts);
	k3d::mesh::blobbies_t::floats_t& target_floats = *k3d::make_unique(target_blobbies.floats);
	k3d::mesh::blobbies_t::operands_t& target_operands = *k3d::make_unique(target_blobbies.operands);

	target_varying_data = k3d::attribute_arrays::clone_types(source_varying_data);
	target_vertex_data = k3d::attribute_arrays::clone_types(source_vertex_data);

	target_first_primitives.push_back(0);
	target_primitive_counts.push_back(0);
	target_first_operators.push_back(0);
	target_operator_counts.push_back(0);
	target_materials.push_back(Material);

	k3d::mesh::indices_t blobby_roots;

	// Iterate over each input mesh ...
	for(mesh_collection::const_iterator mesh = Inputs.begin(); mesh != Inputs.end(); ++mesh)
	{
		const k3d::mesh::indices_t& source_first_primitives = *(*mesh)->blobbies->first_primitives;
		const k3d::mesh::counts_t& source_primitive_counts = *(*mesh)->blobbies->primitive_counts;
		const k3d::mesh::indices_t& source_first_operators = *(*mesh)->blobbies->first_operators;
		const k3d::mesh::counts_t& source_operator_counts = *(*mesh)->blobbies->operator_counts;
		const k3d::mesh::materials_t& source_materials = *(*mesh)->blobbies->materials;
		const k3d::mesh::blobbies_t::primitives_t& source_primitives = *(*mesh)->blobbies->primitives;
		const k3d::mesh::indices_t& source_primitive_first_floats = *(*mesh)->blobbies->primitive_first_floats;
		const k3d::mesh::counts_t& source_primitive_float_counts = *(*mesh)->blobbies->primitive_float_counts;
		const k3d::mesh::attribute_arrays_t& source_varying_data = (*mesh)->blobbies->varying_data;
		const k3d::mesh::attribute_arrays_t& source_vertex_data = (*mesh)->blobbies->vertex_data;
		const k3d::mesh::blobbies_t::operators_t& source_operators = *(*mesh)->blobbies->operators;
		const k3d::mesh::indices_t& source_operator_first_operands = *(*mesh)->blobbies->operator_first_operands;
		const k3d::mesh::counts_t& source_operator_operand_counts = *(*mesh)->blobbies->operator_operand_counts;
		const k3d::mesh::blobbies_t::floats_t& source_floats = *(*mesh)->blobbies->floats;
		const k3d::mesh::blobbies_t::operands_t& source_operands = *(*mesh)->blobbies->operands;

		k3d::attribute_array_copier varying_copier(source_varying_data, target_varying_data, k3d::attribute_array_copier::copy_subset());
		k3d::attribute_array_copier vertex_copier(source_vertex_data, target_vertex_data, k3d::attribute_array_copier::copy_subset());

		const k3d::uint_t blobby_begin = 0;
		const k3d::uint_t blobby_end = blobby_begin + source_first_primitives.size();
		for(k3d::uint_t blobby = blobby_begin; blobby != blobby_end; ++blobby)
		{
			const k3d::uint_t operand_offset = target_primitive_counts[0] + target_operator_counts[0];

			const k3d::uint_t source_primitive_count = source_primitive_counts[blobby];
			const k3d::uint_t source_primitives_begin = source_first_primitives[blobby];
			const k3d::uint_t source_primitives_end = source_primitives_begin + source_primitive_count;
			for(k3d::uint_t source_primitive = source_primitives_begin; source_primitive != source_primitives_end; ++source_primitive)
			{
				target_primitives.push_back(source_primitives[source_primitive]);
				target_primitive_first_floats.push_back(target_floats.size());
				target_primitive_float_counts.push_back(source_primitive_float_counts[source_primitive]);
				varying_copier.push_back(source_primitive);
				vertex_copier.push_back(source_primitive);

				const k3d::uint_t source_floats_begin = source_primitive_first_floats[source_primitive];
				const k3d::uint_t source_floats_end = source_floats_begin + source_primitive_float_counts[source_primitive];
				for(k3d::uint_t source_float = source_floats_begin; source_float != source_floats_end; ++source_float)
					target_floats.push_back(source_floats[source_float]);
			}
			target_primitive_counts[0] += source_primitive_count;

			const k3d::uint_t source_operator_count = source_operator_counts[blobby];
			const k3d::uint_t source_operators_begin = source_first_operators[blobby];
			const k3d::uint_t source_operators_end = source_operators_begin + source_operator_count;
			for(k3d::uint_t source_operator = source_operators_begin; source_operator != source_operators_end; ++source_operator)
			{
				target_operators.push_back(source_operators[source_operator]);
				target_operator_first_operands.push_back(target_operands.size());
				target_operator_operand_counts.push_back(source_operator_operand_counts[source_operator]);

				const k3d::uint_t source_operands_begin = source_operator_first_operands[source_operator];
				const k3d::uint_t source_operands_end = source_operands_begin + source_operator_operand_counts[source_operator];
				for(k3d::uint_t source_operand = source_operands_begin; source_operand != source_operands_end; ++source_operand)
					target_operands.push_back(source_operands[source_operand] + operand_offset);
			}
			target_operator_counts[0] += source_operator_count;

			if(source_primitive_count || source_operator_count)
				blobby_roots.push_back(target_primitive_counts[0] + target_operator_counts[0] - 1);
		}
	}	

	if(!blobby_roots.empty())
	{
		target_operator_counts[0] += 1;
		target_operators.push_back(Operator);
		target_operator_first_operands.push_back(target_operands.size());

		if(VariableArguments)
		{
			target_operator_operand_counts.push_back(blobby_roots.size() + 1);
			target_operands.push_back(blobby_roots.size());
		}
		else
		{
			target_operator_operand_counts.push_back(blobby_roots.size());
		}

		target_operands.insert(target_operands.end(), blobby_roots.begin(), blobby_roots.end());
	}
}

} // namespace detail

} // namespace blobby

} // namespace module
