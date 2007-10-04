// K-3D
// Copyright (c) 1995-2005, Timothy M. Shead
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

#include "cached_triangulation.h"

namespace module
{

namespace development
{

void cached_triangulation::on_execute(const k3d::mesh& Mesh)
{
	m_input_points = Mesh.points;
	if (m_points.empty())
	{
		m_point_links.resize(Mesh.points->size());
		k3d::triangulator::process(Mesh);
	}
	else
	{
		if (m_affected_indices.empty())
		{
			for (k3d::uint_t index = 0; index != m_point_links.size(); ++index)
			{
				k3d::mesh::indices_t& linked_points = m_point_links[index];
				for (k3d::uint_t i = 0; i != linked_points.size(); ++i)
				{
					m_points[linked_points[i]] = m_input_points->at(index);
				} 
			}
		}
		else
		{
			for (k3d::uint_t index = 0; index != m_affected_indices.size(); ++index)
			{
				k3d::mesh::indices_t& linked_points = m_point_links[m_affected_indices[index]];
				for (k3d::uint_t i = 0; i != linked_points.size(); ++i)
				{
					m_points[linked_points[i]] = m_input_points->at(m_affected_indices[index]);
				} 
			}
		}
	}
	m_affected_indices.clear();
}

void cached_triangulation::start_face(const k3d::uint_t Face)
{
	m_point_map.clear();
	m_face_starts.push_back(m_indices.size());
	m_face_points.push_back(k3d::mesh::indices_t());
}

void cached_triangulation::add_vertex(const k3d::point3& Coordinates, k3d::uint_t Vertices[4], k3d::double_t Weights[4], k3d::uint_t& NewVertex)
{
	k3d::log() << debug << "New vertex in triangulated painter: " << Coordinates << ", " << NewVertex << std::endl;
}

void cached_triangulation::add_triangle(const k3d::uint_t Point1, const k3d::uint_t Point2, const k3d::uint_t Point3)
{
	typedef std::pair<point_map_t::iterator, bool> result_t;
	// Create point copies for this face, if they don't exist already
	result_t r1 = m_point_map.insert(std::make_pair(Point1, m_points.size()));
	if (r1.second)
	{
		m_point_links[Point1].push_back(m_points.size());
		m_face_points.back().push_back(m_points.size());
		m_points.push_back(m_input_points->at(Point1));
	}
	result_t r2 = m_point_map.insert(std::make_pair(Point2, m_points.size()));
	if (r2.second)
	{
		m_point_links[Point2].push_back(m_points.size());
		m_face_points.back().push_back(m_points.size());
		m_points.push_back(m_input_points->at(Point2));
	}
	result_t r3 = m_point_map.insert(std::make_pair(Point3, m_points.size()));
	if (r3.second)
	{
		m_point_links[Point3].push_back(m_points.size());
		m_face_points.back().push_back(m_points.size());
		m_points.push_back(m_input_points->at(Point3));
		
	}
		
	// Store corner indices
	m_indices.push_back(r1.first->second);
	m_indices.push_back(r2.first->second);
	m_indices.push_back(r3.first->second);
}

}

}
