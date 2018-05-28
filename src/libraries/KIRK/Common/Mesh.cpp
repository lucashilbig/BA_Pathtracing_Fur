#include "Mesh.h"

std::vector<KIRK::Mesh> KIRK::Mesh::separateByMaterial() const 
{
	std::vector<KIRK::Mesh> separate_meshes;

	if (m_materials.size() <= 1) 
	{
		// Don't run if there is already only one material.
		separate_meshes.push_back(*this);
		return separate_meshes;
	}

	std::vector<unsigned int> indices_per_material;
	std::vector<unsigned int> vertices_per_material;

	//Go through faces and determine the border indices for the materials
	int material_index = 0;
	for (Mesh::face face : m_faces) 
	{
		if (face.material_index >= separate_meshes.size())
		{
			//If the material index is >= the mesh count, we increase it and see it as a "new" material mesh zone.
			Mesh mesh;
			separate_meshes.push_back(mesh);
			indices_per_material.push_back(0);
			material_index = separate_meshes.size() - 1;
		}

		//3 indices per face
		indices_per_material[material_index] += 3;
	}

	//For multiple meshes, we need to have a start index and a start vertex.
	unsigned int index_offset = 0;
	int vertex_offset = 0;

	for (int m = 0; m<indices_per_material.size(); m++) 
	{
		vertices_per_material.push_back(0);

		vertex_offset += m>0 ? vertices_per_material[m - 1] : 0;

		//Reserve needed space.
		separate_meshes[m].m_indices.reserve(indices_per_material[m]);
		separate_meshes[m].m_faces.reserve(indices_per_material[m] / 3);

		for (int f = index_offset; f<index_offset + indices_per_material[m] / 3; f++) 
		{
			//Push indices to separated material.
			unsigned int index_a = m_faces[f].vertex_index_a - vertex_offset;
			unsigned int index_b = m_faces[f].vertex_index_b - vertex_offset;
			unsigned int index_c = m_faces[f].vertex_index_c - vertex_offset;

			vertices_per_material[m] = glm::max(vertices_per_material[m], index_a, index_b, index_c);

			separate_meshes[m].m_faces.push_back(m_faces[f]);
			separate_meshes[m].m_indices.push_back(index_a);
			separate_meshes[m].m_indices.push_back(index_b);
			separate_meshes[m].m_indices.push_back(index_c);
		}

		//Reserve needed space.
		separate_meshes[m].m_vertices.reserve(vertices_per_material[m] + 1);

		for (int v = vertex_offset; v<vertex_offset + vertices_per_material[m] + 1; v++) 
		{
			//Push vertices to separated material.
			separate_meshes[m].m_vertices.push_back(m_vertices[v]);
		}

		//Add material
		separate_meshes[m].m_materials.push_back(m_materials[m]);

		//Then increase start index by the current face count.
		index_offset += separate_meshes[m].m_faces.size();
	}

	return separate_meshes;
}
