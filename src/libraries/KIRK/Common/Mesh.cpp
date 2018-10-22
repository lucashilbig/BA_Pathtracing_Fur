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

	for (int m = 0; m < indices_per_material.size(); m++)
	{
		vertices_per_material.push_back(0);

		vertex_offset += m > 0 ? vertices_per_material[m - 1] : 0;

		//Reserve needed space.
		separate_meshes[m].m_indices.reserve(indices_per_material[m]);
		separate_meshes[m].m_faces.reserve(indices_per_material[m] / 3);

		for (int f = index_offset; f < index_offset + indices_per_material[m] / 3; f++)
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

		for (int v = vertex_offset; v < vertex_offset + vertices_per_material[m] + 1; v++)
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

void KIRK::Mesh::addFurToFaces(unsigned int fibers_per_face, unsigned int num_fiber_verts, float fiber_radius)
{
	//check if radius is > 0
	if (fiber_radius <= 0)
	{
		LOG_ERRORs() << "addFurToFaces() : Radius has to be > 0";
		return;
	}
	//random value generator
	std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<float> dist(0.0f, std::nextafter(1.0f, DBL_MAX));//std::nextafter so we get the [0,1] interval instead of [0,1)

	//Iterate over every face in the mesh
	for (int f = 0; f < m_faces.size(); f++)
	{
		//create as many furFibers per face as fibers_per_face says
		for (int fiber = 1; fiber <= fibers_per_face; fiber++)
		{
			//create furFiber 
			furFiber furFace;
			//get faces corner points
			glm::vec3 a = m_vertices[m_faces[f].vertex_index_a].position;
			glm::vec3 b = m_vertices[m_faces[f].vertex_index_b].position;
			glm::vec3 c = m_vertices[m_faces[f].vertex_index_c].position;
			//compute random values for position calculation
			float r1 = dist(mt);
			float r2 = dist(mt);
			if (r1 + r2 >= 1) { r1 = 1 - r1; r2 = 1 - r2; }//edge-case. Otherwise point would be outside the triangle
			//compute fiber start position 
			glm::vec3 pos = a + r1 * (b - a) + r2 * (c - a);
			//move start position down, otherwise the bottom cylinder will stick out of the ground
			pos.y -= 0.003f;
			//Add fiber start position
			furFace.fiber_positions.push_back(pos);
			//Add radius for start position
			float radius = fiber_radius;
			furFace.fiber_radius.push_back(radius);
			//offset for the z-axis position
			float offset_z = 0.2f;//0.2

			//add positions for a single fiber
			for (int i = num_fiber_verts; i > 1; i--) {
				//calculate the offset (y-axis) for the new point of the fiber. Using log and i values which get smaller every iteration,
				//the distance between two vertices is getting smaller towards the end (top) of the fiber.
				float offset_y = glm::log((float)i) / 90.0f;
				//compute new position of fiber vertice
				glm::vec3 point = pos + glm::vec3(0.0f, offset_y, 0.06f); // Param b ändert Krümmung des Haares an sich. Param c ändert Neigung aller Haare.
				//calculate the new offset (z-axis)
				offset_z -= (offset_z / ((float)i + 5));
				//decrease radius towards the top of the fiber
				radius -= (radius / ((float)i + 5));
				//Add position and radius to furFace;
				furFace.fiber_positions.push_back(point);
				furFace.fiber_radius.push_back(radius);
				//change pos to new vertice pos
				pos = point;
			}
			//Add last vertice with radius of 0 at the end of the fiber
			//furFace.fiber_positions.push_back(pos + glm::vec3(0.0f, 0.003f, 0.01f));
			furFace.fiber_radius[furFace.fiber_radius.size() - 1] = 0.0005f;

			//Add new furFace to m_furFaces member
			m_furFibers.push_back(furFace);
		}
	}
}

