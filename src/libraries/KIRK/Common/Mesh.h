#ifndef KIRK_MESH_H
#define KIRK_MESH_H

#include "SceneNode.h"
#include "Material.h"

namespace KIRK
{
	/**
	Represents a collection of faces and lists of their vertices, indices and all materials assigned to them.
	*/
    class Mesh : public NodeDataObject
    {
    public:
		/**
		A vertex has just a position, normal and texture coordinate
		*/
        struct vertex
		{
            glm::vec3 position;
            glm::vec3 normal;
            glm::vec2 texcoord;
        };

		/**
		A face is a triangle. It has the indices of three vertices in this mesh and the index for it's material.
		*/
        struct face
		{
            unsigned int vertex_index_a;
            unsigned int vertex_index_b;
            unsigned int vertex_index_c;
            unsigned int material_index = 0;
        };

		virtual ~Mesh() {}

		/**
		Creates a separate mesh for each material contained in this mesh. Won't modify the mesh.
		@return A vector of meshes containing one or more separated meshes with one material each.
		*/
		std::vector<Mesh> separateByMaterial() const;

        std::vector<vertex> m_vertices;			//!< All vertices contained in this mesh.
        std::vector<unsigned int> m_indices;	//!< The indices in this mesh, sorted like Triangles as they should be for OpenGL rendering.
        std::vector<std::shared_ptr<Material>> m_materials;		//!< All materials this mesh uses.
        std::vector<face> m_faces;				//!< The triangles this mesh consists of.

    };
}

#endif //KIRK_MESH_H
