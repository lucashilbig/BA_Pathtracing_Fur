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

		/**
		A specific face for fur geometry. Adds a vector of vertices which represents the fur fiber as a Line Strap.
		The first vertice is the start point which is located at the center of the triangle.
		In Adition each vertice has a radius which determines the width of the fur fiber at that point.
		The radius of the last vertice is always 0.
		*/
		struct furFace : KIRK::Mesh::face	//TODO: Fur Geometry implementation
		{
			std::vector<glm::vec3> fiber_position;
			std::vector<unsigned float> fiber_radius;
		};

		virtual ~Mesh() {}

		/**
		Creates a separate mesh for each material contained in this mesh. Won't modify the mesh.
		@return A vector of meshes containing one or more separated meshes with one material each.
		*/
		std::vector<Mesh> separateByMaterial() const;

		/**
		Changes the faces in m_faces from Mesh::face to Mesh::furFace therefor adding fur fibers to the faces.
		This will add fiber_position and fiber_radius at the center of the face pointing in the normals direction.
		@param num_fiber_verts amount of vertices if which the fiber consists. Also determines the length of the fiber since the distance between each vertice is the same.
		@param fiber_radius is the radius of the fur fiber. Has to be > 0.
		*/
		void changeToFurFaces(unsigned int num_fiber_verts, float fiber_radius);

        std::vector<vertex> m_vertices;			//!< All vertices contained in this mesh.
        std::vector<unsigned int> m_indices;	//!< The indices in this mesh, sorted like Triangles as they should be for OpenGL rendering.
        std::vector<std::shared_ptr<Material>> m_materials;		//!< All materials this mesh uses.
        std::vector<std::shared_ptr<face>> m_faces;				//!< The triangles this mesh consists of.

    };
}

#endif //KIRK_MESH_H
