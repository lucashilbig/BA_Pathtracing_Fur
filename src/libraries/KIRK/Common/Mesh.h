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
		A specific fur geometry. Represents a fur fiber as a sequence of cones.
		You can imagine a cone between two consecetive vertex positions with the radiuses of fiber_radius at each vertice.
		@member fiber_positions : Vertex positions of the fur fiber. These are always the center of the cones.
		@member fiber_radius    : Radius of the fur fiber cone at that vertice(same index as fiber_positions).
								  The radius of the last vertice is always 0.				
		*/
		struct furFiber
		{
			std::vector<glm::vec3> fiber_positions;
			std::vector<float> fiber_radius;
		};

		virtual ~Mesh() {}

		/**
		Creates a separate mesh for each material contained in this mesh. Won't modify the mesh.
		@return A vector of meshes containing one or more separated meshes with one material each.
		*/
		std::vector<Mesh> separateByMaterial() const;

		/**
		Adds fur fiber Information for every face in m_faces and saves them in m_furFibers.
		This will add fiber_position and fiber_radius at the center of the face.
		@param num_fiber_verts amount of vertices if which the fiber consists. Also determines the length of the fiber.
		@param fiber_radius is the radius of the fur fiber. Has to be > 0.
		*/
		void addFurToFaces(unsigned int num_fiber_verts, float fiber_radius);


        std::vector<vertex> m_vertices;			//!< All vertices contained in this mesh.
        std::vector<unsigned int> m_indices;	//!< The indices in this mesh, sorted like Triangles as they should be for OpenGL rendering.
        std::vector<std::shared_ptr<Material>> m_materials;		//!< All materials this mesh uses.
        std::vector<face> m_faces;				//!< The triangles this mesh consists of.
		std::vector<furFiber> m_furFibers;		//!< The fur fibers this mesh contains.

    };
}

#endif //KIRK_MESH_H
