#ifndef KIRK_SCENEGRAPH_H
#define KIRK_SCENEGRAPH_H

#include <map>
#include <vector>
#include <memory>
#include <stdexcept>

#include <assimp/scene.h>
#include <assimp/Importer.hpp>
#include <assimp/postprocess.h>
#include <assimp/version.h>

#include "Camera.h"
#include "Light.h"
#include "Mesh.h"
#include "Material.h"
#include "SceneNode.h"
#include "Environment.h"
#include <stack>

namespace KIRK
{
    // Predefine some classes to avoid ring-includes
    class SceneNode;
    class NodeDataObject;
	class SceneNodeIterator;

	/**
	This class is an implementation of a raytracing-oriented SceneGraph.
	*/
    class SceneGraph
    {
    public:
		/**
		Creates a SceneGraph depending on what input path it gets.
		@param file_path If empty, create a default scene with default camera, environment and lights. If not empty, load the given file as a sceneGraph.
		*/
        static std::shared_ptr<KIRK::SceneGraph> makeSceneGraph(const std::string &file_path = "");

		/**Checks, whether the given string ends with the given ending.*/
        static bool endsWith(std::string const &fullString, std::string const &ending);



		/**Default constructor. Just add a root node and no default lights, cameras or environment.*/
        SceneGraph();

        ~SceneGraph();

		/**
		Creates a default camera, positioned at (-3.f, 7.f, -9.f), looking at (0, 0, 0) and having (0, 1, 0) as up vector.
		*/
        void createDefaultCamera();

		/**
		Creates some default lights. (Three point lights with slightly differing colors and an ambient light)
		*/
        void createDefaultLight();

		/**
		Creates a default environment having KIRK::Color::SKY_BLUE as color.
		*/
        void createDefaultEnvironment();

		/**
		Iterates over every Mesh and adds fur fiber information to every face in the meshes.
		@param num_fiber_verts amount of vertices if which the fiber consists. Also determines the length of the fiber since the distance between each vertice is the same.
		@param fiber_radius is the radius of the fur fiber. Has to be > 0.
		*/
		void addFurFibersToAllMeshes(unsigned int num_fiber_verts, float fiber_radius);

		/**
		Imports an .obj file as a SceneNode. The resulting SceneNode has a root-node structure where each child of the returned node is a parent node to a mesh node.
		@param file_path The path to the .obj file.
		@return A shared_ptr to the newly created SceneNode sub-root node.
		*/
        std::shared_ptr<KIRK::SceneNode> importObj(const std::string &file_path);

		/**
		@return The top most SceneNode in the node tree.
		*/
        std::shared_ptr<KIRK::SceneNode> getRootNode() const;

		/**
		@return The currently active camera. Set by setActiveCamera;
		*/
        std::shared_ptr<KIRK::Camera> getActiveCamera() const;

        /**
         * @brief sets the new active camera (the camera should belong to a node of this scene graph)
         */
        void setActiveCamera(std::shared_ptr<Camera> cam) { m_active_camera = cam;};

		/**
		@return A pointer to the currently used Environment.
		*/
        Environment* getEnvironment() const;

		SceneNodeIterator begin();

		SceneNodeIterator end();

    private:
		/**
		Loads an aiNode and assigns the loaded data to a given target_node (SceneNode)
		@param scene The loaded base aiScene
		@param source_node The current aiNode
		@param target_node The node to load the data into
		*/
        void assignAiNode(const aiScene* scene, aiNode* source_node, std::shared_ptr<KIRK::SceneNode> target_node, const std::string &scene_path);

		/**
		Loads a mesh from the given aiNode.
		@param scene The loaded base aiScene
		@param source The current aiNode to load the mesh from
		@return The newly created mesh
		*/
        Mesh* createMeshFromAi(const aiScene* scene, aiNode* source, const std::string &scene_path);

		/**
		Creates a SceneNode and sets it as the root of the node tree.
		*/
        void initRoot();

        std::shared_ptr<SceneNode> m_root_node;					//!< The SceneGraph tree root.
        std::unique_ptr<Environment> m_environment;				//!< The currently used environment.

        std::shared_ptr<KIRK::Camera> m_active_camera;						//!< The currently active camera index.
        unsigned int m_node_count = 0;							//!< The number of nodes currently in the tree.

    };

	class SceneNodeIterator
	{
	public:
		SceneNodeIterator(std::shared_ptr<SceneNode> scenegraph)
			: m_current_node(scenegraph)
		{
			m_child_counter.push(1);
		}

		SceneNodeIterator &operator++()
		{
			if(m_current_node->m_children.size() != 0)
			{
				m_current_node = m_current_node->m_children[m_child_counter.top()-1];
				m_child_counter.top()++;
				m_child_counter.push(1);
			} 
			else
			{
				while(m_child_counter.size() > 0 && m_child_counter.top() > m_current_node->m_children.size() && !m_current_node->m_parent.expired())
				{
					m_child_counter.pop();
					m_current_node = m_current_node->m_parent.lock();
				}

				if (m_child_counter.top() <= m_current_node->m_children.size()) {
					m_current_node = m_current_node->m_children[m_child_counter.top() - 1];
					m_child_counter.top()++;
					m_child_counter.push(1);
				} else
				{
					m_current_node = std::shared_ptr<SceneNode>(nullptr);
				}
			}

			return *this;
		}

		bool operator!=(const SceneNodeIterator & other) { 
			return m_current_node != other.m_current_node;
		}

		std::shared_ptr<SceneNode> &operator*()
		{
			return m_current_node;
		}

		operator bool()
		{
			return bool(m_current_node);
		}

	private:
		std::shared_ptr<SceneNode> m_current_node;
		std::stack<unsigned> m_child_counter;
	};
}

#endif //KIRK_SCENEGRAPH_H
