#ifndef KIRK_SCENENODE_H
#define KIRK_SCENENODE_H

#include <vector>
#include <memory>
#include <algorithm>
#include <glm/glm.hpp>
#include "KIRK/Utils/Gui/Gui.h"

namespace KIRK
{
	//Need to pre-define SceneNode for use in NodeDataObject.
	class SceneNode;
	class Mesh;
	class Light;

	/**
	Extend this class to be able to put an object into a SceneNode.
	*/
    class NodeDataObject
	{
	public:
		/**
		@return The global transformation matrix of this node.
		*/
		virtual glm::mat4 calculateTransform();

		std::weak_ptr<SceneNode> m_parent_node;	 //!< The node this object is being assigned to.
	};

	/**
	A SceneNode belongs in the SceneGraph tree and can contain exactly one data object of a well-defined type.
	*/
    class SceneNode : public std::enable_shared_from_this<SceneNode>, public GuiElement
    {
		// Wheeee :)
        friend class SceneGraph;
    public:
		/**
		Defines of which type the assigned data object is.
		*/
        enum Type{
            EMPTY,
            MESH,
            LIGHT,
            CAMERA
        };

		///<Empty constructor -> Empty node.
        SceneNode() : m_name("<empty " + std::to_string(counter++) + ">") {};

		void setName(const std::string &name)
		{
			m_name = name;
		}

		/**
		Assigns the given data to this Node. You have to properly specify of which type the assigned data is.
		@param Assign this to this Node.
		@param The type the assigned data belongs to.
		*/
        void assignData(std::shared_ptr<NodeDataObject> data, const Type data_type);

		/**
		Removes a specified child from this Node.
		@param child Child to remove.
		*/
        void removeChild(const std::shared_ptr<SceneNode> child);

		/**
		Attaches a specified child to this Node.
		@param child Child to attach. Parent of this child will be set in this method too.
		*/
        void attachChild(const std::shared_ptr<SceneNode> child);

		/**
		Transforms this node locally by the given transformation matrix.
		*/
		void transformWith(glm::mat4 transform);

		void onGui() override;

		std::string m_name;
        std::weak_ptr<SceneNode> m_parent;						//!< The parent of this Node
        std::vector<std::shared_ptr<SceneNode>> m_children;	//!< All children of this Node
        Type m_data_type = EMPTY;							//!< The type of which this Node's data object is
        std::shared_ptr<NodeDataObject> m_data_object;						//!< The data object contained in this Node

		glm::mat4 m_transform = glm::mat4(1.0f);			//!< The current local transformation matrix
	private:
		static int counter;
    };
}

#include "Mesh.h"
#include "Light.h"

#endif //KIRK_SCENENODE_H
