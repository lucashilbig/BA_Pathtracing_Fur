#ifndef __CVK_NODE_H
#define __CVK_NODE_H

#include "CVK_Defs.h"
#include <string>
#include "CVK_Geometry.h"
#include "CVK_Material.h"
#include "CVK_State.h"

#include <memory>
#include <assimp/Importer.hpp> //<- WHY???? The node shouldn't care less about loading anything at all... -.-
#include <assimp/scene.h>
#include <assimp/postprocess.h>

namespace CVK {

// forward declaration
class ShaderMinimal;

/**
 * @brief Class for Nodes in the scene graph.
 *
 * A Node in the scene graph represents the hierarchical structure
 * of the scene. One node consists of a list of subnodes, representing
 * all children of this node. It consists of a 4x4 model matrix, for
 * the current localization, rotation and scale of the represented
 * object and all children.
 *
 * Furthermore a node may have a connection to a geometry, which can
 * be rendered and a material, which can be used while rendering by
 * the current shader. Additionally a node may have a name, which can
 * be used to search this specific node in the hierarchy.
 *
 * When a node is created, standard values for all members are set,
 * which can be set manually later. An example creation of a node is:
 * CVK::Node node();
 * node.setGeometry( geometry);
 * node.setMaterial( material);
 * scene.addChild( &node);
 *
 * Another way to create a node is to load a scene with ASSIMP. This
 * is done by using the path of the file, which has to be loaded. The
 * loaded scene is then represented in this node and its child nodes. An
 * example for this is:
 * CVK::Node* scene = new CVK::Node("Scene", "/pathToFile.obj");
 *
 * While rendering the scene, the current node accumulates the model
 * matrix to the current model matrix and sends it to the bound shader.
 * Also a shader update is called with this node. If the node has a
 * geometry, then it is rendered. At the end the render function is
 * called in every child, with a copy of the current accumulated model
 * matrix. An example of rendering the scene is:
 * scene->render();
 */
class Node
{
public:
    /**
     * Standard Constructor for Node
     */
    Node();
    /**
     * Constructor for Node with given parameters
     * @param name The name for this node
     */
    Node(std::string name);
    /**
     * Constructor for Node with given parameters
     * @param name The name for this node
     * @param path The path to a model, which is loaded with assimp
     */
    Node(std::string name, std::string path);
    /**
     * Standard Destructor for Node
     */
    ~Node();

    /**
     * Is called by the constructors and initializes the Node. Normally it is not necessary to call it manually.
     */
    void init();
    /**
     * Loads the model given by the path with assimp. Creates a Geometry with all necessary VBOs and VAO.
     * @brief Loads the model with assimp
     * @param path The path to the model file
     */
    void load(std::string path);

    /**
     * Main render call for a Node and all subnodes. Calls render with model matrix.
     * @brief Renders Geometry and all subnodes.
     */
    void render();
    /**
     * Renders the node by using the given model matrix. Calls render method for all subnodes.
     * @brief Renders Geometry and all subnodes with given modelmatrix.
     * @param modelMatrix The local coordinate system of the Node
     */
    void render(glm::mat4 modelMatrix);

    /**
     * Searches the scene Graph for a Node with the given name. Returns the first occurence.
     * @brief Returns first Node with given name
     * @param name The name to search for
     * @return The first found Node as Pointer.
     */
    Node *find(std::string name);

    /**
     * @brief Standard Setter for name
     * @param name The new name of this object
     */
    void setName(std::string name);
    /**
     * @brief Standard Getter for name
     * @return The name of this object as pointer
     */
    std::string *getName();

    /**
     * @brief Standard Setter for geometry
     * @param geometry The new geometry of this object as pointer
     */
    void setGeometry(std::shared_ptr<CVK::Geometry> geometry);
    /**
     * @brief Standard Getter for geometry
     * @return The geometry of this object as pointer
     */
	std::shared_ptr<CVK::Geometry> getGeometry();

    /**
     * Returns if the current Node has any Geometry attached to it
     * @brief Return if Node has Geometry
     * @return true, if the node has Geometry, false otherwise
     */
    bool hasGeometry();

    /**
     * @brief Standard Setter for material
     * @param material The new material of this object as pointer
     */
    void setMaterial(std::shared_ptr<CVK::Material> material);
    /**
     * @brief Standard Getter for material
     * @return The material of this object as pointer
     */
	std::shared_ptr<CVK::Material> getMaterial();
    /**
    * Returns if the current Node has any Material attached to it
    * @brief Return if Node has Material
    * @return true, if the node has Material, false otherwise
    */
    bool hasMaterial();

    /**
     * @brief Standard Setter for model matrix
     * @param modelMatrix The new model matrix of this object
     */
    void setModelMatrix(glm::mat4 modelMatrix);
    /**
     * @brief Standard Getter for model matrix
     * @return The model matrix of this object as pointer
     */
    glm::mat4 *getModelMatrix();

    /**
     * @brief Standard Setter for parent node
     * @param node The new parent node of this object as pointer
     */
    void setParent(Node *node);
    /**
     * @brief Standard Getter for parent node
     * @return The parent node of this object as pointer
     */
    Node *getParent();

    /**
     * @brief Adds given subnode to subnodelist
     * @param node The new subnode of this object as pointer
     */
    void addChild(std::shared_ptr<Node> node);
    /**
     * @brief Standard Getter for list of subnodes
     * @return The list of subnodes of this object as ref
     */
    std::vector<std::shared_ptr<Node>> &getChildren();
private:
    std::string m_name; //!< The name of the Node
    std::shared_ptr<CVK::Geometry> m_geometry; //!< The geometry to render of the Node
	std::shared_ptr<CVK::Material> m_material; //!< The material used in the shader of the Node
    glm::mat4 m_modelMatrix; //!< The model matrix for the local coordinate system of the Node
    Node *m_parent; //!< The parent node of the Node
    std::vector<std::shared_ptr<Node>> m_children; //!< The list of subnodes of the Node
};

};

#endif /* __CVK_NODE_H */

