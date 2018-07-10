#ifndef __CVK_RT_SCENE_H
#define __CVK_RT_SCENE_H

#include "KIRK/Common/SceneGraph.h"
#include "KIRK/Common/SceneNode.h"
#include "KIRK/Common/Triangle.h"

//Ray tracing needs a copy of the scene, because...:
//    - only world coordinates (to speed up ray tracing)
//    - smaller objects without VBOs (for better cache coherence)
//    - to substitute "geometry" with many triangles by a list of triangle objects (without VBOs)

//////////////////////////////////////////////////////////////////////
//////////
//////////      RING INCLUDE PREVENTIONS
//////////
//////////////////////////////////////////////////////////////////////
namespace KIRK {
class SceneGraph;
}

namespace KIRK {
    namespace CPU{
        class CPU_DataStructure;
    }
}


//////////////////////////////////////////////////////////////////////
//////////
//////////      SCENE
//////////
//////////////////////////////////////////////////////////////////////

namespace KIRK {
namespace CPU {

class Scene
{
public:
    /**
     * The Constructor for KIRK::SceneGraph
     * @param sceneGraph give a shared pointer to the sceneGraph you want to render
     * @param datastructure create a datastructur of the type you want to use and transfer ownership to the scene
     * */
    Scene(std::weak_ptr <KIRK::SceneGraph> sceneGraph, std::unique_ptr <KIRK::CPU::CPU_DataStructure> datastructure);

    /** non copy and non moveable (add copy constructor yourself if you need one :p)*/
    Scene(Scene &other) = delete;
    Scene(Scene &&other) = delete;

    /** Destructor */
    ~Scene();

    /**
     * Apply changes from the scene graph to the local data structures.
     * @param If updateAll = true all Nodes of the SceneGraph will be flattened again and the
     * datascructure will be rebuild. if false only camera and environment will be updated
     */
    void updateFromSceneGraph(bool updateAll = false);

    /**
     * Change the sceneGraph this scene is build from
     * @param sceneGraph Create the Scene from this SceneGraph.
     */
    void setSceneGraph(std::weak_ptr <KIRK::SceneGraph> sceneGraph);

    /**
     * Set the datastructure by taking ownership of the passed unique_ptr of any CPUDataStructure.
     * @param datastructure Any CPU_Datastructure. After calling the method, this parameter will be an invalid unique_ptr.
     * @param buildNew decides if the datastructure should be autmonatically filled. If you set this to false, you will not be able
     *        to render your scene until you call updateFromSceneGraph().
     */
    void setDataStructure(std::unique_ptr <KIRK::CPU::CPU_DataStructure> datastructure, bool buildNew = true);

	/*
	   When set to true, all fur fibers will be converted to cylinder intersection objects rather than to triangle objects.
	*/
	void setFiberAsCylinder(bool fiberAsCylinder) { m_fiberAsCylinder = fiberAsCylinder; }

    /**
     * [DEPRICATED]
    * Set the datastructure of the scene by creating a copy of the given datastructure.
    * @param datastructure Any CPU_DataStructure.
    * @attention Better use setDataStructure. This function acts as an easy to use test for older implementations.
    */
    template<typename T>
    void copyDataStructure(T *datastructure)
    {
        //make a unique Pointer of a new object using the CPU_DataStructures copy constructor.
        m_datastructure = std::make_unique<T>(&datastructure);
    }

    /**
     * @brief get a reference to the datastructure. DO NOT CHACNGE IT. It should totally be a const ref, but is not
     * supported by all the functions that need the datastructure yet.
     * @return A reference to the datastructure containing the scene geometry.
     */
    CPU_DataStructure &getDataStructure() const
    {
        return *m_datastructure;
    }

    /**
     * Getter for all the objects stored in the scene
     * @return A const reference to this Scene's Objects vector. It chall not be changed outside!
     */
    const std::vector<Object *> &getSceneObjects() const;

    /**
     * Get the boundig box of the scene
     * @return A constant glm::vec3[2]. It shall not be changed outside!
     */
    const glm::vec3 *getBounds() const;

    /**
     * Get the active camera. DO NOT CHACNGE IT. It should totally be a const ref, but is not
     * supported by all the functions that need the camera yet.
     * @return The camera Object the scene should be rendered from
     */
    KIRK::Camera& getActiveCamera() {return *m_camera;}

    /**
     * Get a reference to a std::vector of all lights. DO NOT CHACNGE IT. It should totally be a const ref, but is not
     * supported by all the functions that need the camera yet.
     * @return a const reference to a std::vector of all the lightsources in the scene
     */
    std::vector<std::unique_ptr<KIRK::Light>>& getLights() {return m_lights;}

    /**
     * Get the environment of the scene. DO NOT CHACNGE IT. It should totally be a const ref, but is not
     * supprted by all the functions that need the environment yet.
     * @return The camera Object the scene should be rendered from
     */
    KIRK::Environment& getEnvironment() {return *m_environment;}

private:
    /**
     Makes CVK::Triangles from KIRK::Meshes and transforms them.
     @param sceneNode The current SceneNode to flatten.
     @param base_transform The current SceneNode's world transform.
     */
    void flattenNode(std::shared_ptr <KIRK::SceneNode> sceneNode, glm::mat4 base_transform);

    /**
    Builds the datastructure by calling its addBaseDataStructure method with this scene as parameter.
    */
    void buildDatastructure();

    /**
     Computes and sets the whole scene's bounding box.
     */
    void computeBounds();

	/**
	Converts a fur fiber, which consists of cones, to KIRK::Triangles.
	@param fiber The fur fiber struct from the Mesh which will be converted to triangles.
	@param mesh_transform The current meshes transformation matrix.
	@param resolution Determines the amount of triangles which will be build from the cone geometry of the fiber.
	@return Vector of all Triangles that were build from the fur fiber.
	*/
	std::vector<Triangle *> fiberToTriangles(KIRK::Mesh::furFiber fiber, glm::mat4 mesh_transform, unsigned int resolution);

    glm::vec3 m_bound[2]; //!< Bounding limits of the whole scene.
    std::vector<Object *> m_scene_objects; //!< all Objects of the scene
    std::weak_ptr <KIRK::SceneGraph> m_sceneGraph; //!< Remember what scene you were created from. "Scene, I am your father!"
    std::unique_ptr <KIRK::CPU::CPU_DataStructure> m_datastructure; //!< The datastructure to be created from this scene.

    std::unique_ptr<KIRK::Camera> m_camera; //!< the camera we want to render from
    std::vector <std::unique_ptr<KIRK::Light>> m_lights; //!< the lights of the scene
    std::vector <std::shared_ptr<KIRK::Material>> m_materials; //!< the materials of the scene we dont need to copy so share them
    std::unique_ptr<KIRK::Environment> m_environment; //! the environment of the scene
	bool m_fiberAsCylinder = true;//! (default)true: fibers are converted to cylinders. false : fibers are converted to triangles. 
};
}
}

//Those need to be included down here to prevent ring-includes
#include "KIRK/Common/SceneGraph.h"
#include "KIRK/CPU/CPU_Datastructures/CPU_DataStructure.h"

#endif // __CVK_RT_SCENE_H

