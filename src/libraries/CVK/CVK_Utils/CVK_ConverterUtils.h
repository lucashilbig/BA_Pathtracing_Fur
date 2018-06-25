#ifndef KIRK_CVK_CONVERTER_UTILS_H
#define KIRK_CVK_CONVERTER_UTILS_H

#include <memory>

#include "KIRK/Common/Camera.h"
#include "KIRK/Common/Material.h"
#include "KIRK/Common/SceneGraph.h"
#include <CVK/CVK_2/CVK_Node.h>
#include <CVK/CVK_2/CVK_Camera.h>
#include <CVK/CVK_2/CVK_Material.h>
#include <CVK/CVK_2/CVK_FreeFlight.h>
#include <CVK/CVK_2/CVK_Trackball.h>
#include <CVK/CVK_2/CVK_Cone.h>
#include <CVK/CVK_2/CVK_Sphere.h>
#include "KIRK/Common/Ray.h"

namespace CVK
{
	/**
	This is a wrapper for CVK-KIRK-Camera interaction. We can link a CVK::Camera to a KIRK::Camera to automate their synchronization of values/settings.
	*/
    class CVKCameraSynchronizer
    {
    public:
        CVKCameraSynchronizer(std::weak_ptr<KIRK::SceneGraph> scene, std::shared_ptr<CVK::Camera> source);

		/**
		Updates the CVK::Camera and updates all values of the KIRK::Camera accordingly.
		*/
        void update(GLFWwindow* window);

		/**
		* Sets the current used cvk_camera
		* @param scene the KIRK::SceneGraph of the current scene
		* @param camera The new CVK::Camera we want to use
		*/
		void setCVKCamera(std::weak_ptr<KIRK::SceneGraph> scene, std::shared_ptr<CVK::Camera> camera);

		std::shared_ptr<CVK::Camera> getCVKCamera() const;

    private:
        KIRK::Camera* m_kirk_camera;					//!< The KIRK::Camera contained in the sceneGraph
        std::shared_ptr<CVK::Camera> m_cvk_camera;		//!< The CVK::Camera. Can be a Trackball or something.
    };

    class SceneToCVK
    {
    public:
        /**
         Converts the SceneGraph. Available types are CVK_RT::Scene, CVK_RT::GPU_Scene and CVK::Node.
         @param The base KIRK::SceneGraph
		 @param [Optional]showLightGeometry If set to true a sphere will be rendered at every Light position
         @returns The converted Scene.
         */

        static std::shared_ptr<CVK::Node> exportScene(std::shared_ptr<KIRK::SceneGraph> sceneGraph, bool showLightGeometry = false);
    private:
        /**
         Converts the KIRK::SceneNode recursively to a vector of CVK::Node pointers. That is needed, as KIRK::SceneNodes
         sometimes have to be split because they contain Meshes with more than one Material.
         @param sceneNode Base node to convert the sub-tree of.
         @param abs_base_transform The current parent node transform.
         @return an std::vector<CVK::Node*> of nodes containing either nothing or a with exactly one Material each.
         */
        static std::vector<std::shared_ptr<CVK::Node>> toCVKNode(const std::shared_ptr<KIRK::SceneNode> sceneNode, bool showLightGeometry);
    };
}

#endif //KIRK_CVK_CONVERTER_UTILS_H
