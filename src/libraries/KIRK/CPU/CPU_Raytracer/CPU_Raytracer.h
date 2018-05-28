//
// Created by Maximilian Luzius on 11/05/16.
//

#ifndef RT_CVK_RT_CPU_RAYTRACER_H
#define RT_CVK_RT_CPU_RAYTRACER_H

#include "KIRK/Common/Texture.h"
#include "KIRK/Common/Camera.h"
#include "KIRK/Common/SceneGraph.h"
#include "KIRK/CPU/CPU_Raytracer/CPU_Scene.h"

namespace KIRK {
namespace CPU {

class CPU_Raytracer
{

public:

    CPU_Raytracer()
    {};

    virtual ~CPU_Raytracer()
    {
    };

    virtual void init(std::shared_ptr<KIRK::CPU::Scene> cpuscene)
    {
		setCPUScene(cpuscene);
    };

	void setCPUScene(std::shared_ptr<KIRK::CPU::Scene> cpuscene)
	{
		m_cpuscene = cpuscene;
	}

    /**
     Renders the image out to the currently active kirk camera and copies it over to a texture. In the process, the target texture will be resized accordingly.
     @param texture The target texture to render to.
     */
    void renderToTexture(std::shared_ptr<KIRK::Texture> render_texture)
    {
        m_render_texture = render_texture;
        m_cpuscene->updateFromSceneGraph();
        m_cpuscene->getActiveCamera().setResolution(m_render_texture->getSize());
        m_cpuscene->getActiveCamera().applyParameters();
        render();
    }


	/**
	* Get the scene. DO NOT CHACNGE IT. It should totally be a const ref, but is not
	* supported by all the functions that need the scene yet.
	* @return The scene Object which will be rendered by the raytracer.
	*/
	KIRK::CPU::Scene &getScene() const
	{
		return *m_cpuscene;
	}

	/**
	 * \brief Sets the bounce count aka the depth
	 * \param depth 
	 */
	void setDepth(int depth) { m_depth = depth; }

protected:

    virtual void render() = 0;

    std::shared_ptr<KIRK::Texture> m_render_texture;
	std::shared_ptr<KIRK::CPU::Scene> m_cpuscene;

    int m_depth = 8;


};
};
}

#endif //RT_CVK_RT_CPU_RAYTRACER_H
