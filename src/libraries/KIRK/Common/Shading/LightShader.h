#ifndef KIRK_ENVIRONMENTSHADER_H
#define KIRK_ENVIRONMENTSHADER_H

#include <KIRK/Common/Shading/Shader.h>
#include "KIRK/Common/Ray.h"
#include "ShaderFactory.h"
#include "KIRK/CPU/CPU_Raytracer/CPU_PathTracer.h"

namespace KIRK {

	class LightShader : public Shader
	{
	public:
        LightShader(const std::string name) : Shader(name){}
		void shade(KIRK::CPU::PathTracer& pathtracer, const KIRK::Intersection& hit, KIRK::CPU::Bounce& resultBounce, KIRK::Ray& resultRay) override;
	};

    const ShaderRegistrator <LightShader> lightShaderRegistrator("LightShader");

	inline void LightShader::shade(KIRK::CPU::PathTracer& pathtracer, const KIRK::Intersection& hit, KIRK::CPU::Bounce& resultBounce, KIRK::Ray& resultRay)
	{
		auto light = pathtracer.getScene().getLights()[hit.m_barycentric_coord.y].get();
		resultBounce.color += glm::min(glm::vec4(light->sampleLightSource(hit.m_ray.m_direction, hit.m_ray.m_origin, hit.m_lambda) * resultBounce.radiance, 1.0f), 1.0f);
		resultBounce.radiance = glm::vec3(0);
	}

}

#endif
