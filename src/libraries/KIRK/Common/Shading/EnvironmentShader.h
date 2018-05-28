#ifndef KIRK_ENVIRONMENTSHADER_H
#define KIRK_ENVIRONMENTSHADER_H

#include <KIRK/Common/Shading/Shader.h>
#include "KIRK/Common/Ray.h"
#include "ShaderFactory.h"
#include "KIRK/CPU/CPU_Raytracer/CPU_PathTracer.h"

namespace KIRK {

	class EnvironmentShader : public Shader
	{
	public:
        EnvironmentShader(const std::string name) : Shader(name){}
		void shade(KIRK::CPU::PathTracer& pathtracer, const KIRK::Intersection& hit, KIRK::CPU::Bounce& resultBounce, KIRK::Ray& resultRay) override;
	};

    const ShaderRegistrator <EnvironmentShader> environmentShaderRegistrator("EnvironmentShader");

	inline void EnvironmentShader::shade(KIRK::CPU::PathTracer& pathtracer, const KIRK::Intersection& hit, KIRK::CPU::Bounce& resultBounce, KIRK::Ray& resultRay)
	{
		
		
		resultBounce.color += pathtracer.getScene().getEnvironment().getColor(hit.m_ray.m_direction) * glm::vec4(resultBounce.radiance, 1.0f);
		resultBounce.radiance = glm::vec3(0);
	}		
}

#endif
