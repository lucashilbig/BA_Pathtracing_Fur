#ifndef __KIRK_MARSCHNERHAIRSHADER_H
#define __KIRK_MARSCHNERHAIRSHADER_H

#include <KIRK/Common/Shading/Shader.h>
#include "KIRK/Common/Shading/Bsdf.h"
#include "KIRK/Common/Ray.h"
#include "KIRK/Common/Shading/ShaderFactory.h"
#include "KIRK/CPU/CPU_Raytracer/CPU_PathTracer.h"
#include "KIRK/CPU/CPU_Raytracer/CPU_LightBounce.h"

namespace KIRK {

	class MarschnerHairShader : public Shader
	{
	public:
		MarschnerHairShader(const std::string name) : Shader(name)
		{
			std::random_device rd;
			m_gen = std::mt19937(rd());
		}

		void shade(KIRK::CPU::PathTracer& pathtracer, const KIRK::Intersection& hit, KIRK::CPU::Bounce& resultBounce, KIRK::Ray& resultRay) override;
	private:
		mutable std::mt19937 m_gen;
		mutable std::uniform_real_distribution<> m_dist{ 0, 1 };
		Color::RGBA calcDirectLight(KIRK::CPU::PathTracer& pathtracer, const KIRK::Intersection& hit, const glm::vec2 sample) const;
	};

	const ShaderRegistrator <MarschnerHairShader> MarschnerHairShaderRegistrator("MarschnerHairShader");

	inline void MarschnerHairShader::shade(KIRK::CPU::PathTracer& pathtracer, const KIRK::Intersection& hit, KIRK::CPU::Bounce& resultBounce, KIRK::Ray& resultRay)
	{
		glm::vec3 counter_ray = -normalize(resultRay.m_direction);
		Material* material = hit.m_object->getMaterial();
		std::shared_ptr<BSDF> bsdf = material->m_bsdf;
		glm::vec2 sample(0.f);

		glm::vec3 result_direction; // reflected/refracted ray
		float pdf; // probability distribution function		

		// how much light (r, g, b) gets reflected in this hit point
		glm::vec3 reflectance = bsdf->sample(hit, counter_ray, hit.m_normal, result_direction, pdf, resultBounce.mat_flags, sample); // Sample bsdf
		 
		//offset for new ray, so it doesnt hit the same object directly again
		glm::vec3 offset = 1e-4f * result_direction;

		// If it was not a specular bounce, we optimise the offset for that case 
		if (!(resultBounce.mat_flags & BSDFHelper::MATFLAG_SPECULAR_BOUNCE))
		{
			offset = 1e-4f * hit.m_normal;
			offset = glm::faceforward(-offset, hit.m_normal, result_direction);
		}

		//New Ray to follow
		resultRay = KIRK::Ray(hit.m_location + offset, result_direction);

		//If we have an unfinished TT-/TRT-Path we dont have to calculate the color/radiance yet
		if ((resultBounce.mat_flags & BSDFHelper::MATFLAG_CYLINDER_T_BOUNCE) || (resultBounce.mat_flags & BSDFHelper::MATFLAG_CYLINDER_TR_BOUNCE))
			return;
				

		//Color vectors
		Color::RGBA accumulatedColor(0.0f);
		Color::RGBA ambientLight;
		Color::RGBA directLight;

		//Ambient light
		ambientLight = pathtracer.getScene().getEnvironment().getAmbientLight() * glm::vec4(bsdf->evaluateLight(hit, hit.m_normal, hit.m_normal) * glm::one_over_pi<float>(), 1.0f);
		
		//Calculate direct Light from all light Sources. 
		directLight = calcDirectLight(pathtracer, hit, pathtracer.getSampler().sample2D());

		//Add color from light sources together
		accumulatedColor += directLight * glm::vec4(resultBounce.radiance, 1.0f);
		accumulatedColor += ambientLight * glm::vec4(resultBounce.radiance, 1.0f);
		
		// Scattering Integral of outgoing radiance (Marschner Paper Formel 1)
		if (reflectance == glm::vec3(0.f) || pdf <= 1E-4f || std::max(resultBounce.radiance.x, std::max(resultBounce.radiance.y, resultBounce.radiance.z)) < 0.01f)
			resultBounce.radiance = glm::vec3(0);
		else
			resultBounce.radiance *= reflectance * glm::abs(glm::cos(sample.x));

		resultBounce.color += accumulatedColor;
	}


	inline KIRK::Color::RGBA MarschnerHairShader::calcDirectLight(KIRK::CPU::PathTracer& pathtracer, const KIRK::Intersection& hit, const glm::vec2 sample) const
	{
		Color::RGBA directLight(0.0f);

		if (pathtracer.getScene().getLights().size() == 0)
			return directLight;

		/*******************************************************/
		/******************* LIGHT SAMPLING ********************/
		/*******************************************************/
		int lightIndex = m_dist(m_gen) * pathtracer.getScene().getLights().size();
		auto light = pathtracer.getScene().getLights()[lightIndex].get();

		float attenuation;
		KIRK::Ray hitToLight = light->calcLightdir(hit.m_location, attenuation, true);
		glm::vec3 lightpos = hitToLight.m_origin + hitToLight.m_direction;
		hitToLight.m_origin += 1e-4f * glm::faceforward(hit.m_normal, hitToLight.m_origin - lightpos, hit.m_normal);
		hitToLight.m_direction = glm::normalize(hitToLight.m_direction);

		auto lightColor = light->m_color;

		if (light->m_color.x > 0 || light->m_color.y > 0 || light->m_color.z > 0)
		{
			lightColor *= attenuation
				* glm::vec4(hit.m_object->getMaterial()->m_bsdf->evaluateLight(hit, hitToLight.m_direction, -hit.m_ray.m_direction), 1.0f)
				* std::abs(glm::dot(hitToLight.m_direction, hit.m_normal));

			if (light->m_color.x > 0 || light->m_color.y > 0 || light->m_color.z > 0)
			{
				float t_max = glm::length(lightpos - hitToLight.m_origin);
				bool hitToLight_hasIntersection = pathtracer.getScene().getDataStructure().isIntersection(&hitToLight, t_max);
				if (!hitToLight_hasIntersection)
				{
					// Test for light Intersections
					for (int i = 0; i < pathtracer.getScene().getLights().size(); i++)
					{
						auto lightIntersectionTest = pathtracer.getScene().getLights()[i].get();
						float t;
						if (lightIntersectionTest->isIntersection(hitToLight, t) && (t < t_max))
						{
							hitToLight_hasIntersection = true;
							break;
						}
					}
				}
				lightColor *= !hitToLight_hasIntersection;
			}
			directLight += lightColor;
		}

		return directLight;
	}
}

#endif
