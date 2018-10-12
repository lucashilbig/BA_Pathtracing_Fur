#ifndef __KIRK_CHIANGHAIRSHADER_H
#define __KIRK_CHIANGHAIRSHADER_H

#include <KIRK/Common/Shading/Shader.h>
#include "KIRK/Common/Shading/Bsdf.h"
#include "KIRK/Common/Ray.h"
#include "KIRK/Common/Shading/ShaderFactory.h"
#include "KIRK/CPU/CPU_Raytracer/CPU_PathTracer.h"
#include "KIRK/CPU/CPU_Raytracer/CPU_LightBounce.h"

namespace KIRK {

	class ChiangHairShader : public Shader
	{
	public:
		ChiangHairShader(const std::string name) : Shader(name)
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

	const ShaderRegistrator <ChiangHairShader> ChiangHairShaderRegistrator("ChiangHairShader");

	inline void ChiangHairShader::shade(KIRK::CPU::PathTracer& pathtracer, const KIRK::Intersection& hit, KIRK::CPU::Bounce& resultBounce, KIRK::Ray& resultRay)
	{
		Material* material = hit.m_object->getMaterial();
		std::shared_ptr<ChiangHairBSDF> bsdf = material->m_hairBSDF;
		//initialise bsdf values
		float offset_h = -1.f + 2 * m_dist(m_gen);//fibers h offset calculated with random sample [0,1). Values of h [-1,1]
		bsdf->init(offset_h, material->m_beta_m, material->m_beta_n, material->m_alpha, material->m_ior, glm::vec3(material->m_diffuse.value));

		glm::vec2 sample = pathtracer.getSampler().sample2D();

		glm::vec3 output_ray; //Outgoing ray 
		float pdf; //probability distribution function. Since it is not used before bsdf sample we used it to pass the offset_h to the bsdf		

		//////////////////////////////////////////////
		// Get random Light Source to sample
		//////////////////////////////////////////////
		if (pathtracer.getScene().getLights().size() == 0)//If we have no Light Source we cant see objects anyways
			return;

		//Select one random light source from all lights
		int lightIndex = m_dist(m_gen) * pathtracer.getScene().getLights().size();
		auto light = pathtracer.getScene().getLights()[lightIndex].get();
		float attenuation;
		//calculate ray from intersection point towards light source
		KIRK::Ray hitToLight = light->calcLightdir(hit.m_location, attenuation, true);
		//move light ray to cylinders (fibers) local space.
		/*  Our v-axis(vector through the center of the cylinder(apexpoint - basepoint)) is u-axis in marschner model
		and we used the local_output_ray to store our ray towards the light source we want to sample*/
		glm::vec3 hitToLight_local = Math::worldToLocal(glm::normalize(hitToLight.m_direction), hit.m_object->getV(), hit.m_object->getU(), hit.m_object->getW());

		// Sample bsdf 
		glm::vec3 scattering = bsdf->Sample_f(hitToLight_local, &output_ray, sample, &pdf);
		
		//////////////////////////////////////////////
		//  Final Color and outgoing path calculation
		//////////////////////////////////////////////		
		//Transform output ray to world space
		output_ray = Math::localToWorld(output_ray, hit.m_object->getU(), hit.m_object->getV(), hit.m_object->getW());
		//New Ray to follow
		resultRay = KIRK::Ray(hit.m_location + 1e-4f * output_ray, output_ray);//offset for new ray, so it doesnt hit the same object directly again

		//Color vectors
		Color::RGBA accumulatedColor(0.0f);
		Color::RGBA ambientLight;
		Color::RGBA directLight;

		//Ambient light
		//ambientLight = pathtracer.getScene().getEnvironment().getAmbientLight() * glm::vec4(bsdf->evaluateLight(hit, hit.m_normal, hit.m_normal) * glm::one_over_pi<float>(), 1.0f);
		//accumulatedColor += ambientLight * glm::vec4(resultBounce.radiance, 1.0f);
		
		//Calculate direct Light from all light Sources. 
		directLight = calcDirectLight(pathtracer, hit, sample);
		accumulatedColor += directLight * glm::vec4(resultBounce.radiance, 1.0f);//add direct light to color
		
		// Scattering Integral of outgoing radiance (Marschner Paper Formel 1)
		if (scattering == glm::vec3(0.f) || pdf <= 1E-4f 
			|| std::max(resultBounce.radiance.x, std::max(resultBounce.radiance.y, resultBounce.radiance.z)) < 0.01f)
			resultBounce.radiance = glm::vec3(0);
		else
			resultBounce.radiance *= scattering * glm::abs(glm::dot(glm::normalize(resultRay.m_direction), hit.m_normal)) / pdf;//Part of Rendering Equation. MarschnerPaper Equation 1
		

		resultBounce.color += accumulatedColor;
	}


	inline KIRK::Color::RGBA ChiangHairShader::calcDirectLight(KIRK::CPU::PathTracer& pathtracer, const KIRK::Intersection& hit, const glm::vec2 sample) const
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
		KIRK::Ray hitToLight = light->calcLightdir(hit.m_location, attenuation, true);//calculates ray from intersection point towards light source
		glm::vec3 lightpos = hitToLight.m_origin + hitToLight.m_direction;
		hitToLight.m_origin += 1e-4f * glm::faceforward(hit.m_normal, hitToLight.m_origin - lightpos, hit.m_normal);
		hitToLight.m_direction = glm::normalize(hitToLight.m_direction);
		
		attenuation = 1;//set attenuation to 1 because we already calculated attenuation factor in the bsdf sample

		auto lightColor = light->m_color;

		if (light->m_color.x > 0 || light->m_color.y > 0 || light->m_color.z > 0)
		{
			lightColor *= attenuation//attenuation factor of the light. Calculated in light->calcLightdir()
				* glm::vec4(hit.m_object->getMaterial()->m_bsdf->evaluateLight(hit, hitToLight.m_direction, -hit.m_ray.m_direction), 1.0f)//light influence at hit point. calc by the bsdf
				* std::abs(glm::dot(hitToLight.m_direction, hit.m_normal));//angle between ray towards light and normal

			if (lightColor.x > 0 || lightColor.y > 0 || lightColor.z > 0)
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
