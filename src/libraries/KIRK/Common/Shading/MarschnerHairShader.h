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
		//Cast Object to cylinder
		KIRK::Cylinder *cylinder_obj = dynamic_cast<KIRK::Cylinder*>(hit.m_object);

		//Diameter of the fiber
		float fiber_width;

		//if we have no cylinder as Object we set fiber width to 100
		if (cylinder_obj == NULL)
			fiber_width = 100.f;
		else
			//Calculate mean of the fibers diameter. Diameter is (2 * radius) and convertion from centimeter in micrometer is (* 10000), hence (* 20000)
			fiber_width = 20000 * ((cylinder_obj->getBaseRadius() + cylinder_obj->getApexRadius()) / 2);

		glm::vec3 input_ray_normalized = normalize(resultRay.m_direction);
		Material* material = hit.m_object->getMaterial();
		std::shared_ptr<BSDF> bsdf = material->m_bsdf;
		glm::vec2 sample(0.f);

		glm::vec3 result_direction_r; //Outgoing ray for R-Path
		glm::vec3 result_direction_tt; //Outgoing ray for TT-Path
		glm::vec3 result_direction_trt; //Outgoing ray for TRT-Path
		glm::vec3 result_direction; //Final Outgoing ray after deciding which of the 3 from above we take
		float pdf_r = 1; //probability distribution function for R-Path
		float pdf_tt = 1; //probability distribution function for TT-Path
		float pdf_trt = 1; //probability distribution function for TRT-Path	
		float pdf; //Combined probability distribution function

		//random uniform distribution for lobe alpha shift and beta width. We calculate it here not in the bsdf so we have the same value for every path in the lobe.
		std::uniform_real_distribution<float> dist(5.0f, std::nextafter(10.0f, DBL_MAX));//std::nextafter so we get the [5,10] interval instead of [5,10)
		//We use the sample parameter, because we dont need it in the bsdf and we cant get rid of it cause of inherited method signature.
		sample.x = -1.0f * dist(m_gen);//longitudinal shift. Suggested value from marschner hair paper between -10 and -5 degrees
		sample.y = dist(m_gen);//longitudinal width (stdev.). Suggested value from marschner hair paper between 5 and 10 degrees

		//////////////////////////////////////////////
		// Get random Light Source to sample
		//////////////////////////////////////////////
		if (pathtracer.getScene().getLights().size() == 0)//If we have no Light Source we cant see objects anyways
			return;

		//Select one random light source from all lights
		int lightIndex = m_dist(m_gen) * pathtracer.getScene().getLights().size();
		auto light = pathtracer.getScene().getLights()[lightIndex].get();
		float attenuation;

		//////////////////////////////////////////////
		//  R-PATH: Scattered Light calculation (RGB)
		//////////////////////////////////////////////

		//calculate ray from intersection point towards light source
		KIRK::Ray hitToLight = light->calcLightdir(hit.m_location, attenuation, true);
		//use the output_ray to pass the hit_to_light direction to the bsdf sampling method
		result_direction_r = hitToLight.m_direction;

		// Sample bsdf for r-path
		glm::vec3 scattering_r = bsdf->sample(hit, input_ray_normalized, hit.m_normal, result_direction_r, pdf_r, resultBounce.mat_flags, sample);
		//glm::vec3 scattering_r = SpecularReflectionBSDF::localSample(hit, -input_ray_normalized, hit.m_normal, sample, result_direction_r, pdf_r, resultBounce.mat_flags);

		//////////////////////////////////////////////
		//  TT-PATH: Scattered Light calculation (RGB)
		//////////////////////////////////////////////

		//refract input ray on surface and put it into new intersection
		glm::vec3 t_dir = glm::refract(input_ray_normalized, glm::faceforward(hit.m_normal, input_ray_normalized, hit.m_normal), 1.0f / material->m_ior);
		glm::vec3 t_normal;
		Intersection t_hit(Ray(hit.m_location + 1e-4f * t_dir, t_dir));

		//get the intersection point on the second wall
		if (pathtracer.getScene().getDataStructure().closestIntersection(&t_hit))
		{
			hit.m_object->calcNormal(&t_hit);// calculation of the normal
			t_normal = t_hit.m_normal;
		}
		else //if we dont find an intersection we use the old normal. Other direction cause we are at the second wall
		{
			t_normal = -hit.m_normal;
			t_hit.m_object = hit.m_object;
		}

		//calculate ray from intersection point towards light source
		hitToLight = light->calcLightdir(t_hit.m_location, attenuation, true);
		//use the output_ray to pass the hit_to_light direction to the bsdf sampling method
		result_direction_tt = hitToLight.m_direction;

		//set the material flags to TT-Path
		resultBounce.mat_flags = BSDFHelper::MATFLAG_CYLINDER_T_BOUNCE;
		// Sample bsdf for tt-path
		glm::vec3 scattering_tt = bsdf->sample(t_hit, input_ray_normalized, t_normal, result_direction_tt, pdf_tt, resultBounce.mat_flags, sample);

		//////////////////////////////////////////////
		//  TRT-PATH: Scattered Light calculation (RGB)
		//////////////////////////////////////////////

		//Reflect on the second cylinder wall 
		glm::vec3 tr_dir = glm::reflect(glm::normalize(t_dir), faceforward(t_normal, glm::normalize(t_dir), t_normal));
		glm::vec3 tr_normal;
		//New Intersection for TR path
		Intersection tr_hit(Ray(t_hit.m_location + 1e-4f * tr_dir, tr_dir));
		//get the intersection for the second point on the first wall
		if (pathtracer.getScene().getDataStructure().closestIntersection(&tr_hit))
		{
			hit.m_object->calcNormal(&tr_hit);// calculation of the normal
			tr_normal = tr_hit.m_normal;
		}
		else//if we dont find an intersection we use the old normal
		{
			tr_normal = hit.m_normal;
			tr_hit.m_object = hit.m_object;
		}

		//calculate ray from intersection point towards light source
		hitToLight = light->calcLightdir(tr_hit.m_location, attenuation, true);
		//use the output_ray to pass the hit_to_light direction to the bsdf sampling method
		result_direction_trt = hitToLight.m_direction;

		//set the material flags to TRT-Path
		resultBounce.mat_flags = BSDFHelper::MATFLAG_CYLINDER_TR_BOUNCE;
		// Sample bsdf for trt-path
		glm::vec3 scattering_trt = bsdf->sample(tr_hit, input_ray_normalized, tr_normal, result_direction_trt, pdf_trt, resultBounce.mat_flags, sample);

		//////////////////////////////////////////////
		//  Final Color and outgoing path calculation
		//////////////////////////////////////////////

		//pdf mean
		pdf = (pdf_r + pdf_tt + pdf_trt) / 3.f;

		//decide which direction we follow randomly
		std::uniform_int_distribution<> dis(0, 1);
		int path = dis(m_gen);
		//path = 0;
		if (path == 0)//we take R-Path as outgoing direction
		{
			resultRay = KIRK::Ray(hit.m_location + 1e-4f * result_direction_r, result_direction_r);
		}
		else if (path == 1)//we take TRT-Path as outgoing direction
		{
			resultRay = KIRK::Ray(tr_hit.m_location + 1e-4f * result_direction_trt, result_direction_trt);			
		}
		else//we take TT-Path as outgoing direction
		{
			resultRay = KIRK::Ray(t_hit.m_location + 1e-4f * result_direction_tt, result_direction_tt);			
		}

		// Scattering Integral of outgoing radiance (Marschner Paper Formel 1)
		if ((scattering_r + scattering_tt + scattering_trt) == glm::vec3(0.f) || pdf <= 1E-4f
			|| std::max(resultBounce.radiance.x, std::max(resultBounce.radiance.y, resultBounce.radiance.z)) < 0.01f)
			resultBounce.radiance = glm::vec3(0);
		else//Render-Equation (Marschner Equation 1). Theta_i is stored in sample.x and pdf is to account for radiance fallof
		{
			//glm::vec3 r_radiance = glm::min(scattering_r * glm::abs(glm::dot(result_direction_r, hit.m_normal)) / pdf_r, 1.f);
			glm::vec3 mar_radiance =  ((scattering_r + scattering_trt + scattering_tt) * fiber_width * glm::abs(glm::cos(sample.x)) / pdf);
			if (mar_radiance.x >= 1.f || mar_radiance.y >= 1.f || mar_radiance.z >= 1.f)
				mar_radiance = glm::normalize(mar_radiance);
			resultBounce.radiance *= glm::min(mar_radiance , 1.0f);
		}


		//Color vectors
		Color::RGBA accumulatedColor(0.0f);
		Color::RGBA ambientLight;
		Color::RGBA directLight;

		//Ambient light
		ambientLight = pathtracer.getScene().getEnvironment().getAmbientLight() * glm::vec4(bsdf->evaluateLight(hit, hit.m_normal, hit.m_normal) * glm::one_over_pi<float>(), 1.0f);
		accumulatedColor += ambientLight * glm::vec4(resultBounce.radiance, 1.0f);

		//Calculate direct Light from all light Sources. 
		directLight = calcDirectLight(pathtracer, hit, pathtracer.getSampler().sample2D());
		accumulatedColor += directLight * glm::vec4(resultBounce.radiance, 1.0f);//add direct light to color

		//Add new color to bounce
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
		KIRK::Ray hitToLight = light->calcLightdir(hit.m_location, attenuation, true);//calculates ray from intersection point towards light source
		glm::vec3 lightpos = hitToLight.m_origin + hitToLight.m_direction;
		hitToLight.m_origin += 1e-4f * glm::faceforward(hit.m_normal, hitToLight.m_origin - lightpos, hit.m_normal);
		hitToLight.m_direction = glm::normalize(hitToLight.m_direction);

		auto lightColor = light->m_color;

		if (light->m_color.x > 0 || light->m_color.y > 0 || light->m_color.z > 0)
		{
			lightColor *= attenuation//attenuation factor of the light. Calculated in light->calcLightdir()
				* glm::vec4(hit.m_object->getMaterial()->m_bsdf->evaluateLight(hit, hitToLight.m_direction, hit.m_ray.m_direction), 1.0f)//light influence at hit point. calc by the bsdf
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
