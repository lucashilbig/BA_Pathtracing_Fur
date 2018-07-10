#include "Bsdf.h"
#include "KIRK/Common/Intersection.h"
#define _USE_MATH_DEFINES
#include <math.h>

namespace KIRK {
////////////////////////////////////////////////////////////////////////////////////
//////////////
//////////////	BSDF HELPER FUNCTIONS
//////////////
////////////////////////////////////////////////////////////////////////////////////

float BSDFHelper::schlickFresnel (const glm::vec3& view, const glm::vec3& normal, float ior_in, float ior_out)
{
	// for fresnel we have to calculate the amount of light that is refracted if the angle is 0 (orthogonal to surface)
	float R0 = std::pow((ior_in - ior_out) / (ior_in + ior_out), 2.0); // before default value R0 = 0.04
	// for angles between 0 and 90
	return R0 + std::pow(1.0 - glm::dot(view, normal), 5.0) * (1.0 - R0);
}

glm::vec2 BSDFHelper::concentricSampleDisk (const glm::vec2& randoms)
{
	glm::vec2 offset = 2.f * randoms - glm::vec2(1, 1); // offset random numbers to range [ -1; 1 ]^2

	if(offset.x == 0 && offset.y == 0) { return glm::vec2(0, 0); }

	float theta, random;

	if(glm::abs(offset.x) > glm::abs(offset.y))
	{
		random = offset.x;
		theta = glm::quarter_pi <float>() * (offset.y / offset.x);
	}
	else
	{
		random = offset.y;
		theta = glm::half_pi <float>() - glm::quarter_pi <float>() * (offset.x / offset.y);
	}

	return random * glm::vec2(glm::cos(theta), glm::sin(theta));
}

glm::vec3 BSDFHelper::sampleAngle(const glm::vec2& randoms, float max_angle)
{
	float phi = randoms.x * 2.f * M_PI;
	float cosTheta = 1.0f - randoms.y * (1.0f - std::cos(max_angle));
	float sinTheta = std::sqrt(1.f - cosTheta * cosTheta);
	return{ std::cos(phi) * sinTheta, std::sin(phi) * sinTheta, cosTheta };
}

glm::vec3 BSDFHelper::cosineSampleHemisphere (const glm::vec2& s)
{
	// Malley's method: Uniformly sample points on the unit disk and project them up to the unit sphere.
	glm::vec2 d = concentricSampleDisk(s);
	float z = glm::sqrt(glm::max(0.f, 1 - d.x * d.x - d.y * d.y));

	return glm::vec3(d.x, d.y, z);
}

glm::vec3 BSDFHelper::uniformSphereSample (float u, float v)
{
	float phi = v * 2.f * glm::pi<float>();
	float cosTheta = 2 * u - 1;
	float sinTheta = glm::sqrt(glm::max(0.0f, 1.0f - cosTheta * cosTheta));
	return glm::vec3(sinTheta * glm::cos(phi), sinTheta * glm::sin(phi), cosTheta);

}

float BSDFHelper::dialectricFresnel (const float cos_theta, float eta_i, float eta_t)
{
	float cos_theta_i = glm::clamp(cos_theta, -1.f, 1.f);

	// If we are not entering, we have to swap the indices of refraction:
	if(cos_theta_i <= 0)
	{
		float tmp = eta_i;
		eta_i = eta_t;
		eta_t = tmp;
		cos_theta_i = glm::abs(cos_theta_i);
	}

	// Snell's law
	float sin_theta_i = glm::sqrt(glm::max(0.f, 1 - cos_theta_i * cos_theta_i));
	float sin_theta_t = eta_i / eta_t * sin_theta_i;

	if(sin_theta_t >= 1.0f) { return 1.0f; }

	float cos_theta_t = glm::sqrt(glm::max(0.f, 1.0f - sin_theta_t * sin_theta_t));

	float rparl = ((eta_t * cos_theta_i) - (eta_i * cos_theta_t)) /
	              ((eta_t * cos_theta_i) + (eta_i * cos_theta_t));
	float rperp = ((eta_i * cos_theta_i) - (eta_t * cos_theta_t)) /
	              ((eta_i * cos_theta_i) + (eta_t * cos_theta_t));

	// Valid only for unpolarised light, which is what we assume here.
	return (rparl * rparl + rperp * rperp) / 2.0f;
}

////////////////////////////////////////////////////////////////////////////////////
//////////////
//////////////	DIFFUSE BSDF
//////////////
////////////////////////////////////////////////////////////////////////////////////

glm::vec3 BSDF::sample (const Intersection& hit, const glm::vec3& ray_in, const glm::vec3& normal, glm::vec3& ray_out, float& pdf, int& mat_flags, const glm::vec2& sampling, bool useRadianceOverImportance) const
{
	if(glm::dot(ray_in, normal) == 0) { return glm::vec3(0.f); }

	return m_localSample(hit, ray_in, normal, sampling, ray_out, pdf, mat_flags, useRadianceOverImportance);
}

glm::vec3 LambertianReflectionBSDF::localSample (const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& normal, glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance)
{
	bool entering = glm::dot(local_input_ray, normal) > 0.f;
	local_output_ray = Math::localToWorldNormal(entering ? BSDFHelper::cosineSampleHemisphere(sample) : - BSDFHelper::cosineSampleHemisphere(sample), normal);
	output_pdf = glm::abs(glm::dot(local_output_ray, normal)) * glm::one_over_pi <float>();
	mat_flags = 0;
	if(output_pdf == 0.f) { return glm::vec3(0.f); }

	return glm::vec3(hit.m_object->getMaterial()->fetchParameterColor <MatParamType::DIFFUSE>(hit.m_texcoord)) * glm::one_over_pi <float>();
}

glm::vec3 LambertianReflectionBSDF::evaluateLight (const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray)
{
	bool reflect = glm::dot(local_input_ray, hit.m_normal) * glm::dot(local_output_ray, hit.m_normal) > 0;
	if(reflect) { return glm::vec3(hit.m_object->getMaterial()->fetchParameterColor <MatParamType::DIFFUSE>(hit.m_texcoord)) * glm::one_over_pi <float>(); }
	return glm::vec3(0);
}

////////////////////////////////////////////////////////////////////////////////////
//////////////
//////////////	REFLECTIVE BSDF
//////////////
////////////////////////////////////////////////////////////////////////////////////

glm::vec3 SpecularReflectionBSDF::localSample (const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& normal, glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance)
{
	local_output_ray = glm::reflect(-local_input_ray, glm::faceforward(normal, -local_input_ray, normal));
	output_pdf = 1.f;
	mat_flags |= BSDFHelper::MATFLAG_SPECULAR_BOUNCE;

	return glm::vec3(hit.m_object->getMaterial()->fetchParameterColor <MatParamType::SPECULAR>(hit.m_texcoord)) / glm::abs(glm::dot(local_output_ray, normal));
}

glm::vec3 SpecularReflectionBSDF::evaluateLight (const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray) { return glm::vec3(.0f); }

////////////////////////////////////////////////////////////////////////////////////
//////////////
//////////////	GLOSSY BSDF
//////////////
////////////////////////////////////////////////////////////////////////////////////

glm::vec3 GlossyBSDF::localSample(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& normal, glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance)
{
	float rad = glm::radians(180.0f - (1.0f - hit.m_object->getMaterial()->fetchParameterFloat <MatParamType::ROUGHNESS>(hit.m_texcoord)) * 180.0f);
	glm::vec3 reflected = reflect(-local_input_ray, faceforward(normal, -local_input_ray, normal));
	glm::vec3 sampledPoint = BSDFHelper::sampleAngle(sample, rad);

	local_output_ray = Math::localToWorldNormal(sampledPoint, reflected);

	if (dot(local_output_ray, faceforward(normal, -local_input_ray, normal)) < 0.0f)
	{
		local_output_ray = Math::localToWorldNormal(sampledPoint * glm::vec3(-1, -1, 1), reflected);
	}

	output_pdf = 1;

	mat_flags |= BSDFHelper::MATFLAG_SPECULAR_BOUNCE;

	return glm::vec3(hit.m_object->getMaterial()->fetchParameterColor <MatParamType::SPECULAR>(hit.m_texcoord)) / glm::abs(glm::dot(local_output_ray, normal));
}

glm::vec3 GlossyBSDF::evaluateLight(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray)
{
	return glm::vec3(0);
}

////////////////////////////////////////////////////////////////////////////////////
//////////////
//////////////	REFRACTIVE BSDF
//////////////
////////////////////////////////////////////////////////////////////////////////////

glm::vec3 SpecularTransmissionBSDF::localSample(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& normal, glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance)
{

	Material* material = hit.m_object->getMaterial();

	//Needed parameters for refraction and reflection
	bool entering = glm::dot(local_input_ray, normal) > 0.f;
	float eta_i = entering ? 1.f : material->m_ior;
	float eta_t = entering ? material->m_ior : 1.f;
	float fresnel = BSDFHelper::dialectricFresnel(glm::abs(glm::dot(local_input_ray, normal)), eta_i, eta_t);

	mat_flags |= BSDFHelper::MATFLAG_SPECULAR_BOUNCE;

	//refract ray.
	local_output_ray = glm::refract(glm::normalize(-local_input_ray), glm::faceforward(normal, -local_input_ray, normal), eta_i / eta_t);

	output_pdf = 1.0f;

	// Choose whether to sample refractive or reflective BSDF
	if (local_output_ray != glm::vec3(0) && !std::isnan(local_output_ray.x))
	{
		mat_flags |= BSDFHelper::MATFLAG_TRANSPARENT_BOUNCE;

		glm::vec3 ft = glm::vec3(hit.m_object->getMaterial()->fetchParameterColor <MatParamType::VOLUME>(hit.m_texcoord)) * (glm::vec3(1, 1, 1) - fresnel);

		if (useRadianceOverImportance) { ft *= (eta_i * eta_i) / (eta_t * eta_t); }

		return ft / glm::abs(glm::dot(local_output_ray, normal));
	}
	return glm::vec3(0);
}

glm::vec3 SpecularTransmissionBSDF::evaluateLight(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray) { return glm::vec3(0); }

////////////////////////////////////////////////////////////////////////////////////
//////////////
//////////////	TRANSLUCENT BSDF
//////////////
////////////////////////////////////////////////////////////////////////////////////

glm::vec3 LambertianTransmissionBSDF::localSample (const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& normal, glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance)
{
	bool entering = glm::dot(local_input_ray, normal) > 0.f;
	local_output_ray = Math::localToWorldNormal(entering ? -BSDFHelper::cosineSampleHemisphere(sample) : BSDFHelper::cosineSampleHemisphere(sample), normal);
	output_pdf = glm::abs(glm::dot(local_output_ray, normal)) * glm::one_over_pi <float>();
	mat_flags = BSDFHelper::MATFLAG_TRANSPARENT_BOUNCE;

	if(output_pdf == 0.f) { return glm::vec3(0.f); }

	return glm::vec3(hit.m_object->getMaterial()->fetchParameterColor <MatParamType::VOLUME>(hit.m_texcoord)) * glm::one_over_pi <float>();
}

glm::vec3 LambertianTransmissionBSDF::evaluateLight (const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray)
{
	bool reflect = glm::dot(local_input_ray, hit.m_normal) * glm::dot(local_output_ray, hit.m_normal) > 0;

	if(!reflect)
		return glm::vec3(hit.m_object->getMaterial()->fetchParameterColor <MatParamType::DIFFUSE>(hit.m_texcoord)) * glm::one_over_pi <float>();

	return glm::vec3(0);
}

////////////////////////////////////////////////////////////////////////////////////
//////////////
//////////////	GLASS BSDF
//////////////
////////////////////////////////////////////////////////////////////////////////////

glm::vec3 GlassBSDF::localSample (const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& normal, glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance)
{

	Material* material = hit.m_object->getMaterial();

	//Needed parameters for refraction and reflection
	bool entering = glm::dot(local_input_ray, normal) > 0.f;
	float eta_i = entering ? 1.f : material->m_ior;
	float eta_t = entering ? material->m_ior : 1.f;
	float fresnel = BSDFHelper::dialectricFresnel(glm::abs(glm::dot(glm::normalize(local_input_ray), normal)), eta_i, eta_t);

	mat_flags |= BSDFHelper::MATFLAG_SPECULAR_BOUNCE;

	//refract ray.
	local_output_ray = glm::refract(glm::normalize(-local_input_ray), glm::faceforward(normal, -glm::normalize(local_input_ray), normal), eta_i / eta_t);

	// Choose whether to sample refractive or reflective BSDF
	if(local_output_ray != glm::vec3(0) && sample.y > fresnel && !std::isnan(local_output_ray.x))
	{
		mat_flags |= BSDFHelper::MATFLAG_TRANSPARENT_BOUNCE;

		output_pdf = 1 - fresnel;
		glm::vec3 ft = glm::vec3(hit.m_object->getMaterial()->fetchParameterColor <MatParamType::VOLUME>(hit.m_texcoord)) * (glm::vec3(1, 1, 1) - fresnel);

		if(useRadianceOverImportance) { ft *= (eta_i * eta_i) / (eta_t * eta_t); }

		return ft / glm::abs(glm::dot(local_output_ray, normal));
	}
	local_output_ray = glm::reflect(glm::normalize(-local_input_ray), glm::faceforward(normal, -glm::normalize(local_input_ray), normal));
	output_pdf = fresnel;
	return fresnel * glm::vec3(hit.m_object->getMaterial()->fetchParameterColor <MatParamType::SPECULAR>(hit.m_texcoord)) / glm::abs(glm::dot(local_output_ray, normal));
}

glm::vec3 GlassBSDF::evaluateLight (const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray) { return glm::vec3(0); }

////////////////////////////////////////////////////////////////////////////////////
//////////////
//////////////	MILKGLASS BSDF
//////////////
////////////////////////////////////////////////////////////////////////////////////

glm::vec3 MilkGlassBSDF::localSample(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& normal, glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance)
{

	Material* material = hit.m_object->getMaterial();

	//Needed parameters for refraction and reflection
	bool entering = glm::dot(local_input_ray, normal) > 0.f;
	float eta_i = entering ? 1.f : material->m_ior;
	float eta_t = entering ? material->m_ior : 1.f;
	float fresnel = BSDFHelper::dialectricFresnel(glm::abs(glm::dot(glm::normalize(local_input_ray), normal)), eta_i, eta_t);

	mat_flags |= BSDFHelper::MATFLAG_SPECULAR_BOUNCE;

	//refract ray.
	glm::vec3 refracted = glm::refract(glm::normalize(-local_input_ray), glm::faceforward(normal, -glm::normalize(local_input_ray), normal), eta_i / eta_t);

	// Choose whether to sample refractive or reflective BSDF
	if (refracted != glm::vec3(0) && sample.y > fresnel && !std::isnan(refracted.x))
	{
		float rad = glm::radians(180.0f - (1.0f - hit.m_object->getMaterial()->fetchParameterFloat <MatParamType::ROUGHNESS>(hit.m_texcoord)) * 180.0f);
		glm::vec3 sampledPoint = BSDFHelper::sampleAngle(sample, rad);

		local_output_ray = Math::localToWorldNormal(sampledPoint, refracted);

		if (dot(local_output_ray, faceforward(normal, -local_input_ray, normal)) > 0.0f)
		{
			local_output_ray = Math::localToWorldNormal(sampledPoint * glm::vec3(-1, -1, 1), refracted);
		}
		mat_flags |= BSDFHelper::MATFLAG_TRANSPARENT_BOUNCE;

		output_pdf = 1 - fresnel;
		glm::vec3 ft = glm::vec3(hit.m_object->getMaterial()->fetchParameterColor <MatParamType::VOLUME>(hit.m_texcoord)) * (glm::vec3(1, 1, 1) - fresnel);

		if (useRadianceOverImportance) { ft *= (eta_i * eta_i) / (eta_t * eta_t); }

		return ft / glm::abs(glm::dot(local_output_ray, normal));
	}
	float rad = glm::radians(180.0f - (1.0f - hit.m_object->getMaterial()->fetchParameterFloat <MatParamType::ROUGHNESS>(hit.m_texcoord)) * 180.0f);
	glm::vec3 reflected = reflect(-local_input_ray, faceforward(normal, -local_input_ray, normal));
	glm::vec3 sampledPoint = BSDFHelper::sampleAngle(sample, rad);

	local_output_ray = Math::localToWorldNormal(sampledPoint, reflected);

	if (dot(local_output_ray, faceforward(normal, -local_input_ray, normal)) < 0.0f)
	{
		local_output_ray = Math::localToWorldNormal(sampledPoint * glm::vec3(-1, -1, 1), reflected);
	}
	output_pdf = fresnel;
	return fresnel * glm::vec3(hit.m_object->getMaterial()->fetchParameterColor <MatParamType::SPECULAR>(hit.m_texcoord)) / glm::abs(glm::dot(local_output_ray, normal));
}

glm::vec3 MilkGlassBSDF::evaluateLight(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray) { return glm::vec3(0); }


////////////////////////////////////////////////////////////////////////////////////
//////////////
//////////////	EMISSION BSDF
//////////////
////////////////////////////////////////////////////////////////////////////////////

glm::vec3 EmissionBSDF::localSample (const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& normal, glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance)
{
	// PDF > 0 and a return value length > 0 means emission.
	output_pdf = 1.f;
	local_output_ray = glm::vec3(0);
	mat_flags = BSDFHelper::MATFLAG_EMISSIVE_BOUNCE;

	return glm::vec3(1);
}

glm::vec3 EmissionBSDF::evaluateLight (const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray) { return glm::vec3(.0f); }

////////////////////////////////////////////////////////////////////////////////////
//////////////
//////////////	TRANSPARENT BSDF
//////////////
////////////////////////////////////////////////////////////////////////////////////

glm::vec3 TransparentBSDF::localSample (const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& normal, glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance)
{
	Material* mat = hit.m_object->getMaterial();
	glm::vec2 texcoord = hit.m_texcoord;
	local_output_ray = -local_input_ray;
	mat_flags = BSDFHelper::MATFLAG_TRANSPARENT_BOUNCE;
	mat_flags |= BSDFHelper::MATFLAG_SPECULAR_BOUNCE;
	output_pdf = 1;
	return glm::vec3(mat->fetchParameterColor <MatParamType::VOLUME>(texcoord)) / glm::abs(glm::dot(local_output_ray, normal));
}

glm::vec3 TransparentBSDF::evaluateLight (const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray) { return glm::vec3(.0f); }


////////////////////////////////////////////////////////////////////////////////////
//////////////
//////////////	MARSCHNER HAIR BSDF
//////////////
////////////////////////////////////////////////////////////////////////////////////

glm::vec3 MarschnerHairBSDF::localSample(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& normal, glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance)
{
	//todo: implement
	return glm::vec3(.0f);
}

glm::vec3 MarschnerHairBSDF::evaluateLight(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray) { return glm::vec3(.0f);/*todo: implement*/ }
}
