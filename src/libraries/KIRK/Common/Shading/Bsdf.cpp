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

	int BSDFHelper::SolveP3(double *x, double a, double b, double c) {	// solve cubic equation x^3 + a*x^2 + b*x + c = 0
		const double TwoPi = 6.28318530717958648;
		const double eps = 1e-14;
		double a2 = a * a;
		double q = (a2 - 3 * b) / 9;
		double r = (a*(2 * a2 - 9 * b) + 27 * c) / 54;
		// equation x^3 + q*x + r = 0
		double r2 = r * r;
		double q3 = q * q*q;
		double A, B;
		if (r2 <= (q3 + eps)) {//<<-- FIXED!
			double t = r / sqrt(q3);
			if (t < -1) t = -1;
			if (t > 1) t = 1;
			t = acos(t);
			a /= 3; q = -2 * sqrt(q);
			x[0] = q * cos(t / 3) - a;
			x[1] = q * cos((t + TwoPi) / 3) - a;
			x[2] = q * cos((t - TwoPi) / 3) - a;
			return(3);
		}
		else {
			//A =-pow(fabs(r)+sqrt(r2-q3),1./3); 
			A = -root3(fabs(r) + sqrt(r2 - q3));
			if (r < 0) A = -A;
			B = (A == 0 ? 0 : B = q / A);

			a /= 3;
			x[0] = (A + B) - a;
			x[1] = -0.5*(A + B) - a;
			x[2] = 0.5*sqrt(3.)*(A - B);
			if (fabs(x[2]) < eps) { x[2] = x[1]; return(2); }
			return(1);
		}
	}

	double BSDFHelper::_root3(double x)
	{
		double s = 1.;
		while (x < 1.)
		{
			x *= 8.;
			s *= 0.5;
		}
		while (x > 8.)
		{
			x *= 0.125;
			s *= 2.;
		}
		double r = 1.5;
		r -= 1. / 3. * (r - x / (r * r));
		r -= 1. / 3. * (r - x / (r * r));
		r -= 1. / 3. * (r - x / (r * r));
		r -= 1. / 3. * (r - x / (r * r));
		r -= 1. / 3. * (r - x / (r * r));
		r -= 1. / 3. * (r - x / (r * r));
		return r * s;
	}

	double BSDFHelper::root3(double x)
	{
		if (x > 0) return _root3(x); else
			if (x < 0) return-_root3(-x); else
				return 0.;
	}

	float BSDFHelper::normal_gauss_pdf(float x, float mean, float stddev)
	{
		static const float inv_sqrt_2pi = 0.3989422804014327;
		float a = (x - mean) / stddev;

		return inv_sqrt_2pi / stddev * std::exp(-0.5f * (a * a));
	}

	float BSDFHelper::schlickFresnel(const glm::vec3& view, const glm::vec3& normal, float ior_in, float ior_out)
	{
		// for fresnel we have to calculate the amount of light that is refracted if the angle is 0 (orthogonal to surface)
		float R0 = std::pow((ior_in - ior_out) / (ior_in + ior_out), 2.0); // before default value R0 = 0.04
		// for angles between 0 and 90
		return R0 + std::pow(1.0 - glm::dot(view, normal), 5.0) * (1.0 - R0);
	}

	glm::vec2 BSDFHelper::concentricSampleDisk(const glm::vec2& randoms)
	{
		glm::vec2 offset = 2.f * randoms - glm::vec2(1, 1); // offset random numbers to range [ -1; 1 ]^2

		if (offset.x == 0 && offset.y == 0) { return glm::vec2(0, 0); }

		float theta, random;

		if (glm::abs(offset.x) > glm::abs(offset.y))
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

	glm::vec3 BSDFHelper::cosineSampleHemisphere(const glm::vec2& s)
	{
		// Malley's method: Uniformly sample points on the unit disk and project them up to the unit sphere.
		glm::vec2 d = concentricSampleDisk(s);
		float z = glm::sqrt(glm::max(0.f, 1 - d.x * d.x - d.y * d.y));

		return glm::vec3(d.x, d.y, z);
	}

	glm::vec3 BSDFHelper::uniformSphereSample(float u, float v)
	{
		float phi = v * 2.f * glm::pi<float>();
		float cosTheta = 2 * u - 1;
		float sinTheta = glm::sqrt(glm::max(0.0f, 1.0f - cosTheta * cosTheta));
		return glm::vec3(sinTheta * glm::cos(phi), sinTheta * glm::sin(phi), cosTheta);

	}

	float BSDFHelper::dialectricFresnel(const float cos_theta, float eta_i, float eta_t)
	{
		float cos_theta_i = glm::clamp(cos_theta, -1.f, 1.f);

		// If we are not entering, we have to swap the indices of refraction:
		if (cos_theta_i <= 0)
		{
			float tmp = eta_i;
			eta_i = eta_t;
			eta_t = tmp;
			cos_theta_i = glm::abs(cos_theta_i);
		}

		// Snell's law
		float sin_theta_i = glm::sqrt(glm::max(0.f, 1 - cos_theta_i * cos_theta_i));
		float sin_theta_t = eta_i / eta_t * sin_theta_i;

		if (sin_theta_t >= 1.0f) { return 1.0f; }

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

	glm::vec3 BSDF::sample(const Intersection& hit, const glm::vec3& ray_in, const glm::vec3& normal, glm::vec3& ray_out, float& pdf, int& mat_flags, glm::vec2& sampling, bool useRadianceOverImportance) const
	{
		if (glm::dot(ray_in, normal) == 0) { return glm::vec3(0.f); }

		return m_localSample(hit, ray_in, normal, sampling, ray_out, pdf, mat_flags, useRadianceOverImportance);
	}

	glm::vec3 LambertianReflectionBSDF::localSample(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& normal, glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance)
	{
		bool entering = glm::dot(local_input_ray, normal) > 0.f;
		local_output_ray = Math::localToWorldNormal(entering ? BSDFHelper::cosineSampleHemisphere(sample) : -BSDFHelper::cosineSampleHemisphere(sample), normal);
		output_pdf = glm::abs(glm::dot(local_output_ray, normal)) * glm::one_over_pi <float>();
		mat_flags = 0;
		if (output_pdf == 0.f) { return glm::vec3(0.f); }

		return glm::vec3(hit.m_object->getMaterial()->fetchParameterColor <MatParamType::DIFFUSE>(hit.m_texcoord)) * glm::one_over_pi <float>();
	}

	glm::vec3 LambertianReflectionBSDF::evaluateLight(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray)
	{
		bool reflect = glm::dot(local_input_ray, hit.m_normal) * glm::dot(local_output_ray, hit.m_normal) > 0;
		if (reflect) { return glm::vec3(hit.m_object->getMaterial()->fetchParameterColor <MatParamType::DIFFUSE>(hit.m_texcoord)) * glm::one_over_pi <float>(); }
		return glm::vec3(0);
	}

	////////////////////////////////////////////////////////////////////////////////////
	//////////////
	//////////////	REFLECTIVE BSDF
	//////////////
	////////////////////////////////////////////////////////////////////////////////////

	glm::vec3 SpecularReflectionBSDF::localSample(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& normal, glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance)
	{
		local_output_ray = glm::reflect(-local_input_ray, glm::faceforward(normal, -local_input_ray, normal));
		output_pdf = 1.f;
		mat_flags |= BSDFHelper::MATFLAG_SPECULAR_BOUNCE;

		return glm::vec3(hit.m_object->getMaterial()->fetchParameterColor <MatParamType::SPECULAR>(hit.m_texcoord)) / glm::abs(glm::dot(local_output_ray, normal));
	}

	glm::vec3 SpecularReflectionBSDF::evaluateLight(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray) { return glm::vec3(.0f); }

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

	glm::vec3 LambertianTransmissionBSDF::localSample(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& normal, glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance)
	{
		bool entering = glm::dot(local_input_ray, normal) > 0.f;
		local_output_ray = Math::localToWorldNormal(entering ? -BSDFHelper::cosineSampleHemisphere(sample) : BSDFHelper::cosineSampleHemisphere(sample), normal);
		output_pdf = glm::abs(glm::dot(local_output_ray, normal)) * glm::one_over_pi <float>();
		mat_flags = BSDFHelper::MATFLAG_TRANSPARENT_BOUNCE;

		if (output_pdf == 0.f) { return glm::vec3(0.f); }

		return glm::vec3(hit.m_object->getMaterial()->fetchParameterColor <MatParamType::VOLUME>(hit.m_texcoord)) * glm::one_over_pi <float>();
	}

	glm::vec3 LambertianTransmissionBSDF::evaluateLight(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray)
	{
		bool reflect = glm::dot(local_input_ray, hit.m_normal) * glm::dot(local_output_ray, hit.m_normal) > 0;

		if (!reflect)
			return glm::vec3(hit.m_object->getMaterial()->fetchParameterColor <MatParamType::DIFFUSE>(hit.m_texcoord)) * glm::one_over_pi <float>();

		return glm::vec3(0);
	}

	////////////////////////////////////////////////////////////////////////////////////
	//////////////
	//////////////	GLASS BSDF
	//////////////
	////////////////////////////////////////////////////////////////////////////////////

	glm::vec3 GlassBSDF::localSample(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& normal, glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance)
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
		if (local_output_ray != glm::vec3(0) && sample.y > fresnel && !std::isnan(local_output_ray.x))
		{
			mat_flags |= BSDFHelper::MATFLAG_TRANSPARENT_BOUNCE;

			output_pdf = 1 - fresnel;
			glm::vec3 ft = glm::vec3(hit.m_object->getMaterial()->fetchParameterColor <MatParamType::VOLUME>(hit.m_texcoord)) * (glm::vec3(1, 1, 1) - fresnel);

			if (useRadianceOverImportance) { ft *= (eta_i * eta_i) / (eta_t * eta_t); }

			return ft / glm::abs(glm::dot(local_output_ray, normal));
		}
		local_output_ray = glm::reflect(glm::normalize(-local_input_ray), glm::faceforward(normal, -glm::normalize(local_input_ray), normal));
		output_pdf = fresnel;
		return fresnel * glm::vec3(hit.m_object->getMaterial()->fetchParameterColor <MatParamType::SPECULAR>(hit.m_texcoord)) / glm::abs(glm::dot(local_output_ray, normal));
	}

	glm::vec3 GlassBSDF::evaluateLight(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray) { return glm::vec3(0); }

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

	glm::vec3 EmissionBSDF::localSample(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& normal, glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance)
	{
		// PDF > 0 and a return value length > 0 means emission.
		output_pdf = 1.f;
		local_output_ray = glm::vec3(0);
		mat_flags = BSDFHelper::MATFLAG_EMISSIVE_BOUNCE;

		return glm::vec3(1);
	}

	glm::vec3 EmissionBSDF::evaluateLight(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray) { return glm::vec3(.0f); }

	////////////////////////////////////////////////////////////////////////////////////
	//////////////
	//////////////	TRANSPARENT BSDF
	//////////////
	////////////////////////////////////////////////////////////////////////////////////

	glm::vec3 TransparentBSDF::localSample(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& normal, glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance)
	{
		Material* mat = hit.m_object->getMaterial();
		glm::vec2 texcoord = hit.m_texcoord;
		local_output_ray = -local_input_ray;
		mat_flags = BSDFHelper::MATFLAG_TRANSPARENT_BOUNCE;
		mat_flags |= BSDFHelper::MATFLAG_SPECULAR_BOUNCE;
		output_pdf = 1;
		return glm::vec3(mat->fetchParameterColor <MatParamType::VOLUME>(texcoord)) / glm::abs(glm::dot(local_output_ray, normal));
	}

	glm::vec3 TransparentBSDF::evaluateLight(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray) { return glm::vec3(.0f); }


	////////////////////////////////////////////////////////////////////////////////////
	//////////////
	//////////////	MARSCHNER HAIR BSDF
	//////////////
	////////////////////////////////////////////////////////////////////////////////////

	glm::vec3 MarschnerHairBSDF::localSample(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& normal, glm::vec2& sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance)
	{
		//get our object
		KIRK::Object *cylinder_obj = hit.m_object;

		//needed parameters
		Material* material = hit.m_object->getMaterial();
		glm::vec2 texcoord = hit.m_texcoord;

		//calculation of tangent vector
		glm::vec3 c1 = glm::cross(normal, glm::vec3(0.0, 0.0, 1.0));
		glm::vec3 c2 = glm::cross(normal, glm::vec3(0.0, 1.0, 0.0));
		glm::vec3 tangent = (glm::length(c1) > glm::length(c2)) ? glm::normalize(c1) : glm::normalize(c2);

		//Get lobe alpha and beta which is stored in the sample parameter
		float alpha_r = sample.x;//longitudinal shift. Suggested value from marschner hair paper between -10 and -5 degrees
		float beta_r = sample.y;//longitudinal width (stdev.). Suggested value from marschner hair paper between 5 and 10 degrees		

		//We are inside the cylinder and we want to take Cylinder Path: TT (refractive transmission, p = 1)
		if (mat_flags & BSDFHelper::MATFLAG_CYLINDER_T_BOUNCE)
		{
			//refract on the second wall of the cylinder to get output ray
			local_output_ray = glm::refract(local_input_ray, glm::faceforward(normal, local_input_ray, normal), 1.f / material->m_ior);

			//rotate towards normal to account for the tilted fiber surface
			local_output_ray = glm::vec3(glm::vec4(local_output_ray, 0.f) * glm::rotate(glm::radians(-alpha_r / 2), tangent));

			//reset material flags and set it to finished path since we leave the cylinder
			mat_flags = 0;

			//--------------------------------------------------------------------
			// M_tt(theta_h) : Marschner marginal, longitudinal	scattering function(M)  
			//--------------------------------------------------------------------
			//calculate parameters for M
			float cos_theta_i = glm::dot(local_input_ray, normal);
			float sin_theta_i = glm::dot(local_input_ray, tangent);
			float theta_i = glm::acos(cos_theta_i);//Angle between input ray and fibers normal-plane
			float cos_theta_r = glm::dot(glm::normalize(local_output_ray), normal);
			float sin_theta_r = glm::dot(glm::normalize(local_output_ray), tangent);
			float theta_r = glm::acos(cos_theta_r);//Angle between reflected output ray and fibers normal-plane	
			float theta_h = (theta_r + theta_i) / 2.f;//theta half angle
			float theta_d = (theta_r - theta_i) / 2.f;//theta difference angle

			//M_tt(theta_h) -> gaussian function with zero-mean and our fibers material standart derivation.
			output_pdf = BSDFHelper::normal_gauss_pdf(theta_h - glm::radians(-alpha_r / 2.f), 0.0f, beta_r / 2.f);

			//--------------------------------------------------------------------
			// N_tt(phi) : Marschner conditional, azimuthal scattering function(N)
			//--------------------------------------------------------------------

			//calculate parameters for N
			glm::vec3 in_ray_normplane = glm::normalize(local_input_ray - sin_theta_i * tangent);//Input lightvector, projected onto the normal-plane
			glm::vec3 out_ray_normplane = glm::normalize(glm::normalize(local_output_ray) - sin_theta_r * tangent);;//Output lightvector, projected onto the normal-plane
			float phi = glm::acos(glm::min(1.0f, glm::dot(out_ray_normplane, in_ray_normplane)));

			//calculation of h from marschner but the concrete equation is from d'Eon paper (Above Equation 9)
			float a = 1 / material->m_ior;
			float nenner = glm::sqrt(1 + glm::pow2(a) - 2 * a * glm::sign(phi) * glm::sin(phi / 2.f));
			float h = (glm::sign(phi) * glm::cos(phi / 2.f)) / nenner;
			float gamma_i = glm::asin(h);

			//Bravais (virtual index of reflection eta_one, eta_two) calculation
			float cos_theta_d = glm::cos(theta_d);
			float x1 = glm::sqrt(glm::pow2(material->m_ior) - glm::pow2(glm::sin(theta_d)));
			float bravais_first = x1 / cos_theta_d;//first bravais index (virtual eta)
			float bravais_sec = glm::pow2(material->m_ior) * cos_theta_d / x1;//second bravais index (virtual eta)
			float c = glm::asin(1.f / bravais_first);

			//From ruenz12 bachelorthesis Equation 44, which derived it from marschner Equation 10.
			float dphi_dh = 1.f / (1 / glm::sqrt(1 - h * h) * ((-(24 * c / glm::pow3(glm::pi<float>())) * glm::pow2(gamma_i)) + (6 * c / glm::pi<float>() - 2)));

			//calculate fresnel for attenuation factor
			float fresnel = BSDFHelper::dialectricFresnel(glm::cos(gamma_i), bravais_first, bravais_sec);
			if (fresnel == 1) fresnel = 0.f;//so it doesnt change att_color to 0

			//helper variables for attenuation
			float cos_gamma_t = 2.f * glm::cos(glm::asin(h / bravais_first));
			glm::vec3 new_sigma = glm::vec3(material->fetchParameterColor<MatParamType::DIFFUSE>(texcoord)) / glm::cos(theta_r);//new color absorption coefficient
			//attenuation factor, which contains color absorbtion
			glm::vec3 att_color = glm::pow2(1 - fresnel) * glm::exp(new_sigma * cos_gamma_t);

			//final term for N_tt(phi). Marschner Equation 8.
			glm::vec3 n_tt = 0.5f * att_color / glm::abs(2 * dphi_dh);

			//--------------------------------------------------------------------

			//return final scattering function
			return (output_pdf * n_tt / glm::pow2(cos_theta_d));

		}//We are inside the cylinder and we want to take Cylinder Path: TRT (refractive transmission, p = 2)
		else if (mat_flags & BSDFHelper::MATFLAG_CYLINDER_TR_BOUNCE)//We are inside the cylinder and we want to take Cylinder Path: TRT (p = 2)
		{
			//refract on the first wall of the cylinder to get output ray
			local_output_ray = glm::refract(local_input_ray, glm::faceforward(normal, local_input_ray, normal), 1.f);

			//rotate towards normal to account for the tilted fiber surface
			local_output_ray = glm::vec3(glm::vec4(local_output_ray, 0.f) * glm::rotate(glm::radians(-3.f * alpha_r / 2.f), tangent));

			//reset material flags and set it to finished path since we leave the cylinder
			mat_flags = 0;

			//--------------------------------------------------------------------
			// M_trt(theta_h) : Marschner marginal, longitudinal	scattering function(M)  
			//--------------------------------------------------------------------

			//calculate parameters for M
			float cos_theta_i = glm::dot(local_input_ray, normal);
			float sin_theta_i = glm::dot(local_input_ray, tangent);
			float theta_i = glm::acos(cos_theta_i);//Angle between input ray and fibers normal-plane
			float cos_theta_r = glm::dot(glm::normalize(local_output_ray), normal);
			float sin_theta_r = glm::dot(glm::normalize(local_output_ray), tangent);
			float theta_r = glm::acos(cos_theta_r);//Angle between reflected output ray and fibers normal-plane	
			float theta_h = (theta_r + theta_i) / 2.f;//theta half angle
			float theta_d = (theta_r - theta_i) / 2.f;//theta difference angle
			sample.x = theta_i; sample.y = 0;//We store the theta_i angle in it, because we need it in the MarschnerHairShader method
											 //and this trt-path is the last call of the bsdf. We dont want to override the previously stored alpha/beta values.

			//M_trt(theta_h) -> gaussian function with zero-mean and our fibers material standart derivation.
			output_pdf = BSDFHelper::normal_gauss_pdf(theta_h - glm::radians(-3 * alpha_r / 2.f), 0.0f, 2.f * beta_r);

			//--------------------------------------------------------------------
			// N_trt(phi) : Marschner conditional, azimuthal scattering function(N) [Chapter 5.2.2]
			//--------------------------------------------------------------------

			//calculate parameters for N
			glm::vec3 in_ray_normplane = glm::normalize(local_input_ray - sin_theta_i * tangent);//Input lightvector, projected onto the normal-plane
			glm::vec3 out_ray_normplane = glm::normalize(glm::normalize(local_output_ray) - sin_theta_r * tangent);;//Output lightvector, projected onto the normal-plane
			float phi = glm::acos(glm::min(1.0f, glm::dot(out_ray_normplane, in_ray_normplane)));
			float w_c = 15.f;//azimuthal width of caustic. Between 10 and 25 degrees
			float k_g = 2.f;//glint scale factor. Between 0.5 and 5
			float t, phi_c, h;

			//help variables for bravais calculation
			float cos_theta_d = glm::cos(theta_d);
			float x1 = glm::sqrt(glm::pow2(material->m_ior) - glm::pow2(glm::sin(theta_d)));
			float bravais_first = x1 / cos_theta_d;//first bravais index (virtual eta)
			float bravais_sec = glm::pow2(material->m_ior) * cos_theta_d / x1;//second bravais index (virtual eta)
			float c = glm::asin(1.f / bravais_first);

			if (bravais_first < 2)//procedure from marschner paper chapter 5.2.2
			{
				float h_c = glm::sqrt((4 - glm::pow2(bravais_first)) / 3.f);//Marschner paper Equation 4
				phi_c = glm::radians(4 * glm::asin(h_c / bravais_first) - 2 * glm::asin(h_c) + 2 * glm::pi<float>());
				//From ruenz12 bachelorthesis Equation 44, which derived it from marschner Equation 10.
				float dphi_dh_h_c = 1.f / (1 / glm::sqrt(1 - h_c * h_c) * ((-(48 * c / glm::pow3(glm::pi<float>())) * glm::pow2(glm::asin(h_c)) + (12 * c / glm::pi<float>() - 2))));
				h = glm::min(0.5f, 2 * glm::sqrt(2 * w_c / glm::abs(dphi_dh_h_c)));
				t = 1;
			}
			else
			{
				phi_c = 0;
				h = 0.5f;
				t = glm::smoothstep(2.f, 2.3f, bravais_first);
			}

			float gamma_i = glm::asin(h);
			//From ruenz12 bachelorthesis Equation 44, which derived it from marschner Equation 10.
			float dphi_dh = 1.f / (1 / glm::sqrt(1 - h * h) * ((-(48 * c / glm::pow3(glm::pi<float>())) * glm::pow2(gamma_i) + (12 * c / glm::pi<float>() - 2))));

			//calculate fresnel part of attenuation
			float fresnel = BSDFHelper::dialectricFresnel(glm::cos(gamma_i), bravais_first, bravais_sec);
			if (fresnel == 1.f) fresnel = 0.f;
			float gamma_t = glm::asin(h / bravais_first);
			float cos_gamma_t = glm::cos(gamma_t);
			float fresnel_exit = BSDFHelper::dialectricFresnel(cos_gamma_t, 1 / bravais_first, 1 / bravais_sec);

			//new color absorption coefficient
			glm::vec3 new_sigma = glm::vec3(material->fetchParameterColor<MatParamType::DIFFUSE>(texcoord)) / glm::cos(theta_r);
			//full attenuation
			glm::vec3 att_color = glm::pow2(1 - fresnel) * fresnel_exit * glm::pow2(glm::exp(new_sigma * (2.f * cos_gamma_t)));

			//calculate gauss values for power apporximation from section 5.2.2
			float gauss_0 = BSDFHelper::normal_gauss_pdf(0.f, 0.f, w_c);
			float gauss_phi_diff = BSDFHelper::normal_gauss_pdf(phi - phi_c, 0.f, w_c);
			float gauss_phi_sum = BSDFHelper::normal_gauss_pdf(phi + phi_c, 0.f, w_c);

			//final term for N_trt(phi). Marschner Equation 8
			glm::vec3 n_trt = 0.5f * att_color / glm::abs(2 * dphi_dh);
			n_trt *= (1 - t * gauss_phi_diff / gauss_0);
			n_trt *= (1 - t * gauss_phi_sum / gauss_0);
			n_trt += t * k_g * att_color * h * (gauss_phi_diff + gauss_phi_sum);

			//--------------------------------------------------------------------

			//return final scattering function
			return (output_pdf * n_trt / glm::pow2(cos_theta_d));
		}
		else//We take Cylinder Path: R (Reflection on Surface)
		{
			//reflect input ray on the surface
			local_output_ray = glm::reflect(local_input_ray, glm::faceforward(normal, local_input_ray, normal));
			//rotate towards normal to account for the tilted fiber surface
			local_output_ray = glm::vec3(glm::vec4(local_output_ray, 0.f) * glm::rotate(glm::radians(2 * alpha_r), tangent));
			//Reset the mat_flags and set it to finished path since we leave the cylinder
			mat_flags = BSDFHelper::MATFLAG_SPECULAR_BOUNCE;

			//--------------------------------------------------------------------
			// M_r(theta_h) : Marschner marginal, longitudinal	scattering function(M)  			
			//--------------------------------------------------------------------

			//calculate parameters for M
			float cos_theta_i = glm::dot(local_input_ray, normal);
			float sin_theta_i = glm::dot(local_input_ray, tangent);
			float theta_i = glm::acos(cos_theta_i);//Angle between input ray and fibers normal-plane
			float cos_theta_r = glm::dot(glm::normalize(local_output_ray), normal);
			float sin_theta_r = glm::dot(glm::normalize(local_output_ray), tangent);
			float theta_r = glm::acos(cos_theta_r);//Angle between reflected output ray and fibers normal-plane	
			float theta_h = (theta_r + theta_i) / 2;//theta half angle
			float theta_d = (theta_r - theta_i) / 2;//theta difference angle			

			//M_r(theta_h) -> gaussian function with zero-mean and our fibers material standart derivation. TODO: Try Logistic function as suggested in pixar paper
			output_pdf = BSDFHelper::normal_gauss_pdf(theta_h - glm::radians(alpha_r), 0.0f, beta_r);

			//--------------------------------------------------------------------
			// N_r(phi) : Marschner conditional, azimuthal scattering function(N)
			//--------------------------------------------------------------------

			//calculate parameters for N
			glm::vec3 in_ray_normplane = glm::normalize(local_input_ray - sin_theta_i * tangent);//Input lightvector, projected onto the normal-plane
			glm::vec3 out_ray_normplane = glm::normalize(glm::normalize(local_output_ray) - sin_theta_r * tangent);;//Output lightvector, projected onto the normal-plane
			float phi = glm::acos(glm::min(1.0f, glm::dot(out_ray_normplane, in_ray_normplane)));
			float h = glm::sin(phi) * -0.5f;//root of the approximation for phi^ in marschner paper (Equation 10) with p = 0

			//From ruenz12 bachelorthesis Equation 44, which derived it from marschner Equation 10.
			float dh_dphi = -2.f / glm::sqrt(1 - h * h);

			//help variables for bravais calculation
			float cos_theta_d = glm::cos(theta_d);
			float x1 = glm::sqrt(glm::pow2(material->m_ior) - glm::pow2(glm::sin(theta_d)));
			//Bravais (virtual index of reflection eta_one, eta_two) calculation
			float bravais_first = x1 / cos_theta_d;//first bravais index (virtual eta)
			float bravais_sec = glm::pow2(material->m_ior) * cos_theta_d / x1;//second bravais index (virtual eta)

			//calculate attenuation factor with fresnel
			float fresnel = BSDFHelper::dialectricFresnel(glm::cos(glm::asin(h)), bravais_first, bravais_sec);

			//final term for N_r(phi). Marschner Equation 8
			float n_r = 0.5f * fresnel / glm::abs(2 * dh_dphi);

			//--------------------------------------------------------------------

			//Final value for the combined scattering function
			float scatt_func = output_pdf * n_r / glm::pow2(cos_theta_d);

			return glm::vec3(scatt_func);
		}
		return glm::vec3(.0f);
	}

	glm::vec3 MarschnerHairBSDF::evaluateLight(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray)
	{
		return glm::vec3(1);
	}

	////////////////////////////////////////////////////////////////////////////////////
	//////////////
	//////////////	D'EON HAIR BSDF
	//////////////
	////////////////////////////////////////////////////////////////////////////////////

	glm::vec3 DEonHairBSDF::localSample(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& normal, glm::vec2& sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance)
	{
		//get our object
		KIRK::Object *cylinder_obj = hit.m_object;

		//needed parameters
		Material* material = hit.m_object->getMaterial();
		glm::vec2 texcoord = hit.m_texcoord;

		//calculation of tangent vector
		glm::vec3 c1 = glm::cross(normal, glm::vec3(0.0, 0.0, 1.0));
		glm::vec3 c2 = glm::cross(normal, glm::vec3(0.0, 1.0, 0.0));
		glm::vec3 tangent = (glm::length(c1) > glm::length(c2)) ? glm::normalize(c1) : glm::normalize(c2);

		//Get lobe alpha and beta which is stored in the sample parameter
		float alpha_r = sample.x;//longitudinal shift. Suggested value from marschner hair paper between -10 and -5 degrees
		float beta_r = sample.y;//longitudinal width (stdev.). Suggested value from marschner hair paper between 5 and 10 degrees		

		//We are inside the cylinder and we want to take Cylinder Path: TT (refractive transmission, p = 1)
		if (mat_flags & BSDFHelper::MATFLAG_CYLINDER_T_BOUNCE)
		{
			//refract on the second wall of the cylinder to get output ray
			local_output_ray = glm::refract(local_input_ray, glm::faceforward(normal, local_input_ray, normal), 1.f / material->m_ior);

			//rotate towards normal to account for the tilted fiber surface
			local_output_ray = glm::vec3(glm::vec4(local_output_ray, 0.f) * glm::rotate(glm::radians(-alpha_r / 2), tangent));

			//reset material flags and set it to finished path since we leave the cylinder
			mat_flags = 0;

			//--------------------------------------------------------------------
			// M_tt(v, theta_i, theta_r) : d'Eon marginal, longitudinal	scattering function(M)  			
			//--------------------------------------------------------------------

			//calculate parameters for M
			float cos_theta_i = glm::dot(local_input_ray, normal);
			float sin_theta_i = glm::dot(local_input_ray, tangent);
			float theta_i = glm::acos(cos_theta_i);//Angle between input ray and fibers normal-plane
			float cos_theta_r = glm::dot(local_output_ray, normal);
			float sin_theta_r = glm::dot(glm::normalize(local_output_ray), tangent);
			float theta_r = glm::acos(cos_theta_r);//Angle between reflected output ray and fibers normal-plane	
			float theta_d = (theta_r - theta_i) / 2;//theta difference angle

			//M_tt -> d'Eon Equation 7
			output_pdf = M_p(glm::cos(-theta_i - glm::radians(alpha_r)), cos_theta_r,
				glm::sin(-theta_i - glm::radians(alpha_r)), sin_theta_r, glm::radians(beta_r * beta_r));

			//--------------------------------------------------------------------
			// N_tt(phi) : Marschner conditional, azimuthal scattering function(N)
			//--------------------------------------------------------------------

			//calculate parameters for N
			glm::vec3 in_ray_normplane = glm::normalize(local_input_ray - sin_theta_i * tangent);//Input lightvector, projected onto the normal-plane
			glm::vec3 out_ray_normplane = glm::normalize(glm::normalize(local_output_ray) - sin_theta_r * tangent);;//Output lightvector, projected onto the normal-plane
			float phi = glm::acos(glm::min(1.0f, glm::dot(out_ray_normplane, in_ray_normplane)));

			//calculation of h from marschner but the concrete equation is from d'Eon paper (Above Equation 9)
			float a = 1 / material->m_ior;
			float nenner = glm::sqrt(1 + glm::pow2(a) - 2 * a * glm::sign(phi) * glm::sin(phi / 2.f));
			float h = (glm::sign(phi) * glm::cos(phi / 2.f)) / nenner;
			float gamma_i = glm::asin(h);

			//help variables for bravais calculation
			float cos_theta_d = glm::cos(theta_d);
			float x1 = glm::sqrt(glm::pow2(material->m_ior) - glm::pow2(glm::sin(theta_d)));
			//Bravais (virtual index of reflection) calculation
			float bravais = x1 / cos_theta_d;

			//second part of N equation (first part is attenuation factor)
			//Gaussian detector with 20 iterations. d'Eon Equation 11
			float d_tt = 0;
			for (int i = -10; i <= 10; i++)
				d_tt += BSDFHelper::normal_gauss_pdf(phi - 2.f * glm::pi<float>() * i, 0.f, beta_r / 2);

			//calculate fresnel which is part of attenuation. d'Eon Equation 14
			float fresnel = BSDFHelper::dialectricFresnel(glm::acos(cos_theta_d * glm::cos(gamma_i)), material->m_ior, 1.f);
			if (fresnel == 1) fresnel = 0.f;//so it doesnt change att_color to 0

			//helper variables for attenuation
			float cos_2gamma_t = glm::cos(2 * glm::asin(h / bravais));
			glm::vec3 new_sigma = glm::vec3(material->fetchParameterColor<MatParamType::DIFFUSE>(texcoord)) / cos_theta_r;//new color absorption coefficent
			//full attenuation factor. d'Eon Equation 13 with p = 1
			glm::vec3 att_color = glm::pow2(1 - fresnel) * glm::exp(-2.f * new_sigma * (1 + cos_2gamma_t));

			//final term for N_tt(phi). d'Eon Equation 10
			glm::vec3 n_tt = 0.5f * att_color * d_tt;

			//--------------------------------------------------------------------

			//return final scattering function
			return (output_pdf * n_tt);

		}//We are inside the cylinder and we want to take Cylinder Path: TRT (refractive transmission, p = 2)
		else if (mat_flags & BSDFHelper::MATFLAG_CYLINDER_TR_BOUNCE)//We are inside the cylinder and we want to take Cylinder Path: TRT (p = 2)
		{
			//refract on the first wall of the cylinder to get output ray
			local_output_ray = glm::refract(local_input_ray, glm::faceforward(normal, local_input_ray, normal), 1.f);

			//rotate towards normal to account for the tilted fiber surface
			local_output_ray = glm::vec3(glm::vec4(local_output_ray, 0.f) * glm::rotate(glm::radians(-3.f * alpha_r / 2.f), tangent));

			//reset material flags and set it to finished path since we leave the cylinder
			mat_flags = 0;

			//--------------------------------------------------------------------
			// M_trt(v, theta_i, theta_r) : d'Eon marginal, longitudinal scattering function(M)  			
			//--------------------------------------------------------------------

			//calculate parameters for M
			float cos_theta_i = glm::dot(local_input_ray, normal);
			float sin_theta_i = glm::dot(local_input_ray, tangent);
			float theta_i = glm::acos(cos_theta_i);//Angle between input ray and fibers normal-plane
			float cos_theta_r = glm::dot(local_output_ray, normal);
			float sin_theta_r = glm::dot(glm::normalize(local_output_ray), tangent);
			float theta_r = glm::acos(cos_theta_r);//Angle between reflected output ray and fibers normal-plane	
			float theta_d = (theta_r - theta_i) / 2;//theta difference angle
			sample.x = theta_i; sample.y = 0;//We store the theta_i angle in it, because we need it in the MarschnerHairShader method
											 //and this trt-path is the last call of the bsdf. We dont want to override the previously stored alpha/beta values.

			//M_trt -> d'Eon Equation 7
			output_pdf = M_p(glm::cos(-theta_i - glm::radians(4 * alpha_r)), cos_theta_r,
				glm::sin(-theta_i - glm::radians(4 * alpha_r)), sin_theta_r, glm::radians(beta_r * beta_r));

			//--------------------------------------------------------------------
			// N_trt(phi) : Marschner conditional, azimuthal scattering function(N)
			//--------------------------------------------------------------------

			//calculate parameters for N
			glm::vec3 in_ray_normplane = glm::normalize(local_input_ray - sin_theta_i * tangent);//Input lightvector, projected onto the normal-plane
			glm::vec3 out_ray_normplane = glm::normalize(glm::normalize(local_output_ray) - sin_theta_r * tangent);;//Output lightvector, projected onto the normal-plane
			float phi = glm::acos(glm::min(1.0f, glm::dot(out_ray_normplane, in_ray_normplane)));
			float gamma_i = glm::angle(local_input_ray, glm::normalize(normal));//angle between input ray and surface normal in radians
			float h = glm::sin(gamma_i);

			//help variables for bravais calculation
			float cos_theta_d = glm::cos(theta_d);
			float x1 = glm::sqrt(glm::pow2(material->m_ior) - glm::pow2(glm::sin(theta_d)));
			//Bravais (virtual index of reflection) calculation
			float bravais = x1 / cos_theta_d;

			//second part of N equation (first part is attenuation factor)
			//Gaussian detector with 20 iterations. d'Eon Equation 11
			float d_trt = 0;
			for (int i = -10; i <= 10; i++)
				d_trt += BSDFHelper::normal_gauss_pdf(phi - 2.f * glm::pi<float>() * i, 0.f, beta_r * 2);

			//calculate fresnel which is part of attenuation. d'Eon Equation 14
			float fresnel = BSDFHelper::dialectricFresnel(glm::acos(cos_theta_d * glm::cos(gamma_i)), material->m_ior, 1.f);

			//helper variables for attenuation
			float cos_2gamma_t = glm::cos(2 * glm::asin(h / bravais));
			glm::vec3 new_sigma = glm::vec3(material->fetchParameterColor<MatParamType::DIFFUSE>(texcoord)) / cos_theta_r;//new color absorption coefficent
			//full attenuation factor. d'Eon Equation 13 with p = 2
			glm::vec3 att_color = glm::pow2(1 - fresnel) * fresnel * glm::pow2(glm::exp(-2.f * new_sigma * (1 + cos_2gamma_t)));

			//final term for N_trt(phi). d'Eon Equation 10
			glm::vec3 n_trt = 0.5f * att_color * d_trt;

			//--------------------------------------------------------------------

			//return final scattering function
			return (output_pdf * n_trt);

		}
		else//We take Cylinder Path: R (Reflection on Surface)
		{
			//reflect input ray on the surface
			local_output_ray = glm::reflect(local_input_ray, glm::faceforward(normal, local_input_ray, normal));
			//rotate towards normal to account for the tilted fiber surface
			local_output_ray = glm::vec3(glm::vec4(local_output_ray, 0.f) * glm::rotate(glm::radians(2 * alpha_r), tangent));
			//Reset the mat_flags and set it to finished path since we leave the cylinder
			mat_flags = BSDFHelper::MATFLAG_SPECULAR_BOUNCE;

			//--------------------------------------------------------------------
			// M_r(v, theta_i, theta_r) : d'Eon marginal, longitudinal	scattering function(M)  			
			//--------------------------------------------------------------------

			//calculate parameters for M
			float cos_theta_i = glm::dot(local_input_ray, normal);
			float sin_theta_i = glm::dot(local_input_ray, tangent);
			float theta_i = glm::acos(cos_theta_i);//Angle between input ray and fibers normal-plane
			float cos_theta_r = glm::dot(local_output_ray, normal);
			float sin_theta_r = glm::dot(glm::normalize(local_output_ray), tangent);
			float theta_r = glm::acos(cos_theta_r);//Angle between reflected output ray and fibers normal-plane
			//M_r -> d'Eon Equation 7
			output_pdf = M_p(glm::cos(-theta_i + 2 * alpha_r), cos_theta_r, glm::sin(-theta_i + 2 * alpha_r), sin_theta_r, glm::radians(beta_r * beta_r));

			//--------------------------------------------------------------------
			// N_r(phi) : d'Eon conditional, azimuthal scattering function(N)
			//--------------------------------------------------------------------

			//calculate parameters for N
			glm::vec3 in_ray_normplane = glm::normalize(local_input_ray - sin_theta_i * tangent);//Input lightvector, projected onto the normal-plane
			glm::vec3 out_ray_normplane = glm::normalize(glm::normalize(local_output_ray) - sin_theta_r * tangent);;//Output lightvector, projected onto the normal-plane
			float phi = glm::acos(glm::min(1.0f, glm::dot(out_ray_normplane, in_ray_normplane)));
			float d_r = 0.25f * glm::abs(glm::cos(phi / 2.f));//d'Eon Equation 6

			//calculate attenuation factor with fresnel. d'Eon Equation 12
			float fresnel = BSDFHelper::dialectricFresnel(0.5f * glm::acos(glm::dot(local_input_ray, glm::normalize(local_output_ray))), 1.f, material->m_ior);

			//final term for N_r(phi). d'Eon Equation 10
			float n_r = 0.5f * fresnel * d_r;

			//--------------------------------------------------------------------

			//Final value for the combined scattering function
			return glm::vec3(output_pdf * n_r);
		}
		return glm::vec3(.0f);
	}

	glm::vec3 DEonHairBSDF::evaluateLight(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray)
	{
		return glm::vec3(1);
	}

	float DEonHairBSDF::M_p(float cos_theta_c, float cos_theta_r, float sin_theta_c, float sin_theta_r, float v)
	{
		//Calculation of M_p -> d'Eon Equation 7
		float csch = 2 * v * glm::sinh(1.f / v);//Cosecant hyperbolic function with 2v factor
		float x = cos_theta_c * cos_theta_r / v;//parameter for bessel calculation
		float y = sin_theta_c * sin_theta_r / v;//parameter for e function calculation
		float m_p = (v <= 0.1f) ? glm::exp(logBessel(x) - y - (1 / v) + 0.6931f + glm::log(1 / (2 * v))) : (glm::exp(-y) * DEonHairBSDF::bessel(x)) / csch;//compensation for low roughness v values from d'eon 2013 paper
		return m_p;
	}

	inline float DEonHairBSDF::bessel(float x)
	{
		float val = 0;
		float x2i = 1;
		int ifact = 1;
		int i4 = 1;
		// I0(x) \approx Sum_i x^(2i) / (4^i (i!)^2)
		for (int i = 0; i < 10; ++i) {
			if (i > 1) ifact *= i;
			val += x2i / (i4 * (ifact * ifact));
			x2i *= x * x;
			i4 *= 4;
		}
		return val;
	}

	inline float DEonHairBSDF::logBessel(float x)
	{
		if (x > 12)
			return x + 0.5 * (-std::log(2 * glm::pi<float>()) + std::log(1 / x) + 1 / (8 * x));
		else
			return std::log(_j0(x));
	}
}
