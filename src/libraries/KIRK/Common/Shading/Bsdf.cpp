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
			if (t<-1) t = -1;
			if (t> 1) t = 1;
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
			if (r<0) A = -A;
			B = (A == 0 ? 0 : B = q / A);

			a /= 3;
			x[0] = (A + B) - a;
			x[1] = -0.5*(A + B) - a;
			x[2] = 0.5*sqrt(3.)*(A - B);
			if (fabs(x[2])<eps) { x[2] = x[1]; return(2); }
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

	float BSDFHelper::normal_gauss_pdf(float x, float mean, float stdder)
	{
		static const float inv_sqrt_2pi = 0.3989422804014327;
		float a = (x - mean) / stdder;

		return inv_sqrt_2pi / stdder * std::exp(-0.5f * a * a);
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

	glm::vec3 BSDF::sample(const Intersection& hit, const glm::vec3& ray_in, const glm::vec3& normal, glm::vec3& ray_out, float& pdf, int& mat_flags, const glm::vec2& sampling, bool useRadianceOverImportance) const
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

	glm::vec3 MarschnerHairBSDF::localSample(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& normal, glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance)
	{
		//cast object to cylinder
		KIRK::Cylinder *cylinder_obj = dynamic_cast<KIRK::Cylinder*>(hit.m_object);
		//if we have no cylinder as object we cant use marschner hair
		if (cylinder_obj == NULL)
			return glm::vec3(0.0f);

		//needed parameters
		Material* material = hit.m_object->getMaterial();
		glm::vec2 texcoord = hit.m_texcoord;
		float alpha_r = material->fetchParameterFloat<MatParamType::ALPHA_LOBESHIFT>(texcoord);//longitudinal shift: R lobe in radians
		float beta_r = material->fetchParameterFloat<MatParamType::BETA_LOBEWIDTH>(texcoord);//longitudinal width (stdev.): R lobe in radians
		glm::vec3 norm_in_ray = glm::normalize(local_input_ray);//normalized input ray
		//move input ray to cylinders local space
		/*our v-axis(vector through the center of the cylinder(apexpoint - basepoint)) is u-axis in marschner model */
		glm::vec3 in_cyl_space = Math::worldToLocal(local_input_ray, cylinder_obj->getV(), cylinder_obj->getU(), cylinder_obj->getW());

		

		//We are inside the cylinder and we want to take Cylinder Path: TT (refractive transmission, p = 1)
		if ((mat_flags & BSDFHelper::MATFLAG_CYLINDER_T_BOUNCE) && !(mat_flags & BSDFHelper::MATFLAG_CYLINDER_TR_BOUNCE))
		{
			//refract on the second wall of the cylinder to get output ray
			local_output_ray = glm::refract(-norm_in_ray, glm::faceforward(normal, -norm_in_ray, normal), material->m_ior);

			//rotate towards normal to account for the tilted fiber surface
			local_output_ray = glm::vec3(glm::vec4(local_output_ray, 0.f) * glm::rotate(-alpha_r / 2, cylinder_obj->getV()));

			//reset material flags since we leave the cylinder
			mat_flags = 0;

			//move output ray to cylinders local space
			/*our v-axis(vector through the center of the cylinder(apexpoint - basepoint)) is u-axis in marschner model */
			glm::vec3 out_cyl_space = Math::worldToLocal(local_output_ray, cylinder_obj->getV(), cylinder_obj->getU(), cylinder_obj->getW());

			//--------------------------------------------------------------------
			// M_tt(theta_h) : Marschner marginal, longitudinal	scattering function(M)  
			//--------------------------------------------------------------------
			//calculate parameters for M
			//float theta_i = glm::asin(glm::abs(glm::dot(cylinder_obj->getV(), norm_in_ray)) / glm::length(local_input_ray));//Angle between input ray and fibers u-w(normal-) plane (In Marschner paper the tangent of the fiber is u-axis but in our cylinder class its the v-axis)
			//float theta_r = glm::asin(glm::abs(glm::dot(cylinder_obj->getV(), local_output_ray)) / glm::length(local_output_ray));//Angle between reflected output ray and fibers u-w(normal-) plane

			float theta_i = std::atan2(std::hypot(in_cyl_space.x, in_cyl_space.z), in_cyl_space.y);//Angle between input ray and fibers u-w(normal-) plane (In Marschner paper the tangent of the fiber is u-axis but in our cylinder class its the v-axis)
			float theta_r = std::atan2(std::hypot(out_cyl_space.x, out_cyl_space.z), out_cyl_space.y);//Angle between reflected output ray and fibers u-w(normal-) plane
			float theta_h = (theta_r + theta_i) / 2.f;//theta half angle
			float theta_d = (theta_r - theta_i) / 2.f;//theta difference angle
			float gaussian_x = theta_h - (-alpha_r / 2.f);//longitudinal shift for TT-Path is smaller than for R-Path

			//M_tt(theta_h) -> gaussian function with zero-mean and our fibers material standart derivation.
			output_pdf = 10.f * BSDFHelper::normal_gauss_pdf(gaussian_x, 0.0f, glm::degrees(beta_r / 2.f));

			//--------------------------------------------------------------------
			// N_tt(phi) : Marschner conditional, azimuthal scattering function(N)
			//--------------------------------------------------------------------

			//calculate parameters for N
			float phi_i = std::atan2(in_cyl_space.x, in_cyl_space.y);
			float phi_r = std::atan2(out_cyl_space.x, out_cyl_space.y);
			float bravais = glm::cos(theta_d) * glm::pow(glm::pow2(material->m_ior) - glm::pow2(glm::sin(theta_d)), -0.5f);//Bravais index
			float c = glm::asin(1 / bravais);

			/* Since we have the normal of the intersection we can calculate gamma_i instead of the approximation that
			*  marschner has done with the cubic function
			*
			//parameters for the cubic function delta_phi(1, gamma_i) in marschner (Equation 10)
			double pi_pow_3 = glm::pow3(std::_Pi);
			double a1 = -8 * c / pi_pow_3;
			double a2 = 6 * c / std::_Pi - 2;
			double a3 = std::_Pi - (phi_r - phi_i);
			double x[3];//result of delta_phi(1, gamma_i)
			BSDFHelper::SolveP3(x, 0, a2 / a1, a3 / a1);
			*/

			float gamma_i = glm::angle(norm_in_ray, glm::normalize(normal));//angle between input ray and surface normal in radians
			float h = glm::sin(gamma_i);
			//second part of N equation (first part is attenuation factor)
			float dh_dphi = 1 / glm::abs((1 / glm::sqrt(1 - h * h)) * ((-(24 * c / glm::pow3(std::_Pi)) * glm::pow2(gamma_i)) + (6 * c / std::_Pi - 2)));

			//calculate fresnel for attenuation factor
			bool entering = glm::dot(local_input_ray, normal) > 0.f;
			float eta_i = entering ? 1.f : material->m_ior;
			float eta_t = entering ? material->m_ior : 1.f;
			float fresnel = BSDFHelper::dialectricFresnel(glm::abs(glm::dot(norm_in_ray, normal)), eta_i, eta_t);
			//attenuation factor, which contains color absorbtion
			float cos_gamma_t = -2.f * glm::cos(glm::asin(h / bravais));
			glm::vec3 att_color = glm::pow2(1 - fresnel) * glm::exp(glm::vec3(material->fetchParameterColor<MatParamType::DIFFUSE>(texcoord)) * cos_gamma_t);

			//final term for N_tt(phi)
			glm::vec3 n_tt = 0.5f * att_color * dh_dphi;

			//--------------------------------------------------------------------

			//return final scattering function
			return (output_pdf * n_tt / glm::pow2(glm::cos(theta_d)));

		}//We are inside the cylinder and we want to take Cylinder Path: TRT (refractive transmission, p = 2)
		else if(mat_flags & BSDFHelper::MATFLAG_CYLINDER_TR_BOUNCE)
		{
			//If we have only made T-Path so far, we reflect on the second cylinder wall to get the TR-Path
			if (!(mat_flags & BSDFHelper::MATFLAG_CYLINDER_T_BOUNCE))
			{
				//Reflect on the second cylinder wall 
				local_output_ray = glm::reflect(-norm_in_ray, faceforward(normal, -norm_in_ray, normal));
				//set material flag, so we know we already have the TR-Path in the TRT-Section
				mat_flags |= BSDFHelper::MATFLAG_CYLINDER_T_BOUNCE;
				return glm::vec3(0);
			}
			else//We already have made the TR-Path. Now we have the final intersection for the TRT-Path
			{
				//refract on the first wall of the cylinder to get output ray
				local_output_ray = glm::refract(-norm_in_ray, glm::faceforward(normal, -norm_in_ray, normal), material->m_ior);

				//rotate towards normal to account for the tilted fiber surface
				local_output_ray = glm::vec3(glm::vec4(local_output_ray, 0.f) * glm::rotate(-3.f * alpha_r / 2.f, cylinder_obj->getV()));

				//reset material flags since we leave the cylinder
				mat_flags = 0;

				//move output ray to cylinders local space
				/*our v-axis(vector through the center of the cylinder(apexpoint - basepoint)) is u-axis in marschner model  */
				glm::vec3 out_cyl_space = Math::worldToLocal(local_output_ray, cylinder_obj->getV(), cylinder_obj->getU(), cylinder_obj->getW());

				//--------------------------------------------------------------------
				// M_trt(theta_h) : Marschner marginal, longitudinal	scattering function(M)  
				//--------------------------------------------------------------------

				//calculate parameters for M
				//float theta_i = glm::asin(glm::abs(glm::dot(cylinder_obj->getV(), norm_in_ray)) / glm::length(local_input_ray));//Angle between input ray and fibers u-w(normal-) plane (In Marschner paper the tangent of the fiber is u-axis but in our cylinder class its the v-axis)
				//float theta_r = glm::asin(glm::abs(glm::dot(cylinder_obj->getV(), local_output_ray)) / glm::length(local_output_ray));//Angle between reflected output ray and fibers u-w(normal-) plane

				float theta_i = std::atan2(std::hypot(in_cyl_space.x, in_cyl_space.z), in_cyl_space.y);//Angle between input ray and fibers u-w(normal-) plane (In Marschner paper the tangent of the fiber is u-axis but in our cylinder class its the v-axis)
				float theta_r = std::atan2(std::hypot(out_cyl_space.x, out_cyl_space.z), out_cyl_space.y);//Angle between reflected output ray and fibers u-w(normal-) plane
				float theta_h = (theta_r + theta_i) / 2.f;//theta half angle
				float theta_d = (theta_r - theta_i) / 2.f;//theta difference angle
				float gaussian_x = theta_h - (-3 * alpha_r / 2.f);//longitudinal shift for TT-Path is smaller than for R-Path

				//M_trt(theta_h) -> gaussian function with zero-mean and our fibers material standart derivation.
				output_pdf = 10.f * BSDFHelper::normal_gauss_pdf(gaussian_x, 0.0f, glm::degrees(2.f * beta_r));

				//--------------------------------------------------------------------
				// N_trt(phi) : Marschner conditional, azimuthal scattering function(N)
				//--------------------------------------------------------------------

				//calculate parameters for N
				float phi_i = std::atan2(in_cyl_space.x, in_cyl_space.y);
				float phi_r = std::atan2(out_cyl_space.x, out_cyl_space.y);
				float bravais = glm::cos(theta_d) * glm::pow(glm::pow2(material->m_ior) - glm::pow2(glm::sin(theta_d)), -0.5f);//Bravais index
				float c = glm::asin(1 / bravais);

				/* Since we have the normal of the intersection we can calculate gamma_i instead of the approximation that
				*  marschner has done with the cubic function
				*
				//parameters for the cubic function delta_phi(1, gamma_i) in marschner (Equation 10)
				double pi_pow_3 = glm::pow3(std::_Pi);
				double a1 = -16 * c / pi_pow_3;
				double a2 = 12 * c / std::_Pi - 2;
				double a3 = std::_Pi - (phi_r - phi_i);
				double x[3];//result of delta_phi(1, gamma_i)
				BSDFHelper::SolveP3(x, 0, a2 / a1, a3 / a1);
				*/

				float gamma_i = glm::angle(norm_in_ray, glm::normalize(normal));//angle between input ray and surface normal in radians
				float h = glm::sin(gamma_i);
				//second part of N equation (first part is attenuation factor)
				float dh_dphi = 1 / glm::abs((1 / glm::sqrt(1 - h * h)) * ((-(48 * c / glm::pow3(std::_Pi)) * glm::pow2(gamma_i) + (12 * c / std::_Pi - 2))));

				//calculate fresnel for attenuation
				bool entering = glm::dot(local_input_ray, normal) > 0.f;
				float eta_i = entering ? 1.f : material->m_ior;
				float eta_t = entering ? material->m_ior : 1.f;
				float fresnel = BSDFHelper::dialectricFresnel(glm::abs(glm::dot(norm_in_ray, normal)), eta_i, eta_t);
				//attenuation factor, which contains color absorbtion
				float gamma_t = glm::cos(glm::asin(h / bravais));
				float cos_gamma_t = -2 * gamma_t;
				float fresnel_exit = BSDFHelper::dialectricFresnel(gamma_t, 1 / eta_i, 1 / eta_t);
				glm::vec3 att_color = glm::pow2(1 - fresnel) * glm::pow2(fresnel_exit) * glm::pow2(glm::exp(glm::vec3(material->fetchParameterColor<MatParamType::DIFFUSE>(texcoord)) * cos_gamma_t));

				//final term for N_trt(phi)
				glm::vec3 n_trt = 0.5f * att_color * dh_dphi;

				//--------------------------------------------------------------------

				//return final scattering function
				return (output_pdf * n_trt / glm::pow2(glm::cos(theta_d)));
			}
		}
		else//We hit the cylinder for the first time and decide which Cylinder Path we want to take
		{
			//random value generator to select our Cylinder Path randomly. p = 0 is R-Path, p = 1 is TT-Path, p = 2 is TRT-Path
			std::random_device rd;
			std::mt19937 mt(rd());
			std::uniform_int_distribution<int> dist(0, 2);//std::nextafter so we get the [0,1] interval instead of [0,1)
			//int p = dist(mt);
			int p = 1;

			//We take Cylinder Path: R (surface reflection)
			if (p == 0)
			{
				//reflect input ray on the surface
				local_output_ray = glm::reflect(-norm_in_ray, glm::faceforward(normal, -norm_in_ray, normal));
				//rotate towards normal to account for the tilted fiber surface
				local_output_ray = glm::vec3(glm::vec4(local_output_ray, 0.f) * glm::rotate(alpha_r, cylinder_obj->getV()));
				//Reset the mat_flags since we leave the cylinder
				mat_flags = 0;

				//move output ray to cylinders local space
				/*our v-axis(vector through the center of the cylinder(apexpoint - basepoint)) is u-axis in marschner model  */
				glm::vec3 out_cyl_space = Math::worldToLocal(local_output_ray, cylinder_obj->getV(), cylinder_obj->getU(), cylinder_obj->getW());

				//--------------------------------------------------------------------
				// M_r(theta_h) : Marschner marginal, longitudinal	scattering function(M)  			
				//--------------------------------------------------------------------

				//calculate parameters for M

				//float theta_i = glm::asin(glm::abs(glm::dot(cylinder_obj->getV(), norm_in_ray)) / glm::length(local_input_ray));//Angle between input ray and fibers u-w(normal-) plane (In Marschner paper the tangent of the fiber is u-axis but in our cylinder class its the v-axis)
				//float theta_r = glm::asin(glm::abs(glm::dot(cylinder_obj->getV(), local_output_ray)) / glm::length(local_output_ray));//Angle between reflected output ray and fibers u-w(normal-) plane

				float theta_i = std::atan2(std::hypot(in_cyl_space.x, in_cyl_space.z), in_cyl_space.y);//Angle between input ray and fibers u-w(normal-) plane (In Marschner paper the tangent of the fiber is u-axis but in our cylinder class its the v-axis)
				float theta_r = std::atan2(std::hypot(out_cyl_space.x, out_cyl_space.z), out_cyl_space.y);//Angle between reflected output ray and fibers u-w(normal-) plane
				float theta_h = (theta_r + theta_i) / 2;//theta half angle
				float theta_d = (theta_r - theta_i) / 2;//theta difference angle
				float gaussian_x = theta_h - alpha_r;

				//M_r(theta_h) -> gaussian function with zero-mean and our fibers material standart derivation. TODO: Try Logistic function as suggested in pixar paper
				output_pdf = 10.f * BSDFHelper::normal_gauss_pdf(gaussian_x, 0.0f, glm::degrees(beta_r));

				//--------------------------------------------------------------------
				// N_r(phi) : Marschner conditional, azimuthal scattering function(N)
				//--------------------------------------------------------------------

				//calculate parameters for N

				/* Since we have the normal of the intersection we can calculate gamma_i instead of the approximation that
				*  marschner has done with the cubic function
				*
				//float phi_i = std::atan2(in_cyl_space.x, in_cyl_space.y);
				//float phi_r = std::atan2(out_cyl_space.x, out_cyl_space.y);
				//float h = glm::sin(phi_r - phi_i) * -0.5f;//root of the approximation for phi^ in marschner paper (Equation 10) with p = 0
				*/
				float gamma_i = glm::angle(norm_in_ray, glm::normalize(normal));//angle between input ray and surface normal in radians
				float h = glm::sin(gamma_i);
				float dh_dphi = glm::abs(-2 / glm::sqrt(1 - h * h));

				//calculate attenuation factor with fresnel
				bool entering = glm::dot(local_input_ray, normal) > 0.f;
				float eta_i = entering ? 1.f : material->m_ior;
				float eta_t = entering ? material->m_ior : 1.f;
				float fresnel = BSDFHelper::dialectricFresnel(glm::abs(glm::dot(norm_in_ray, normal)), eta_i, eta_t);

				//final term for N_r(phi)
				float n_r = 0.5f * fresnel * dh_dphi;

				//--------------------------------------------------------------------

				//Final value for the combined scattering function
				float scatt_func = output_pdf * n_r / glm::pow2(glm::cos(theta_d));

				return glm::vec3(scatt_func);

			}//We take Cylinder Path: TT (refractive transmission)
			else if (p == 1)
			{
				/*
				 * We calculate the refraction on the first cylinder wall (after this we have a T-Path). 
				 * Because we dont know the intersection point on the second wall of the cylinder, we cant calculate 
				 * the components for the scattering function yet. We only set the local_output_ray so it can be traced towards the
				 * TT-Path in the next call of this function.
				*/
				//refract input ray on surface
				local_output_ray = glm::refract(-norm_in_ray, glm::faceforward(normal, -norm_in_ray, normal), 1.0f / material->m_ior);
				//set the material flags to TT-Path
				mat_flags = BSDFHelper::MATFLAG_CYLINDER_T_BOUNCE;
				return glm::vec3(0);
			}//We take Cylinder Path: TRT (internal reflection)
			else if (p == 2)
			{
				/*
				* We calculate the refraction on the first cylinder wall (after this we have a T-Path).
				* Because we dont know the intersection point on the second wall of the cylinder, we cant calculate
				* the components for the scattering function yet. We only set the local_output_ray so it can be traced towards the
				* TRT-Path in the next call of this function.
				*/
				//refract input ray on surface
				local_output_ray = glm::refract(-norm_in_ray, glm::faceforward(normal, -norm_in_ray, normal), 1.0f / material->m_ior);
				//set the material flags to TRT-Path. 
				//When only the TR_BOUNCE is active, we know that we want to take the TRT-Path but we havent done the T-Section of this Path yet.
				mat_flags = BSDFHelper::MATFLAG_CYLINDER_TR_BOUNCE;			
				return glm::vec3(0);
			}
		}
		return glm::vec3(.0f);
	}

	glm::vec3 MarschnerHairBSDF::evaluateLight(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray) { return glm::vec3(0);/*todo: implement*/ }
}
