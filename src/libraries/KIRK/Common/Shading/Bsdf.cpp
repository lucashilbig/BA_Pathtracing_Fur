#include "Bsdf.h"
#include "KIRK/Common/Intersection.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <numeric>

namespace KIRK {
	////////////////////////////////////////////////////////////////////////////////////
	//////////////
	//////////////	BSDF HELPER FUNCTIONS
	//////////////
	////////////////////////////////////////////////////////////////////////////////////

	//Static const member initilization
	const float BSDFHelper::SQRT_PI_OVER_8 = 0.626657069f;

	//Member Functions
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

	float BSDFHelper::atan2ThetaToMarschner(const float theta)
	{
		float pi = glm::pi<float>();
		if (theta > pi / 2.f)
			return (pi - theta);
		else if (theta < -pi / 2.f)
			return (-pi - theta);
		else
			return theta;
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

	float BSDFHelper::bessel(float x)
	{
		float val = 0;
		float x2i = 1;
		long ifact = 1;
		long i4 = 1;
		// I0(x) \approx Sum_i x^(2i) / (4^i (i!)^2)
		for (int i = 0; i < 5; ++i) {
			if (i > 1) ifact *= i;
			val += x2i / (i4 * (ifact * ifact));
			x2i *= x * x;
			i4 *= 4;
		}
		return val;
	}

	float BSDFHelper::logBessel(float x)
	{
		if (x > 12)
			return x + 0.5 * (-std::log(2 * glm::pi<float>()) + std::log(1 / x) + 1 / (8 * x));
		else
			return std::log(bessel(x));
	}

	float BSDFHelper::Logistic(float x, float s) {
		x = std::abs(x);
		return glm::exp(-x / s) / (s * glm::pow2(1 + glm::exp(-x / s)));
	}

	float BSDFHelper::LogisticCDF(float x, float s) {
		return 1 / (1 + std::exp(-x / s));
	}

	float BSDFHelper::TrimmedLogistic(float x, float s, float a, float b) {
		return Logistic(x, s) / (LogisticCDF(b, s) - LogisticCDF(a, s));
	}

	float BSDFHelper::SampleTrimmedLogistic(float u, float s, float a, float b) {
		float k = BSDFHelper::LogisticCDF(b, s) - BSDFHelper::LogisticCDF(a, s);
		float x = -s * std::log(1 / (u * k + BSDFHelper::LogisticCDF(a, s)) - 1);
		return glm::clamp(x, a, b);
	}

	float BSDFHelper::Phi(int p, float gammaO, float gammaT) {
		return 2 * p * gammaT - 2 * gammaO + p * glm::pi<float>();
	}

	glm::vec2 BSDFHelper::DemuxFloat(float f) {
		uint64_t v = f * (1ull << 32);
		uint32_t bits[2] = { BSDFHelper::Compact1By1(v),  BSDFHelper::Compact1By1(v >> 1) };
		return { bits[0] / float(1 << 16), bits[1] / float(1 << 16) };
	}

	uint32_t BSDFHelper::Compact1By1(uint32_t x) {
		// TODO: as of Haswell, the PEXT instruction could do all this in a
		// single instruction.
		// x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
		x &= 0x55555555;
		// x = --fe --dc --ba --98 --76 --54 --32 --10
		x = (x ^ (x >> 1)) & 0x33333333;
		// x = ---- fedc ---- ba98 ---- 7654 ---- 3210
		x = (x ^ (x >> 2)) & 0x0f0f0f0f;
		// x = ---- ---- fedc ba98 ---- ---- 7654 3210
		x = (x ^ (x >> 4)) & 0x00ff00ff;
		// x = ---- ---- ---- ---- fedc ba98 7654 3210
		x = (x ^ (x >> 8)) & 0x0000ffff;
		return x;
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

	glm::vec3 SpecularReflectionBSDF::evaluateLight(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray) 
	{
		/*bool reflect = glm::dot(local_input_ray, hit.m_normal) * glm::dot(local_output_ray, hit.m_normal) > 0;
		if (reflect) { return glm::vec3(hit.m_object->getMaterial()->fetchParameterColor <MatParamType::SPECULAR>(hit.m_texcoord)) * glm::one_over_pi <float>(); }*/
		return glm::vec3(0);
	}

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
		/* Not need atm, because we calculate cylindric uvw axis for triangle in the cpu_scene convertion method
		//Cast Object to cylinder
		KIRK::Cylinder *cylinder_obj = dynamic_cast<KIRK::Cylinder*>(hit.m_object);		
		
		if (cylinder_obj == NULL)
		{
			//calculation of tangent vector, which will be the local u-axis
			glm::vec3 c1 = glm::cross(normal, glm::vec3(0.0, 0.0, 1.0));
			glm::vec3 c2 = glm::cross(normal, glm::vec3(0.0, 1.0, 0.0));

			u_axis = (glm::length(c1) > glm::length(c2)) ? glm::normalize(c1) : glm::normalize(c2);
			v_axis = normal;
			w_axis = glm::normalize(glm::cross(u_axis, v_axis));
		}
		else
		{
			//  Our v-axis(vector through the center of the cylinder(apexpoint - basepoint)) is u-axis in marschner model
			//and we used the local_output_ray to store our ray towards the light source we want to sample
			u_axis = cylinder_obj->getV();
			v_axis = cylinder_obj->getU();
			w_axis = cylinder_obj->getW();
		}*/

		//objects local axis
		glm::vec3 u_axis = hit.m_object->getV();
		glm::vec3 v_axis = hit.m_object->getU();
		glm::vec3 w_axis = hit.m_object->getW();

		//needed parameters
		Material* material = hit.m_object->getMaterial();
		glm::vec2 texcoord = hit.m_texcoord;

		//If we accidently hit an object without hair bsdf we can return
		if (material->m_current_bsdf < 6 || material->m_current_bsdf > 8)
			return glm::vec3(0);

		//move light ray to cylinders (fibers) local space.		
		glm::vec3 hit_to_light = Math::worldToLocal(glm::normalize(local_output_ray), u_axis, v_axis, w_axis);

		//calculate spherical coordinates of light ray. NOTE: std::atan2 is atan(y/x) using the signs of arguments to determine the correct quadrant. range [-pi, pi], so we wrap it to marschner range [-pi/2, pi/2]
		float theta_i = BSDFHelper::atan2ThetaToMarschner(std::atan2(hit_to_light.x, hit_to_light.z));//Angle between light ray and fibers v-w(normal-)plane(Inclinations with respect of the fibers normalplane). Around v-axis
		float phi_i = std::atan2(std::hypot(hit_to_light.x, hit_to_light.z), hit_to_light.y);//azimuth angle. Around x(fibers u-)-axis
		phi_i = (hit_to_light.z >= 0.0f) ? phi_i : -phi_i;//wrap phi to account for std::hypot() erasing minus of z(local w)-coordinate, because marschner phi range is [-pi, pi]

		//Get lobe alpha and beta which is stored in the sample parameter
		float alpha_r = sample.x;//longitudinal shift. Suggested value from marschner hair paper between -10 and -5 degrees
		float beta_r = sample.y;//longitudinal width (stdev.). Suggested value from marschner hair paper between 5 and 10 degrees		

		//We are inside the cylinder and we want to take Cylinder Path: TT (refractive transmission, p = 1)
		if (mat_flags & BSDFHelper::MATFLAG_CYLINDER_T_BOUNCE)
		{
			//refract on the second wall of the cylinder to get output ray
			local_output_ray = glm::refract(local_input_ray, glm::faceforward(normal, local_input_ray, normal), 1.f / material->m_ior);

			//rotate towards normal to account for the tilted fiber surface
			local_output_ray = glm::vec3(glm::vec4(local_output_ray, 0.f) * glm::rotate(glm::radians(-alpha_r / 2), u_axis));

			//reset material flags and set it to finished path since we leave the cylinder
			mat_flags = 0;

			//move output ray to cylinders (fibers) local space.
			glm::vec3 out_ray_cyl = Math::worldToLocal(glm::normalize(local_output_ray), u_axis, v_axis, w_axis);

			//--------------------------------------------------------------------
			// M_tt(theta_h) : Marschner marginal, longitudinal	scattering function(M)  
			//--------------------------------------------------------------------
			//calculate parameters for M
			float theta_r = BSDFHelper::atan2ThetaToMarschner(std::atan2(out_ray_cyl.x, out_ray_cyl.z));//Angle between reflected output ray and fibers u-w(normal-) plane
			float theta_h = (theta_r + theta_i) / 2.f;//theta half angle
			float theta_d = (theta_r - theta_i) / 2.f;//theta difference angle

			//M_tt(theta_h) -> gaussian function with zero-mean and our fibers material standart derivation.
			output_pdf = BSDFHelper::normal_gauss_pdf(theta_h - glm::radians(-alpha_r / 2.f), 0.0f, beta_r / 2.f);

			//--------------------------------------------------------------------
			// N_tt(phi) : Marschner conditional, azimuthal scattering function(N)
			//--------------------------------------------------------------------

			//calculate parameters for N
			float phi_r = std::atan2(std::hypot(out_ray_cyl.x, out_ray_cyl.z), out_ray_cyl.y);
			phi_r = (out_ray_cyl.z >= 0.0f) ? phi_r : -phi_r;//wrap phi to account for std::hypot() erasing minus of z(local w)-coordinate, because marschner phi range is [-pi, pi]
			float phi = phi_r - phi_i;

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
			float dphi_dh = 1 / glm::sqrt(1 - h * h) * ((-(24 * c / glm::pow3(glm::pi<float>())) * glm::pow2(gamma_i)) + (6 * c / glm::pi<float>() - 2));

			//calculate fresnel for attenuation factor
			float fresnel = BSDFHelper::dialectricFresnel(glm::cos(gamma_i), bravais_first, bravais_sec);

			//new color absorption coefficient
			//glm::vec3 sigma_a = calcSigmaFromColor(material->fetchParameterColor<MatParamType::SIGMA_A>(texcoord));
			//glm::vec3 sigma_a = glm::vec3(material->fetchParameterColor<MatParamType::SIGMA_A>(texcoord));
			glm::vec3 sigma_a = glm::vec3(material->fetchParameterColor<MatParamType::DIFFUSE>(texcoord));

			//attenuation factor, which contains color absorbtion. Marschner above Equation 8
			glm::vec3 att_color = glm::pow2(1 - fresnel) * T(sigma_a, glm::asin(h / bravais_first));

			//final term for N_tt(phi). Marschner Equation 8.
			glm::vec3 n_tt = 0.5f * att_color / glm::abs(2 * dphi_dh);

			//--------------------------------------------------------------------

			//If we have a refraction (One vector is in positiv hemisphere and the other in negativ) we return our scat value
			if (!(glm::dot(local_input_ray, hit.m_normal) * glm::dot(local_output_ray, hit.m_normal) > 0))
				return (output_pdf * n_tt / glm::pow2(cos_theta_d));
			else
				return glm::vec3(0.f);

		}//We are inside the cylinder and we want to take Cylinder Path: TRT (refractive transmission, p = 2)
		else if (mat_flags & BSDFHelper::MATFLAG_CYLINDER_TR_BOUNCE)//We are inside the cylinder and we want to take Cylinder Path: TRT (p = 2)
		{
			//refract on the first wall of the cylinder to get output ray
			local_output_ray = glm::refract(local_input_ray, glm::faceforward(normal, local_input_ray, normal), 1.f);

			//rotate towards normal to account for the tilted fiber surface
			local_output_ray = glm::vec3(glm::vec4(local_output_ray, 0.f) * glm::rotate(glm::radians(-3.f * alpha_r / 2.f), u_axis));

			//reset material flags and set it to finished path since we leave the cylinder
			mat_flags = 0;

			//move output ray to cylinders (fibers) local space.
			glm::vec3 out_ray_cyl = Math::worldToLocal(glm::normalize(local_output_ray), u_axis, v_axis, w_axis);

			//--------------------------------------------------------------------
			// M_trt(theta_h) : Marschner marginal, longitudinal	scattering function(M)  
			//--------------------------------------------------------------------

			//calculate parameters for M
			float theta_r = BSDFHelper::atan2ThetaToMarschner(std::atan2(out_ray_cyl.x, out_ray_cyl.z));//Angle between reflected output ray and fibers u-w(normal-) plane			
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
			float phi_r = std::atan2(std::hypot(out_ray_cyl.x, out_ray_cyl.z), out_ray_cyl.y);
			phi_r = (out_ray_cyl.z >= 0.0f) ? phi_r : -phi_r;//wrap phi to account for std::hypot() erasing minus of z(local w)-coordinate, because marschner phi range is [-pi, pi]
			float phi = phi_r - phi_i;
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
				float dphi_dh_h_c = glm::pow2(1 / glm::sqrt(1 - h_c * h_c)) *
					glm::pow2(((-(48 * c / glm::pow3(glm::pi<float>())) * glm::pow2(glm::asin(h_c)) + (12 * c / glm::pi<float>() - 2))));
				h = glm::clamp(glm::min(0.5f, 2 * glm::sqrt(2 * w_c / glm::abs(dphi_dh_h_c))), -1.f, 1.f);
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
			float dphi_dh = 1 / glm::sqrt(1 - h * h) * ((-(48 * c / glm::pow3(glm::pi<float>())) * glm::pow2(gamma_i)) + (12 * c / glm::pi<float>() - 2));

			//calculate fresnel part of attenuation
			float fresnel = BSDFHelper::dialectricFresnel(glm::cos(gamma_i), bravais_first, bravais_sec);
			float gamma_t = glm::asin(h / bravais_first);
			float fresnel_exit = BSDFHelper::dialectricFresnel(glm::cos(gamma_t), 1 / bravais_first, 1 / bravais_sec);

			//new color absorption coefficient
			//glm::vec3 sigma_a = calcSigmaFromColor(material->fetchParameterColor<MatParamType::SIGMA_A>(texcoord));
			//glm::vec3 sigma_a = glm::vec3(material->fetchParameterColor<MatParamType::SIGMA_A>(texcoord));
			glm::vec3 sigma_a = glm::max(glm::vec3(material->fetchParameterColor<MatParamType::DIFFUSE>(texcoord)), 0.2f);//Min Value from Marschner Paper for sigma is 0.2

			//full attenuation
			glm::vec3 att_color = glm::pow2(1 - fresnel) * fresnel_exit * glm::pow2(T(sigma_a, gamma_t));

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

			//final scattering function. Marschner Equation 9
			glm::vec3 scat = output_pdf * (n_trt) / glm::pow2(cos_theta_d);

			//Check for NaN values
			if (std::isnan(scat.x) || std::isnan(scat.y) || std::isnan(scat.z))
				scat = glm::vec3(0.f);

			//If we have a reflection(Both vectors are in positiv hemisphere) we return our scat value
			if ((glm::dot(local_input_ray, hit.m_normal) * glm::dot(local_output_ray, hit.m_normal) > 0))
				return scat;
			else
				return glm::vec3(0.f);
		}
		else//We take Cylinder Path: R (Reflection on Surface)
		{
			//reflect input ray on the surface
			local_output_ray = glm::reflect(local_input_ray, glm::faceforward(normal, local_input_ray, normal));
			//rotate towards normal to account for the tilted fiber surface
			local_output_ray = glm::vec3(glm::vec4(local_output_ray, 0.f) * glm::rotate(glm::radians(alpha_r), u_axis));
			//Reset the mat_flags and set it to finished path since we leave the cylinder
			mat_flags = BSDFHelper::MATFLAG_SPECULAR_BOUNCE;

			//move output ray to cylinders (fibers) local space.
			/*  Our v-axis(vector through the center of the cylinder(apexpoint - basepoint)) is u-axis in marschner model*/
			glm::vec3 out_ray_cyl = Math::worldToLocal(glm::normalize(local_output_ray), u_axis, v_axis, w_axis);

			//--------------------------------------------------------------------
			// M_r(theta_h) : Marschner marginal, longitudinal	scattering function(M)  			
			//--------------------------------------------------------------------

			//calculate parameters for M
			float theta_r = BSDFHelper::atan2ThetaToMarschner(std::atan2(out_ray_cyl.x, out_ray_cyl.z));//Angle between reflected output ray and fibers u-w(normal-) plane
			float theta_h = (theta_r + theta_i) / 2;//theta half angle
			float theta_d = (theta_r - theta_i) / 2;//theta difference angle

			//M_r(theta_h) -> gaussian function with zero-mean and our fibers material standart derivation. TODO: Try Logistic function as suggested in pixar paper
			output_pdf = BSDFHelper::normal_gauss_pdf(theta_h - glm::radians(alpha_r), 0.0f, beta_r);
			//output_pdf = BSDFHelper::Logistic(theta_h - glm::radians(alpha_r), beta_r);

			//--------------------------------------------------------------------
			// N_r(phi) : Marschner conditional, azimuthal scattering function(N)
			//--------------------------------------------------------------------

			//calculate parameters for N
			float phi_r = std::atan2(std::hypot(out_ray_cyl.x, out_ray_cyl.z), out_ray_cyl.y);
			phi_r = (out_ray_cyl.z >= 0.0f) ? phi_r : -phi_r;//wrap phi to account for std::hypot() erasing minus of z(local w)-coordinate, because marschner phi range is [-pi, pi]
			float phi = phi_r - phi_i;

			//root of the approximation for phi^ in marschner paper (Equation 10) with p = 0. Concrete equation is from d'Eon paper (Above Equation 9)
			float h = glm::sin(phi) * -0.5f;

			//From ruenz12 bachelorthesis Equation 44, which derived it from marschner Equation 10.
			float dphi_dh = -2.f / glm::sqrt(1 - h * h);

			//help variables for bravais calculation
			float cos_theta_d = glm::cos(theta_d);
			float x1 = glm::sqrt(glm::pow2(material->m_ior) - glm::pow2(glm::sin(theta_d)));
			//Bravais (virtual index of reflection eta_one, eta_two) calculation
			float bravais_first = x1 / cos_theta_d;//first bravais index (virtual eta)
			float bravais_sec = glm::pow2(material->m_ior) * cos_theta_d / x1;//second bravais index (virtual eta)

			//calculate attenuation factor with fresnel
			//float fresnel = BSDFHelper::dialectricFresnel(glm::cos(glm::asin(h)), bravais_first, bravais_sec);
			float fresnel = BSDFHelper::dialectricFresnel(glm::dot(local_input_ray, hit.m_normal), bravais_first, bravais_sec);

			//final term for N_r(phi). Marschner Equation 8
			float n_r = 0.5f * fresnel / glm::abs(2 * dphi_dh);

			//--------------------------------------------------------------------

			//Final value for the combined scattering function
			float scatt_func = output_pdf * n_r / glm::pow2(cos_theta_d);
			return glm::vec3(scatt_func);

			/*//If we have a reflection(Both vectors are in positiv hemisphere) we return our scat value
			if ((glm::dot(local_input_ray, hit.m_normal) * glm::dot(local_output_ray, hit.m_normal) > 0))
				return glm::vec3(scatt_func);
			else
				return glm::vec3(0.f);*/
		}
		return glm::vec3(.0f);
	}

	glm::vec3 MarschnerHairBSDF::evaluateLight(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray)
	{
		//Normally we would check if we really have an reflection or refraction, but since we dont know if we have to check for R, TT or TRT
		//Path we do this directly in the sample_bsdf function
		return glm::vec3(0);
		return glm::vec3(1.f) * glm::one_over_pi <float>();
	}

	glm::vec3 MarschnerHairBSDF::calcSigmaFromColor(KIRK::Color::RGBA color, float beta_n)
	{
		// Compute sigma_a from the rgb color of the material. Chiang paper equation 9
		glm::vec3 sigma_a = glm::pow2(glm::log(glm::vec3(color)) /
			(5.969f - 0.215f * beta_n + 2.532f * glm::pow2(beta_n) -
				10.73f * glm::pow3(beta_n) + 5.574f * glm::pow4(beta_n) +
				0.245f * glm::pow(beta_n, 5)));
		return sigma_a;
	}

	glm::vec3 MarschnerHairBSDF::T(glm::vec3 sigma_a, float gamma_t)
	{
		return glm::exp(2.f * sigma_a * (1 + glm::cos(2.f * gamma_t)));
		//Note: Changed sigma factor from -2 to 2, because we can use RGB Diffuse colors [0, 1] this way.
		//Its easier to controll than the -2 version with an arbitrary absorption factor
	}

	////////////////////////////////////////////////////////////////////////////////////
	//////////////
	//////////////	D'EON HAIR BSDF
	//////////////
	////////////////////////////////////////////////////////////////////////////////////

	glm::vec3 DEonHairBSDF::localSample(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& normal, glm::vec2& sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance)
	{
		//Cast Object to cylinder
		KIRK::Cylinder *cylinder_obj = dynamic_cast<KIRK::Cylinder*>(hit.m_object);

		//if we have no cylinder as Object we return 0
		if (cylinder_obj == NULL)
			return glm::vec3(0.f);

		//needed parameters
		Material* material = hit.m_object->getMaterial();
		glm::vec2 texcoord = hit.m_texcoord;

		//calculation of tangent vector
		glm::vec3 c1 = glm::cross(normal, glm::vec3(0.0, 0.0, 1.0));
		glm::vec3 c2 = glm::cross(normal, glm::vec3(0.0, 1.0, 0.0));
		glm::vec3 tangent = (glm::length(c1) > glm::length(c2)) ? glm::normalize(c1) : glm::normalize(c2);

		//get light ray
		glm::vec3 hit_to_light = glm::normalize(local_output_ray);

		//calculate spherical coordinates of light ray. NOTE: std::atan2 is atan(y/x) using the signs of arguments to determine the correct quadrant
		float theta_i = std::atan2(hit_to_light.x, hit_to_light.z);//Angle between light ray and fibers v-w(normal-)plane(Inclinations with respect of the fibers normalplane). Around v-axis
		float phi_i = std::atan2(std::hypot(hit_to_light.x, hit_to_light.z), hit_to_light.y);//azimuth angle. Around x(fibers u-)-axis

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

			//get ourgoing ray
			glm::vec3 out_ray_cyl = glm::normalize(local_output_ray);

			//--------------------------------------------------------------------
			// M_tt(v, theta_i, theta_r) : d'Eon marginal, longitudinal	scattering function(M)  			
			//--------------------------------------------------------------------

			//calculate parameters for M
			float theta_r = std::atan2(out_ray_cyl.x, out_ray_cyl.z);//Angle between reflected output ray and fibers u-w(normal-) plane
			float theta_d = (theta_r - theta_i) / 2;//theta difference angle

			//M_tt -> d'Eon Equation 7
			output_pdf = M_p(glm::cos(-theta_i - glm::radians(alpha_r)), glm::cos(theta_r),
				glm::sin(-theta_i - glm::radians(alpha_r)), glm::sin(theta_r), glm::radians(beta_r * beta_r));

			//--------------------------------------------------------------------
			// N_tt(phi) : Marschner conditional, azimuthal scattering function(N)
			//--------------------------------------------------------------------

			//calculate parameters for N
			float phi_r = std::atan2(std::hypot(out_ray_cyl.x, out_ray_cyl.z), out_ray_cyl.y);
			float phi = phi_r - phi_i;

			//calculation of h from marschner but the concrete equation is from d'Eon paper (Above Equation 9)
			float a = 1 / material->m_ior;
			float nenner = glm::sqrt(1 + glm::pow2(a) - 2 * a * glm::sign(phi) * glm::sin(phi / 2.f));
			float h = (glm::sign(phi) * glm::cos(phi / 2.f)) / nenner;
			float gamma_i = glm::asin(h);

			//second part of N equation (first part is attenuation factor)
			//Gaussian detector with 20 iterations. d'Eon Equation 11
			float d_tt = 0;
			for (int i = -10; i <= 10; i++)
				d_tt += BSDFHelper::normal_gauss_pdf(phi - 2.f * glm::pi<float>() * i, 0.f, beta_r / 2);

			//calculate fresnel which is part of attenuation. d'Eon Equation 14
			float fresnel = BSDFHelper::dialectricFresnel(glm::acos(glm::cos(theta_d) * glm::cos(gamma_i)), material->m_ior, 1.f);
			if (fresnel == 1) fresnel = 0.f;//so it doesnt change att_color to 0

			//new color absorption coefficient
			//glm::vec3 sigma_a = calcSigmaFromColor(material->fetchParameterColor<MatParamType::SIGMA_A>(texcoord));
			//glm::vec3 sigma_a = glm::vec3(material->fetchParameterColor<MatParamType::SIGMA_A>(texcoord));
			glm::vec3 sigma_a = glm::max(glm::vec3(material->fetchParameterColor<MatParamType::DIFFUSE>(texcoord)), 0.2f) / glm::cos(theta_r);//Min Value from Marschner Paper for sigma is 0.2

			//full attenuation factor. d'Eon Equation 13 with p = 1
			glm::vec3 att_color = glm::pow2(1 - fresnel) * KIRK::MarschnerHairBSDF::T(sigma_a, glm::asin(h / material->m_ior));

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

			//get ourgoing ray
			glm::vec3 out_ray_cyl = glm::normalize(local_output_ray);

			//--------------------------------------------------------------------
			// M_trt(v, theta_i, theta_r) : d'Eon marginal, longitudinal scattering function(M)  			
			//--------------------------------------------------------------------

			//calculate parameters for M
			float theta_r = std::atan2(out_ray_cyl.x, out_ray_cyl.z);//Angle between reflected output ray and fibers u-w(normal-) plane
			float theta_d = (theta_r - theta_i) / 2;//theta difference angle
			sample.x = theta_i; sample.y = 0;//We store the theta_i angle in it, because we need it in the MarschnerHairShader method
											 //and this trt-path is the last call of the bsdf. We dont want to override the previously stored alpha/beta values.

			//M_trt -> d'Eon Equation 7
			output_pdf = M_p(glm::cos(-theta_i - glm::radians(4 * alpha_r)), glm::cos(theta_r),
				glm::sin(-theta_i - glm::radians(4 * alpha_r)), glm::sin(theta_r), glm::radians(beta_r * beta_r));

			//--------------------------------------------------------------------
			// N_trt(phi) : Marschner conditional, azimuthal scattering function(N)
			//--------------------------------------------------------------------

			//calculate parameters for N
			float phi_r = std::atan2(std::hypot(out_ray_cyl.x, out_ray_cyl.z), out_ray_cyl.y);
			float phi = phi_r - phi_i;
			float gamma_i = glm::angle(local_input_ray, glm::normalize(normal));//angle between input ray and surface normal in radians
			float h = glm::clamp(glm::sin(gamma_i), -1.f, 1.f);

			//second part of N equation (first part is attenuation factor)
			//Gaussian detector with 20 iterations. d'Eon Equation 11
			float d_trt = 0;
			for (int i = -10; i <= 10; i++)
				d_trt += BSDFHelper::normal_gauss_pdf(phi - 2.f * glm::pi<float>() * i, 0.f, beta_r * 2);

			//calculate fresnel which is part of attenuation. d'Eon Equation 14
			float fresnel = BSDFHelper::dialectricFresnel(glm::acos(glm::cos(theta_d) * glm::cos(gamma_i)), material->m_ior, 1.f);

			//new color absorption coefficient
			//glm::vec3 sigma_a = calcSigmaFromColor(material->fetchParameterColor<MatParamType::SIGMA_A>(texcoord));
			//glm::vec3 sigma_a = glm::vec3(material->fetchParameterColor<MatParamType::SIGMA_A>(texcoord));
			glm::vec3 sigma_a = glm::max(glm::vec3(material->fetchParameterColor<MatParamType::DIFFUSE>(texcoord)), 0.2f) / glm::cos(theta_r);//Min Value from Marschner Paper for sigma is 0.2

			//full attenuation factor. d'Eon Equation 13 with p = 2
			glm::vec3 att_color = glm::pow2(1 - fresnel) * fresnel * glm::pow2(KIRK::MarschnerHairBSDF::T(sigma_a, glm::asin(h / material->m_ior)));

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

			//get ourgoing ray
			glm::vec3 out_ray_cyl = glm::normalize(local_output_ray);

			//--------------------------------------------------------------------
			// M_r(v, theta_i, theta_r) : d'Eon marginal, longitudinal	scattering function(M)  			
			//--------------------------------------------------------------------

			//calculate parameters for M
			float theta_r = std::atan2(out_ray_cyl.x, out_ray_cyl.z);//Angle between reflected output ray and fibers u-w(normal-) plane

			//M_r -> d'Eon Equation 7
			output_pdf = M_p(glm::cos(-theta_i + 2 * alpha_r), glm::cos(theta_r), glm::sin(-theta_i + 2 * alpha_r), glm::sin(theta_r), glm::radians(beta_r * beta_r));

			//--------------------------------------------------------------------
			// N_r(phi) : d'Eon conditional, azimuthal scattering function(N)
			//--------------------------------------------------------------------

			//calculate parameters for N
			float phi_r = std::atan2(std::hypot(out_ray_cyl.x, out_ray_cyl.z), out_ray_cyl.y);
			float phi = phi_r - phi_i;
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
		bool reflect = glm::dot(local_input_ray, hit.m_normal) * glm::dot(local_output_ray, hit.m_normal) > 0;
		if (reflect) { return glm::vec3(hit.m_object->getMaterial()->fetchParameterColor <MatParamType::DIFFUSE>(hit.m_texcoord)) * glm::one_over_pi <float>(); }
		return glm::vec3(0);
	}

	float DEonHairBSDF::M_p(float cos_theta_c, float cos_theta_r, float sin_theta_c, float sin_theta_r, float v)
	{
		//Calculation of M_p -> d'Eon Equation 7
		float csch = 2 * v * glm::sinh(1.f / v);//Cosecant hyperbolic function with 2v factor
		float x = cos_theta_c * cos_theta_r / v;//parameter for bessel calculation
		float y = sin_theta_c * sin_theta_r / v;//parameter for e function calculation

		//compensation for low roughness v values from d'eon 2013 paper
		float m_p = (v <= 0.1f) ? glm::exp(BSDFHelper::logBessel(x) - y - (1 / v) + 0.6931f + glm::log(1 / (2 * v))) : (glm::exp(-y) * BSDFHelper::bessel(x)) / csch;
		return m_p;
	}


	////////////////////////////////////////////////////////////////////////////////////
	//////////////
	//////////////	CHIANG HAIR BSDF
	//////////////
	////////////////////////////////////////////////////////////////////////////////////

	glm::vec3 ChiangHairBSDF::localSample(const Intersection & hit, const glm::vec3 & local_input_ray, const glm::vec3 & normal, glm::vec2 & sample, glm::vec3 & local_output_ray, float & output_pdf, int & mat_flags, bool useRadianceOverImportance)
	{
		return glm::vec3();
	}
	glm::vec3 ChiangHairBSDF::evaluateLight(const Intersection & hit, const glm::vec3 & local_input_ray, const glm::vec3 & local_output_ray)
	{
		return glm::vec3();
	}

	void ChiangHairBSDF::init(float offset_h, float mat_beta_m, float mat_beta_n, float mat_alpha, float mat_eta, glm::vec3 color)
	{
		h = offset_h;
		gammaO = glm::asin(glm::clamp(h, -1.f, 1.f));

		//parameters from material
		beta_m = mat_beta_m;
		beta_n = mat_beta_n;
		alpha = mat_alpha;
		eta = mat_eta;

		// Compute longitudinal variance v. Chiang paper equation 7
		v[0] = glm::pow2(0.726f * beta_m + 0.812f * glm::pow2(beta_m) + 3.7f * glm::pow(beta_m, 20));
		v[1] = 0.25f * v[0];
		v[2] = 4 * v[0];
		for (int p = 3; p <= pMax; ++p)
			v[p] = v[2];

		// Compute azimuthal logistic scale factor. Chiang paper equation 8
		s = BSDFHelper::SQRT_PI_OVER_8 *
			(0.265f * beta_n + 1.194f * glm::pow2(beta_n) + 5.372f * glm::pow(beta_n, 22));

		// Compute hair scale alpha terms. 
		sin2kAlpha[0] = glm::sin(glm::radians(alpha));
		cos2kAlpha[0] = glm::sqrt(std::max(0.f, 1 - glm::pow2(sin2kAlpha[0])));
		for (int i = 1; i < 3; ++i) {
			sin2kAlpha[i] = 2 * cos2kAlpha[i - 1] * sin2kAlpha[i - 1];
			cos2kAlpha[i] = glm::pow2(cos2kAlpha[i - 1]) - glm::pow2(sin2kAlpha[i - 1]);
		}

		// Compute sigma_a from the rgb color of the material. Chiang paper equation 9
		sigma_a = glm::pow2(glm::log(color) /
			(5.969f - 0.215f * beta_n + 2.532f * glm::pow2(beta_n) -
				10.73f * glm::pow3(beta_n) + 5.574f * glm::pow4(beta_n) +
				0.245f * glm::pow(beta_n, 5)));
	}

	glm::vec3 ChiangHairBSDF::f(const glm::vec3 & wi, const glm::vec3 & wo) const
	{
		// Needed variables
		float Pi = glm::pi<float>();

		// Compute hair coordinate system terms related to _wi_
		float sinThetaO = wi.x;
		float cosThetaO = glm::sqrt(std::max(0.f, 1 - glm::pow2(sinThetaO)));
		float phiO = std::atan2(wi.z, wi.y);

		// Compute hair coordinate system terms related to _wo_
		float sinThetaI = wo.x;
		float cosThetaI = glm::sqrt(std::max(0.f, 1 - glm::pow2(sinThetaI)));
		float phiI = std::atan2(wo.z, wo.y);

		// Compute $\cos \thetat$ for refracted ray
		float sinThetaT = sinThetaO / eta;
		float cosThetaT = glm::sqrt(std::max(0.f, 1 - glm::pow2(sinThetaT)));

		// Compute $\gammat$ for refracted ray
		float etap = std::sqrt(eta * eta - glm::pow2(sinThetaO)) / cosThetaO;
		float sinGammaT = h / etap;
		float cosGammaT = glm::sqrt(std::max(0.f, 1 - glm::pow2(sinGammaT)));
		float gammaT = glm::sqrt(std::max(0.f, sinGammaT));

		// Compute the transmittance _T_ of a single path through the cylinder
		glm::vec3 T = glm::exp(-sigma_a * (2 * cosGammaT / cosThetaT));

		// Evaluate hair BSDF
		float phi = phiI - phiO;
		std::vector<glm::vec3> ap = Ap(cosThetaO, eta, h, T);
		glm::vec3 fsum(0.f);
		for (int p = 0; p < pMax; ++p) {
			// Compute $\sin \thetai$ and $\cos \thetai$ terms accounting for scales
			float sinThetaIp, cosThetaIp;
			if (p == 0) {
				sinThetaIp = sinThetaI * cos2kAlpha[1] + cosThetaI * sin2kAlpha[1];
				cosThetaIp = cosThetaI * cos2kAlpha[1] - sinThetaI * sin2kAlpha[1];
			}

			// Handle remainder of $p$ values for hair scale tilt
			else if (p == 1) {
				sinThetaIp = sinThetaI * cos2kAlpha[0] - cosThetaI * sin2kAlpha[0];
				cosThetaIp = cosThetaI * cos2kAlpha[0] + sinThetaI * sin2kAlpha[0];
			}
			else if (p == 2) {
				sinThetaIp = sinThetaI * cos2kAlpha[2] - cosThetaI * sin2kAlpha[2];
				cosThetaIp = cosThetaI * cos2kAlpha[2] + sinThetaI * sin2kAlpha[2];
			}
			else {
				sinThetaIp = sinThetaI;
				cosThetaIp = cosThetaI;
			}

			// Handle out-of-range $\cos \thetai$ from scale adjustment
			cosThetaIp = std::abs(cosThetaIp);
			fsum += DEonHairBSDF::M_p(cosThetaIp, cosThetaO, sinThetaIp, sinThetaO, v[p]) * ap[p] *
				Np(phi, p, s, gammaO, gammaT);
		}

		// Compute contribution of remaining terms after _pMax_
		fsum += DEonHairBSDF::M_p(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[pMax]) * ap[pMax] /
			(2.f * Pi);
		if (std::abs(wo.z) > 0) fsum /= std::abs(wo.z);
		return fsum;
	}

	glm::vec3 ChiangHairBSDF::Sample_f(const glm::vec3 & wi, glm::vec3 * wo, const glm::vec2 & u2, float * pdf) const
	{
		//needed variables
		float Pi = glm::pi<float>();

		// Compute hair coordinate system terms related to w_i (input_ray)
		float sinThetaO = wi.x;
		float cosThetaO = glm::sqrt(std::max(0.f, 1 - glm::pow2(sinThetaO)));
		float phiO = std::atan2(wi.z, wi.y);

		// Derive four random samples from _u2_
		glm::vec2 u[2] = { BSDFHelper::DemuxFloat(u2[0]), BSDFHelper::DemuxFloat(u2[1]) };

		// Determine which term $p$ to sample for hair scattering
		std::vector<float> apPdf = ComputeApPdf(cosThetaO);
		int p;
		for (p = 0; p < pMax; ++p) {
			if (u[0][0] < apPdf[p]) break;
			u[0][0] -= apPdf[p];
		}

		// Sample $M_p$ to compute $\thetai$
		u[1][0] = std::max(u[1][0], float(1e-5));
		float cosTheta =
			1 + v[p] * std::log(u[1][0] + (1 - u[1][0]) * std::exp(-2 / v[p]));
		float sinTheta = glm::sqrt(std::max(0.f, 1 - glm::pow2(cosTheta)));
		float cosPhi = std::cos(2 * Pi * u[1][1]);
		float sinThetaI = -cosTheta * sinThetaO + sinTheta * cosPhi * cosThetaO;
		float cosThetaI = glm::sqrt(std::max(0.f, 1 - glm::pow2(sinThetaI)));

		// Update sampled $\sin \thetai$ and $\cos \thetai$ to account for scales
		float sinThetaIp = sinThetaI, cosThetaIp = cosThetaI;
		if (p == 0) {
			sinThetaIp = sinThetaI * cos2kAlpha[1] - cosThetaI * sin2kAlpha[1];
			cosThetaIp = cosThetaI * cos2kAlpha[1] + sinThetaI * sin2kAlpha[1];
		}
		else if (p == 1) {
			sinThetaIp = sinThetaI * cos2kAlpha[0] + cosThetaI * sin2kAlpha[0];
			cosThetaIp = cosThetaI * cos2kAlpha[0] - sinThetaI * sin2kAlpha[0];
		}
		else if (p == 2) {
			sinThetaIp = sinThetaI * cos2kAlpha[2] + cosThetaI * sin2kAlpha[2];
			cosThetaIp = cosThetaI * cos2kAlpha[2] - sinThetaI * sin2kAlpha[2];
		}
		sinThetaI = sinThetaIp;
		cosThetaI = cosThetaIp;

		// Sample $N_p$ to compute $\Delta\phi$

		// Compute $\gammat$ for refracted ray
		float etap = std::sqrt(eta * eta - glm::pow2(sinThetaO)) / cosThetaO;
		float sinGammaT = h / etap;
		float gammaT = glm::asin(glm::clamp(sinGammaT, -1.f, 1.f));
		float dphi;
		if (p < pMax)
			dphi =
			BSDFHelper::Phi(p, gammaO, gammaT) + BSDFHelper::SampleTrimmedLogistic(u[0][1], s, -Pi, Pi);
		else
			dphi = 2 * Pi * u[0][1];

		// Compute _wo_ from sampled hair scattering angles
		float phiI = phiO + dphi;
		*wo = glm::vec3(sinThetaI, cosThetaI * std::cos(phiI),
			cosThetaI * std::sin(phiI));

		// Compute PDF for sampled hair scattering direction _wi_
		*pdf = 0;
		for (int p = 0; p < pMax; ++p) {
			// Compute $\sin \thetai$ and $\cos \thetai$ terms accounting for scales
			float sinThetaIp, cosThetaIp;
			if (p == 0) {
				sinThetaIp = sinThetaI * cos2kAlpha[1] + cosThetaI * sin2kAlpha[1];
				cosThetaIp = cosThetaI * cos2kAlpha[1] - sinThetaI * sin2kAlpha[1];
			}

			// Handle remainder of $p$ values for hair scale tilt
			else if (p == 1) {
				sinThetaIp = sinThetaI * cos2kAlpha[0] - cosThetaI * sin2kAlpha[0];
				cosThetaIp = cosThetaI * cos2kAlpha[0] + sinThetaI * sin2kAlpha[0];
			}
			else if (p == 2) {
				sinThetaIp = sinThetaI * cos2kAlpha[2] - cosThetaI * sin2kAlpha[2];
				cosThetaIp = cosThetaI * cos2kAlpha[2] + sinThetaI * sin2kAlpha[2];
			}
			else {
				sinThetaIp = sinThetaI;
				cosThetaIp = cosThetaI;
			}

			// Handle out-of-range $\cos \thetai$ from scale adjustment
			cosThetaIp = std::abs(cosThetaIp);
			float mp = DEonHairBSDF::M_p(cosThetaIp, cosThetaO, sinThetaIp, sinThetaO, v[p]);
			float np = Np(dphi, p, s, gammaO, gammaT);
			*pdf += mp * apPdf[p] * np;
		}
		float mp = DEonHairBSDF::M_p(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[pMax]);
		*pdf += mp * apPdf[pMax] * (1 / (2 * Pi));

		return f(wi, *wo);
	}

	float ChiangHairBSDF::Pdf(const glm::vec3 & wo, const glm::vec3 & wi) const
	{
		// Needed variables
		float Pi = glm::pi<float>();

		// Compute hair coordinate system terms related to _wi_
		float sinThetaO = wi.x;
		float cosThetaO = glm::sqrt(std::max(0.f, 1 - glm::pow2(sinThetaO)));
		float phiO = std::atan2(wi.z, wi.y);

		// Compute hair coordinate system terms related to _wo_
		float sinThetaI = wo.x;
		float cosThetaI = glm::sqrt(std::max(0.f, 1 - glm::pow2(sinThetaI)));
		float phiI = std::atan2(wo.z, wo.y);

		// Compute $\gammat$ for refracted ray
		float etap = std::sqrt(eta * eta - glm::pow2(sinThetaO)) / cosThetaO;
		float sinGammaT = h / etap;
		float gammaT = glm::sqrt(std::max(0.f, sinGammaT));

		// Compute PDF for $A_p$ terms
		std::vector<float> apPdf = ChiangHairBSDF::ComputeApPdf(cosThetaO);

		// Compute PDF sum for hair scattering events
		float phi = phiI - phiO;
		float pdf = 0;
		for (int p = 0; p < pMax; ++p) {
			// Compute $\sin \thetai$ and $\cos \thetai$ terms accounting for scales
			float sinThetaIp, cosThetaIp;
			if (p == 0) {
				sinThetaIp = sinThetaI * cos2kAlpha[1] + cosThetaI * sin2kAlpha[1];
				cosThetaIp = cosThetaI * cos2kAlpha[1] - sinThetaI * sin2kAlpha[1];
			}

			// Handle remainder of $p$ values for hair scale tilt
			else if (p == 1) {
				sinThetaIp = sinThetaI * cos2kAlpha[0] - cosThetaI * sin2kAlpha[0];
				cosThetaIp = cosThetaI * cos2kAlpha[0] + sinThetaI * sin2kAlpha[0];
			}
			else if (p == 2) {
				sinThetaIp = sinThetaI * cos2kAlpha[2] - cosThetaI * sin2kAlpha[2];
				cosThetaIp = cosThetaI * cos2kAlpha[2] + sinThetaI * sin2kAlpha[2];
			}
			else {
				sinThetaIp = sinThetaI;
				cosThetaIp = cosThetaI;
			}

			// Handle out-of-range $\cos \thetai$ from scale adjustment
			cosThetaIp = std::abs(cosThetaIp);
			pdf += DEonHairBSDF::M_p(cosThetaIp, cosThetaO, sinThetaIp, sinThetaO, v[p]) *
				apPdf[p] * Np(phi, p, s, gammaO, gammaT);
		}
		pdf += DEonHairBSDF::M_p(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[pMax]) *
			apPdf[pMax] * (1 / (2 * Pi));
		return pdf;
	}

	std::vector<float> ChiangHairBSDF::ComputeApPdf(float cosThetaO) const
	{
		// Compute array of $A_p$ values for _cosThetaO_
		float sinThetaO = glm::sqrt(std::max(0.f, 1 - cosThetaO * cosThetaO));

		// Compute $\cos \thetat$ for refracted ray
		float sinThetaT = sinThetaO / eta;
		float cosThetaT = glm::sqrt(std::max(0.f, 1 - glm::pow2(sinThetaT)));

		// Compute $\gammat$ for refracted ray
		float etap = glm::sqrt(eta * eta - glm::pow2(sinThetaO)) / cosThetaO;
		float sinGammaT = h / etap;
		float cosGammaT = glm::sqrt(glm::max(0.f, 1 - glm::pow2(sinGammaT)));

		// Compute the transmittance _T_ of a single path through the cylinder
		glm::vec3 T = glm::exp(-sigma_a * (2 * cosGammaT / cosThetaT));
		std::vector<glm::vec3> ap = Ap(cosThetaO, eta, h, T);

		// Compute $A_p$ PDF from individual $A_p$ terms
		std::vector<float> apPdf(pMax + 1); //vector with pMax + 1 elements with value 0
		float sumY =
			std::accumulate(ap.begin(), ap.end(), float(0),
				[](float s, const glm::vec3 &ap) { return s + ap.y; });
		for (int i = 0; i <= pMax; ++i) apPdf[i] = ap[i].y / sumY;
		return apPdf;
	}

	std::vector<glm::vec3> ChiangHairBSDF::Ap(float cosThetaO, float eta, float h, const glm::vec3 &T)
	{
		std::vector<glm::vec3> ap(pMax + 1);
		// Compute $p=0$ attenuation at initial cylinder intersection
		float cosGammaO = glm::sqrt(std::max(0.f, 1 - h * h));
		float cosTheta = cosThetaO * cosGammaO;
		float f = BSDFHelper::dialectricFresnel(cosTheta, 1.f, eta);
		ap[0] = glm::vec3(f);

		// Compute $p=1$ attenuation term
		ap[1] = glm::pow2(1 - f) * T;

		// Compute attenuation terms up to $p=_pMax_$
		for (int p = 2; p < pMax; ++p) ap[p] = ap[p - 1] * T * f;

		// Compute attenuation term accounting for remaining orders of scattering
		ap[pMax] = ap[pMax - 1] * f * T / (glm::vec3(1.f) - T * f);

		return ap;
	}

	float ChiangHairBSDF::Np(float phi, int p, float s, float gammaO, float gammaT) {
		float dphi = phi - BSDFHelper::Phi(p, gammaO, gammaT);
		float Pi = glm::pi<float>();
		// Remap _dphi_ to $[-\pi,\pi]$
		while (dphi > Pi) dphi -= 2 * Pi;
		while (dphi < -Pi) dphi += 2 * Pi;
		return BSDFHelper::TrimmedLogistic(dphi, s, -Pi, Pi);
	}

}
