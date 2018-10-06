#ifndef __CVK_PT_BSDF_H
#define __CVK_PT_BSDF_H

#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>

#include "KIRK/Utils/Log.h"
#include "KIRK/Utils/Threading.h"
#include "KIRK/Utils/Math.h"
#include "BsdfFactory.h"
#include "KIRK/Common/Intersection.h"

namespace KIRK {

class BSDFHelper
{
public:
	static const int MATFLAG_TRANSPARENT_BOUNCE = 1 << 0;
	static const int MATFLAG_SPECULAR_BOUNCE = 1 << 1;
	static const int MATFLAG_EMISSIVE_BOUNCE = 1 << 2;
	static const int MATFLAG_CYLINDER_T_BOUNCE = 1 << 3;
	static const int MATFLAG_CYLINDER_TR_BOUNCE = 1 << 4;


	/**
	* @brief solve cubic equation x^3 + a*x^2 + b*x + c = 0
	* @param *x array of size 3
	* @return In case 3 real roots: => x[0], x[1], x[2], return 3
	*         2 real roots: x[0], x[1],          return 2
	*         1 real root : x[0], x[1] ± i*x[2], return 1
	*
	// Credit for this Method goes to: poly34.h : solution of cubic and quartic equation
	// (c) Khashin S.I. http://math.ivanovo.ac.ru/dalgebra/Khashin/index.html
	// khash2 (at) gmail.com
	*/
	static int SolveP3(double *x, double a, double b, double c);

	/**
	* @brief Calculate normalized gaussian probability density function
	* @param x The functions x value
	* @param mean The functions mean value mu
	* @param stdder The functions standart deviation sigma
	* @return float between 0 and 1
	*/
	static float normal_gauss_pdf(float x, float mean, float stddev);


	/**
	* @brief Calculate the Fresnel value with Schlick's formula which is much easier to understand and should be short and therefore pretty fast.
	* @param view The negative ray direction
	* @param normal The surface normal at the hit point
	* @param ior_in The ior of the material the ray comes from. (air = 1.0f)
	* @param ior_out The ior of the material the ray goes to. (air = 1.0f)
	*/
	static float schlickFresnel (const glm::vec3& view, const glm::vec3& normal, float ior_in, float ior_out);

	/**
	* @brief Samples a unit disk using "concentric" mapping from unit square
	* @param randoms random numbers in range [0;1)
	* @return a uniformly sampled point on the unit disk
	*/
	static glm::vec2 concentricSampleDisk (const glm::vec2& randoms);

	static glm::vec3 sampleAngle(const glm::vec2& randoms, float max_angle);

	/**
	* @brief Cosine-weighted hemisphere sampling, centered around the z-axis
	* @param randoms random numbers in range [0;1)
	* @return a cosine-weighted random point in the hemisphere around (0, 0, 1)
	*
	* This is used for lambertian reflection. We assume that the outgoing radiant flux per angle is strongest for rays closer
	* to the surface normal. Hence there should be less rays sent out in directions where the angle is small.
	*/
	static glm::vec3 cosineSampleHemisphere (const glm::vec2& randoms);
	static glm::vec3 uniformSphereSample(float u, float v);
	static float dialectricFresnel (const float cos_theta, float eta_i, float eta_t);

private:
	//=============================================================================
	// _root3, root3 from http://prografix.narod.ru
	//=============================================================================
	static double _root3(double x);
	static double root3(double x);
};

typedef std::function <glm::vec3 (const Intersection&, const glm::vec3&, const glm::vec3&, glm::vec2&, glm::vec3&, float&, int&, bool)> localSample_func;
typedef std::function <glm::vec3 (const Intersection&, const glm::vec3&, const glm::vec3&)> evaluateLight_func;

class BSDF
{
public:
    BSDF(const std::string &name, localSample_func ls, evaluateLight_func evLight) : m_name(name), m_evaluateLight(evLight), m_localSample(ls)
    {}

	/**
	* @brief Sample a bsdf
	* @param hit The intersection which caused the method call.
	* @param ray_in negative incident ray direction in world coordinates
	* @param normal The normal at the impact point.
	* @param ray_out a new ray, sampled according to the material's bsdf (OUTPUT)
	* @param pdf probability distribution function of this bsdf (OUTPUT)
	* @param u random numbers in range [0;1]
	* @return the amount of reflected/transmitted energy
	*/
	glm::vec3 sample (const Intersection& hit, const glm::vec3& ray_in, const glm::vec3& normal, glm::vec3& ray_out, float& pdf, int& mat_flags, glm::vec2& u, bool useRadianceOverImportance = true) const;

    /** evaluateLight
	* Evaluates the light influence on the hit point.
	* @param hit The intersection which caused the method call.
	* @param sample random numbers in range [0;1]
	* @return The amount of influence a light has on the hit point color.
	*/
    glm::vec3 evaluateLight (const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray) { return m_evaluateLight(hit, local_input_ray, local_output_ray); }

    std::string getName()
	{ return m_name; }

private:
	std::string m_name;
    evaluateLight_func m_evaluateLight;
    localSample_func m_localSample;
};

////////////////////////////////////////////////////////////////////////////////////
class LambertianReflectionBSDF
{
public:
	static glm::vec3 localSample (const Intersection& hit, const glm::vec3& local_space_ray, const glm::vec3& normal, const glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance = true);

	static glm::vec3 evaluateLight (const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray);
};

const BsdfRegistrator <LambertianReflectionBSDF> diffuseBSDFRegistrator("LambertianReflectionBSDF");

////////////////////////////////////////////////////////////////////////////////////
class SpecularReflectionBSDF
{
public:
	static glm::vec3 localSample (const Intersection& hit, const glm::vec3& local_space_ray, const glm::vec3& normal, const glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance = true);

	static glm::vec3 evaluateLight (const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray);
};

const BsdfRegistrator <SpecularReflectionBSDF> reflectiveBSDFRegistrator("SpecularReflectionBSDF");

class SpecularTransmissionBSDF
{
public:
	static glm::vec3 localSample(const Intersection& hit, const glm::vec3& local_space_ray, const glm::vec3& normal, const glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance = true);

	static glm::vec3 evaluateLight(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray);
};

const BsdfRegistrator <SpecularTransmissionBSDF> refractiveBSDFRegistrator("SpecularTransmissionBSDF");

class GlossyBSDF
{
public:
	static glm::vec3 localSample(const Intersection& hit, const glm::vec3& local_space_ray, const glm::vec3& normal, const glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance = true);

	static glm::vec3 evaluateLight(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray);
};

const BsdfRegistrator <GlossyBSDF> glossyBSDFRegistrator("GlossyBSDF");

////////////////////////////////////////////////////////////////////////////////////
class GlassBSDF
{
public:
	static glm::vec3 localSample (const Intersection& hit, const glm::vec3& local_space_ray, const glm::vec3& normal, const glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance = true);

	static glm::vec3 evaluateLight (const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray);
};

const BsdfRegistrator <GlassBSDF> glassBSDFRegistrator("GlassBSDF");

////////////////////////////////////////////////////////////////////////////////////
class MilkGlassBSDF
{
public:
	static glm::vec3 localSample(const Intersection& hit, const glm::vec3& local_space_ray, const glm::vec3& normal, const glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance = true);

	static glm::vec3 evaluateLight(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray);
};

const BsdfRegistrator <MilkGlassBSDF> milkGlassBSDFRegistrator("MilkGlassBSDF");

////////////////////////////////////////////////////////////////////////////////////
class LambertianTransmissionBSDF
{
public:
	static glm::vec3 localSample (const Intersection& hit, const glm::vec3& local_space_ray, const glm::vec3& normal, const glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance = true);

	static glm::vec3 evaluateLight (const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray);
};

const BsdfRegistrator <LambertianTransmissionBSDF> translucentBSDFRegistrator("LambertianTransmissionBSDF");

////////////////////////////////////////////////////////////////////////////////////
class EmissionBSDF
{
public:
	static glm::vec3 localSample (const Intersection& hit, const glm::vec3& local_space_ray, const glm::vec3& normal, const glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance = true);

	static glm::vec3 evaluateLight (const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray);
};

const BsdfRegistrator <EmissionBSDF> emissionBSDFRegistrator("EmissionBSDF");

////////////////////////////////////////////////////////////////////////////////////
class TransparentBSDF
{
public:
	static glm::vec3 localSample (const Intersection& hit, const glm::vec3& local_space_ray, const glm::vec3& normal, const glm::vec2 sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance = true);

	static glm::vec3 evaluateLight (const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray);
};

const BsdfRegistrator <TransparentBSDF> transparentBSDFRegistrator("TransparentBSDF");

////////////////////////////////////////////////////////////////////////////////////
class MarschnerHairBSDF
{
public:
	static glm::vec3 localSample(const Intersection& hit, const glm::vec3& local_space_ray, const glm::vec3& normal, glm::vec2& sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance = true);

	static glm::vec3 evaluateLight(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray);
};

const BsdfRegistrator <MarschnerHairBSDF> marschnerHairBSDFRegistrator("MarschnerHairBSDF");

////////////////////////////////////////////////////////////////////////////////////
class DEonHairBSDF
{
public:
	static glm::vec3 localSample(const Intersection& hit, const glm::vec3& local_space_ray, const glm::vec3& normal, glm::vec2& sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance = true);

	static glm::vec3 evaluateLight(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray);
private:
	/*@brief Calculates d'Eon marginal, longitudinal scattering function(M_p) after d'Eon Equation 7
	* @param cos_theta_c theta cone angle from d'Eon 2013 paper
	*/
	static float M_p(float cos_theta_c, float cos_theta_r, float sin_theta_c, float sin_theta_r, float v);

	/*@brief Calculates the bessel function (first kind, 0 order).
	* Function is taken from PBRT: Copyright(c) 1998-2016 Matt Pharr, Greg Humphreys, and Wenzel Jakob.
	*/
	inline static float bessel(float x);

	/*@brief Calculates the logarithm of the bessel function (first kind, 0 order). For x values > 12 it will use different calc than bessel.
	* Function is taken from PBRT: Copyright(c) 1998-2016 Matt Pharr, Greg Humphreys, and Wenzel Jakob.
	*/
	inline static float logBessel(float x);
};

const BsdfRegistrator <DEonHairBSDF> dEonHairBSDFRegistrator("DEonHairBSDF");

}
#endif
