#ifndef __CVK_PT_BSDF_H
#define __CVK_PT_BSDF_H

#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <vector>

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

	static const float SQRT_PI_OVER_8;


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

	/*@brief Calculates the bessel function (first kind, 0 order).
	* Function is taken from PBRT: Copyright(c) 1998-2016 Matt Pharr, Greg Humphreys, and Wenzel Jakob.
	*/
	static float bessel(float x);

	/*@brief Calculates the logarithm of the bessel function (first kind, 0 order). For x values > 12 it will use different calc than bessel.
	* Function is taken from PBRT: Copyright(c) 1998-2016 Matt Pharr, Greg Humphreys, and Wenzel Jakob.
	*/
	static float logBessel(float x);

	/*@brief Calculates the logistic function
	* Function is taken from PBRT: Copyright(c) 1998-2016 Matt Pharr, Greg Humphreys, and Wenzel Jakob.
	*/
	static float Logistic(float x, float s);

	/*@brief 
	* Function is taken from PBRT: Copyright(c) 1998-2016 Matt Pharr, Greg Humphreys, and Wenzel Jakob.
	*/
	static float LogisticCDF(float x, float s);

	/*@brief
	* Function is taken from PBRT: Copyright(c) 1998-2016 Matt Pharr, Greg Humphreys, and Wenzel Jakob.
	*/
	static float TrimmedLogistic(float x, float s, float a, float b);

	/*@brief Samples the trimmed logistic function
	* Function is taken from PBRT: Copyright(c) 1998-2016 Matt Pharr, Greg Humphreys, and Wenzel Jakob.
	*/
	static float SampleTrimmedLogistic(float u, float s, float a, float b);

	/*@brief
	* Function is taken from PBRT: Copyright(c) 1998-2016 Matt Pharr, Greg Humphreys, and Wenzel Jakob.
	*/
	static float Phi(int p, float gammaO, float gammaT);

	/*@brief
	* Function is taken from PBRT: Copyright(c) 1998-2016 Matt Pharr, Greg Humphreys, and Wenzel Jakob.
	*/
	static glm::vec2 DemuxFloat(float f);

	// https://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/
	static uint32_t Compact1By1(uint32_t x);

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
	/* Calculates the sigma color absorption from an rgb color using Chiang paper Equation 9
	@param color. RGB Color which the hair should have
	@param beta_n. azimuthal roughness parameter. Values [0, 1]
	*/
	static glm::vec3 calcSigmaFromColor(KIRK::Color::RGBA color, float beta_n = 0.3f);

	/* Calculates the color absorption factor using Equation of Marschner Paper section 4.3
	@param sigma_a. The color absorption coefficient sigma_a
	@param gamma_t. fibers internal offset angle gamma_t
	*/
	static glm::vec3 T(glm::vec3 sigma_a, float gamma_t);
	
};

const BsdfRegistrator <MarschnerHairBSDF> marschnerHairBSDFRegistrator("MarschnerHairBSDF");

////////////////////////////////////////////////////////////////////////////////////
class DEonHairBSDF
{
public:
	static glm::vec3 localSample(const Intersection& hit, const glm::vec3& local_space_ray, const glm::vec3& normal, glm::vec2& sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance = true);

	static glm::vec3 evaluateLight(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray);

	/*@brief Calculates d'Eon marginal, longitudinal scattering function(M_p) after d'Eon Equation 7
	* @param cos_theta_c theta cone angle from d'Eon 2013 paper
	*/
	static float M_p(float cos_theta_c, float cos_theta_r, float sin_theta_c, float sin_theta_r, float v);
	
};

const BsdfRegistrator <DEonHairBSDF> dEonHairBSDFRegistrator("DEonHairBSDF");

////////////////////////////////////////////////////////////////////////////////////

class ChiangHairBSDF
{
	/* Implements "A Practical and Controllable Hair and Fur Model for Production Path Tracing" from Chiang et. al. (2016)
	Most parts of this class implementation are taken from PBRTv3 therefor Credit and Copyright(c) goes to 1998-2016 Matt Pharr, Greg Humphreys, and Wenzel Jakob.
	*/
public:
	//Sample functions
	static glm::vec3 localSample(const Intersection& hit, const glm::vec3& local_space_ray, const glm::vec3& normal, glm::vec2& sample, glm::vec3& local_output_ray, float& output_pdf, int& mat_flags, bool useRadianceOverImportance = true);

	static glm::vec3 evaluateLight(const Intersection& hit, const glm::vec3& local_input_ray, const glm::vec3& local_output_ray);

	glm::vec3 Sample_f(const glm::vec3 &wi, glm::vec3 *wo, const glm::vec2 &u2, float *pdf) const;

	//Extracts needed variables from material and sets members of bsdf. NEEDS TO BE CALLED BEFORE USING BSDF
	void init(float offset_h, float mat_beta_m, float mat_beta_n, float mat_alpha, float mat_eta, glm::vec3 color);
	
private:
	glm::vec3 f(const glm::vec3 &wi, const glm::vec3 &wo) const;
	float Pdf(const glm::vec3 &wi, const glm::vec3 &wo) const;

	std::vector<float> ComputeApPdf(float cosThetaO) const;

	static std::vector<glm::vec3> Ap(float cosThetaO, float eta, float h, const glm::vec3 &T);

	static float Np(float phi, int p, float s, float gammaO, float gammaT);
	
	//////////////////
	////// Member: Hair private data
	/////////////////
	float h;//offset along the cylinder width where the ray intersected the circular cross section. Values [-1, 1]
	float gammaO;//angle between the outgoing ray and normal

	glm::vec3 sigma_a;//Color absorbtion coefficent
	float sin2kAlpha[3], cos2kAlpha[3];//Alpha terms for hair scale

	static const int pMax = 3;//maximal internal fiber pathes that will be taken by the bsdf

	//Controllability parameters from Chiang paper section 4
	float eta;//Index of Refraction
	float v[pMax + 1];//longitudinal surface roughness variance v. Depending on beta_m
	float s;//logistic scale factor. Depending on beta_n
	float beta_m;//longitudinal roughness parameter. Values [0, 1]
	float beta_n;//azimuthal roughness parameter Values [0, 1]
	float alpha;//angle that the hair scales are tilted at the surface. Stored in degrees	
};

const BsdfRegistrator <ChiangHairBSDF> chiangHairBSDFRegistrator("ChiangHairBSDF");


}
#endif
