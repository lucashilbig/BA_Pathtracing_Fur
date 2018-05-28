#ifndef KIRK_LIGHT_H
#define KIRK_LIGHT_H

#include <random>

#include "SceneNode.h"
#include "Color.h"
#include "Ray.h"

#include "Shading/ShaderFactory.h"
#include "Shading/Shader.h"
#include "KIRK/Utils/Gui/Gui.h"

namespace KIRK {

////////////////////////////////////////////////////////////////////////////
///////////////////
///////////////////		ABSTRACT BASE LIGHT CLASS
///////////////////
////////////////////////////////////////////////////////////////////////////

class Light : public NodeDataObject, public GuiNamedElement
{
public:
	static int next_id;
	/**
	 * \brief Abstract Light-Class Constructor
	 * \param position 
	 * \param color 
	 * \param direction 
	 * \param radius 
	 * \param att_lin Linear distance-attenuation factor
	 * \param att_quad  Quadratic distance-attenuation factor
	 * \param name The name of the Light
	 */
	Light(const glm::vec3 position, const Color::RGBA color, const glm::vec3 direction, const float radius, const float att_const, const float att_lin, const float att_quad, const std::string name) :
		GuiNamedElement(name + std::to_string(next_id)), m_color{color}, m_position{position}, m_direction{glm::normalize(direction)}, m_gui_direction(direction), m_radius{radius}, m_lin_att(att_lin), m_quad_att(att_quad), m_const(att_const), m_gui_id(next_id++), m_lshader(KIRK::ShaderFactory::getInstance().getShader("LightShader"))
	{
		std::random_device rd;
		m_gen = std::mt19937(rd());
	}

	virtual ~Light() {}

	/**
	* Calculates the direction towards the light from a given sample Position.
	* If randomize is true, a randomized direction is returned, dependant on the lights size and type.
	* The attenuation will be stored at the reference given as second argument.
	* @param samplePosition the position of the sample
	* @param attenuation The function will calculate the attenuation and save it to this reference
	* @param randomize whether to randomize the direction, defaults to false
	* @return the (randomized) Ray towards the light from samplePosition
	*/
	virtual KIRK::Ray calcLightdir(const glm::vec3& samplePosition, float& attenuation, const bool randomize = false) const = 0;

	/**
	* Calculates a random Photon ray (e.g. for photonmapper and bidirectional pathtracer)
	* @return the random Photon ray
	*/
	virtual KIRK::Ray calcRandomPhotonRay() const = 0;

	virtual KIRK::Ray calcLightBounce(glm::vec3& importance, float& att_lin, float& att_quad) const = 0;

	virtual void transform(glm::mat4 transform)
	{
		m_position = glm::vec3(transform * glm::vec4(m_position, 1.0f));
		m_direction = glm::normalize(glm::mat3(glm::transpose(glm::inverse(transform))) * m_direction);
	};

	float distanceAttenuation(float distance) const
	{
		return m_const > 0.0f || (m_lin_att > 0.0f && m_quad_att > 0.0f) ? 1.f / (m_const + m_lin_att * distance + m_quad_att * glm::pow2(distance)) : 1.0f;
	}

	/**
	 * \brief This is for BDPT and thus may or may not be correct.
	 * \param dir 
	 * \return 
	 */
	virtual float angularAttenuation(const glm::vec3 dir) const = 0;

	virtual bool isIntersection(const Ray ray, float& t) const = 0;

	virtual glm::vec3 sampleLightSource(const glm::vec3 dir, const glm::vec3 origin, const float t) const = 0;

	static bool intersectTriangle(const Ray ray, const glm::vec3 v1, const glm::vec3 v2, const glm::vec3 v3, float& t, float& u, float& v);

	glm::vec3 uniformSampleSphere() const;
	glm::vec3 cosineHemisphereSample() const;
	glm::vec3 sampleAngle(float max_angle) const;
	glm::vec3 sampleDisk(glm::vec3 normal, float radius) const;
	static void orthonormalBase(const glm::vec3& normal, glm::vec3& s, glm::vec3& t);

	Color::RGBA m_color;

	glm::vec3 m_position;
	glm::vec3 m_direction;
	glm::vec3 m_gui_direction;

	float m_radius;
	float m_lin_att;
	float m_quad_att;
	float m_const;
	float m_gui_id;

	mutable std::mt19937 m_gen;
	mutable std::uniform_real_distribution<> m_dist{ 0, 1 };

	std::shared_ptr<KIRK::Shader> m_lshader;// = std::make_shared<LightShader>("LightShader");
protected:
	void checkDirection();
};

////////////////////////////////////////////////////////////////////////////
///////////////////
///////////////////		POINT LIGHT
///////////////////
////////////////////////////////////////////////////////////////////////////

class PointLight : public Light
{
public:
	/**
	 * \brief Constructor for a PointLight.
	 * \param position 
	 * \param color 
	 * \param radius 
	 * \param att_lin 
	 * \param att_quad 
	 */
	explicit PointLight(const glm::vec3 position = {0.0f, 0.0f, 0.0f}, const Color::RGBA color = Color::WHITE, const float radius = 0.8f, const float att_const = 1.0f, const float att_lin = 0.0f, const float att_quad = 0.001f)
		: Light(position, color, {0.0f, 0.0f, 0.0f}, radius, att_const, att_lin, att_quad, "PointLight") {}

	~PointLight() {}

	KIRK::Ray calcLightdir(const glm::vec3& samplePosition, float& attenuation, const bool randomize = false) const override;

	KIRK::Ray calcRandomPhotonRay() const override;
	KIRK::Ray calcLightBounce(glm::vec3& importance, float& att_lin, float& att_quad) const override;

	bool isIntersection(const Ray ray, float& t) const override;
	float angularAttenuation(const glm::vec3 dir) const override;
	glm::vec3 sampleLightSource(const glm::vec3 dir, const glm::vec3 origin, const float t) const override;

	void onGui() override;
};

////////////////////////////////////////////////////////////////////////////
///////////////////
///////////////////		QUAD LIGHT
///////////////////
////////////////////////////////////////////////////////////////////////////

class QuadLight : public Light
{
public:

	/**
	 * \brief QuadLight Constructor.
	 * \param position 
	 * \param color 
	 * \param direction
	 * \param size 
	 * \param att_lin 
	 * \param att_quad 
	 */
	explicit QuadLight(const glm::vec3 position = {0.0f, 0.0f, 0.0f}, const Color::RGBA color = Color::WHITE, const glm::vec3 direction = {1.0f, 0.0f, 0.0f}, const glm::vec2 size = {1.0f, 1.0f}, const float att_const = 0.0f, const float att_lin = 0.0f, const float att_quad = 0.001f)
		: Light(position, color, direction, 0.0f, att_const, att_lin, att_quad, "QuadLight"), m_size(size)
	{
		calcParams();
	}

	~QuadLight() {}

	Ray calcLightdir(const glm::vec3& samplePosition, float& attenuation, const bool randomize = false) const override;

	Ray calcRandomPhotonRay() const override;
	KIRK::Ray calcLightBounce(glm::vec3& importance, float& att_lin, float& att_quad) const override;

	void transform(glm::mat4 transform) override;

	float angularAttenuation(const glm::vec3 dir) const override;

	const std::vector<glm::vec3>& getVertices() const { return m_vertices; }

	bool isIntersection(const Ray ray, float& t) const override;
	glm::vec3 sampleLightSource(const glm::vec3 dir, const glm::vec3 origin, const float t) const override;

	void onGui() override;

	glm::vec2 m_size;

private:
	void calcParams();

	std::vector<glm::vec3> m_vertices; //!< Not needed for JSON, but needed in GPU and internally in the QuadLight
};

////////////////////////////////////////////////////////////////////////////
///////////////////
///////////////////		SPOT LIGHT
///////////////////
////////////////////////////////////////////////////////////////////////////

class SpotLight : public Light
{
public:
	/**
	 * \brief SpotLight constructor.
	 * \param position 
	 * \param color 
	 * \param direction 
	 * \param outer_angle 
	 * \param inner_angle 
	 * \param att_lin 
	 * \param att_quad 
	 */
	explicit SpotLight(const glm::vec3 position = {0.0f, 0.0f, 0.0f}, const Color::RGBA color = Color::WHITE, const glm::vec3 direction = {1.0f, 0.0f, 0.0f}, const float radius = 0.5f, const float outer_angle = 22.5f, const float inner_angle = -1.f, const float att_const = 0.0f, const float att_lin = 0.0f, const float att_quad = 0.0f)
		: Light(position, color, direction, radius, att_const, att_lin, att_quad, "SpotLight"), m_inner_angle(inner_angle < 0.0f || inner_angle > outer_angle ? outer_angle : inner_angle), m_outer_angle(outer_angle) {}

	~SpotLight() {}

	KIRK::Ray calcLightdir(const glm::vec3& samplePosition, float& attenuation, const bool randomize = false) const override;

	KIRK::Ray calcRandomPhotonRay() const override;
	KIRK::Ray calcLightBounce(glm::vec3& importance, float& att_lin, float& att_quad) const override;

	bool isIntersection(const Ray ray, float& t) const override;
	float angularAttenuation(const glm::vec3 dir) const override;
	glm::vec3 sampleLightSource(const glm::vec3 dir, const glm::vec3 origin, const float t) const override;

	void onGui() override;

	float m_inner_angle;
	float m_outer_angle;
};

////////////////////////////////////////////////////////////////////////////
///////////////////
///////////////////		SUN LIGHT
///////////////////
////////////////////////////////////////////////////////////////////////////

class SunLight : public Light
{
public:
	static constexpr float infinity = 1e16;

	SunLight(const Color::RGBA color = Color::WHITE, const glm::vec3 direction = {10.f, -3.f, 3.f}, const float radius = 0.01f) :
		Light({0.0f, 0.0f, 0.0f}, color, direction, radius, 0.0f, 0.0f, 0.0f, "SunLight") {}

	~SunLight() {}

	/**
	* Calculates the direction towards the light from a given sample Position.
	* If randomize is true, a randomized direction is returned, dependant on the lights size and type.
	* @param samplePosition the position of the sample
	* @param randomize whether to randomize the direction
	* @return the (randomized) direction towards the light from samplePosition
	*/
	KIRK::Ray calcLightdir(const glm::vec3& samplePosition, float& attenuation, const bool randomize = false) const override;

	/**
	* Calculates a random Photon ray (e.g. for photonmapper and bidirectional pathtracer)
	* @return the random Photon ray
	*/
	KIRK::Ray calcRandomPhotonRay() const override;
	KIRK::Ray calcLightBounce(glm::vec3& importance, float& att_lin, float& att_quad) const override;

	bool isIntersection(const Ray ray, float& t) const override;
	float angularAttenuation(const glm::vec3 dir) const override;
	glm::vec3 sampleLightSource(const glm::vec3 dir, const glm::vec3 origin, const float t) const override;

	void onGui() override;
};
}

#endif //KIRK_LIGHT_H
