#include "Light.h"
#include "Shading/Bsdf.h"
#include "Shading/LightShader.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <numeric>

namespace KIRK {

int Light::next_id = 0;


bool Light::intersectTriangle(const Ray ray, const glm::vec3 v1, const glm::vec3 v2, const glm::vec3 v3, float& t, float& u, float& v)
{
	float EPSILON = FLT_EPSILON;

	glm::vec3 e1, e2; //Edge1, Edge2
	glm::vec3 P, Q, T;

	float det, inv_det;

	//Find vectors for two edges sharing V1
	e1 = v2 - v1;
	e2 = v3 - v1;

	//Begin calculating determinant - also used to calculate u parameter
	P = glm::cross(ray.m_direction, e2);

	//if determinant is near zero, ray lies in plane of triangle
	det = glm::dot(e1, P);
	//NOT CULLING
	if(det > -EPSILON && det < EPSILON)
		return false;
	inv_det = 1.f / det;

	//calculate distance from V1 to ray origin
	T = ray.m_origin - v1;

	//Calculate u parameter and test bound
	u = glm::dot(T, P) * inv_det;

	//The intersection lies outside of the triangle
	if(u < 0.f || u > 1.f)
		return false;

	//Prepare to test v parameter
	Q = glm::cross(T, e1);

	//Calculate V parameter and test bound
	v = glm::dot(ray.m_direction, Q) * inv_det;
	//The intersection lies outside of the triangle
	if(v < 0.f || u + v > 1.f)
		return false;

	t = glm::dot(e2, Q) * inv_det;

	if(t > EPSILON)
	{ //ray intersection
		return true;
	}

	// No hit, no win
	return false;
}

glm::vec3 Light::uniformSampleSphere() const
{
	float phi = m_dist(m_gen) * 2.f * M_PI;
	float cosTheta = 2 * m_dist(m_gen) - 1;
	float sinTheta = std::sqrt(std::max(0.0f, 1.0f - cosTheta * cosTheta));
	return {sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta};
}

glm::vec3 Light::cosineHemisphereSample() const
{
	float r = std::sqrt(m_dist(m_gen));
	float theta = 2 * M_PI * m_dist(m_gen);
	float x = r * std::cos(theta);
	float y = r * std::sin(theta);

	// Project point up to the unit sphere
	float z = std::sqrt(std::max(0.f, 1 - x * x - y * y));
	return {x, y, z};
}

glm::vec3 Light::sampleAngle(float max_angle) const
{
	float phi = m_dist(m_gen) * 2.f * M_PI;
	float cosTheta = 1.0f - m_dist(m_gen) * (1.0f - std::cos(max_angle));
	float sinTheta = std::sqrt(1.f - cosTheta * cosTheta);
	return {std::cos(phi) * sinTheta, std::sin(phi) * sinTheta, cosTheta};
}

glm::vec3 Light::sampleDisk(glm::vec3 normal, float radius) const
{
	// Sample a point on a unit disk
	float r = std::sqrt(m_dist(m_gen));
	float theta = 2 * M_PI * m_dist(m_gen);
	float x = r * std::cos(theta);
	float y = r * std::sin(theta);

	// multiply radius
	glm::vec3 dir = {radius * x, radius * y, 0.0f};

	glm::vec3 s, t;
	orthonormalBase(normal, s, t);

	// Transform into local shading coordinate system
	return dir.x * s + dir.y * t;
}

void Light::orthonormalBase(const glm::vec3& normal, glm::vec3& s, glm::vec3& t)
{
	s = std::abs(normal.x) > std::abs(normal.y) ?
		glm::vec3(-normal.z, 0, normal.x) / std::sqrt(normal.x * normal.x + normal.z * normal.z) :
		glm::vec3(0, normal.z, -normal.y) / std::sqrt(normal.y * normal.y + normal.z * normal.z);
	t = glm::cross(normal, s);
}

void Light::checkDirection()
{
	if (glm::length2(m_gui_direction) == 0)
		m_gui_direction = { 1.0f, 0.0f, 0.0f };
	m_direction = glm::normalize(m_gui_direction);
}

KIRK::Ray PointLight::calcLightdir(const glm::vec3& samplePosition, float& attenuation, const bool randomize) const
{
	glm::vec3 position = m_position;
	glm::vec3 direction = glm::normalize(position - samplePosition);

	float dirDotDir = 1.0f;
	if(randomize)
	{
		glm::vec3 point1 = uniformSampleSphere();
		position += point1 * m_radius;
		dirDotDir = glm::clamp(glm::dot(point1, -direction), 0.0f, 1.0f);
	}

	float dist = glm::length(position - samplePosition);

	attenuation = dirDotDir * distanceAttenuation(dist);

	return KIRK::Ray(samplePosition, position - samplePosition);
}

KIRK::Ray PointLight::calcRandomPhotonRay() const
{
	glm::vec3 position = m_position + m_radius * uniformSampleSphere();
	glm::vec3 normal = position - m_position;

	glm::vec3 s, t;
	orthonormalBase(normal, s, t);

	glm::vec3 hemi = cosineHemisphereSample();
	glm::vec3 direction = hemi.x * s + hemi.y * t + hemi.z * normal;

	return KIRK::Ray(position, direction);
}

KIRK::Ray PointLight::calcLightBounce(glm::vec3& importance, float& att_lin, float& att_quad) const
{
	importance = glm::vec3(m_color);
	att_lin = m_lin_att;
	att_quad = m_quad_att;
	return calcRandomPhotonRay();
}

bool PointLight::isIntersection(const Ray ray, float& t) const
{
	float radius_sq = m_radius * m_radius;
	if(radius_sq == 0)
		return false;
	if(glm::dot(ray.m_direction, ray.m_origin - m_position) > 0)
		return false;

	float a = glm::dot(ray.m_direction, ray.m_direction);
	float b = glm::dot(ray.m_direction, 2 * (ray.m_origin - m_position));
	float c = glm::dot(m_position, m_position) + glm::dot(ray.m_origin, ray.m_origin) - 2 * glm::dot(ray.m_origin, m_position) - radius_sq;
	float d = b * b - 4 * a * c;

	if(d < 0)
		return false;

	d = sqrt(d);
	t = (-0.5f) * (b + d) / a;

	return true;
}

float PointLight::angularAttenuation(const glm::vec3 dir) const
{
	return 1.0f;
}

glm::vec3 PointLight::sampleLightSource(const glm::vec3 dir, const glm::vec3 origin, const float t) const
{
	return glm::min(glm::one_over_pi<float>() * glm::vec3(m_color) / (m_const > 0 ? m_const : 1), 1.0f);
}

void PointLight::onGui()
{
	ImGui::PushID(this->m_title.c_str());
	if (ImGui::CollapsingHeader(m_title.c_str())) {
		ImGui::ColorEdit4("Color", &m_color[0], false);
		ImGui::DragFloat3("Position", &m_position[0], 0.01f, -100.0f, 100.0f);
		ImGui::DragFloat("Radius", &m_radius, 0.01f, 0.0f, 100.0f);
		ImGui::DragFloat("Constant", &m_const, 0.01f, 0.0f, 100.0f, "%.2f");
		ImGui::DragFloat("Linear", &m_lin_att, 0.01f, 0.0f, 100.0f, "%.2f");
		ImGui::DragFloat("Quadratic", &m_quad_att, 0.00001f, 0.0f, 10.0f, "%.5f");
	}
	ImGui::PopID();

}

void QuadLight::transform(glm::mat4 transform)
{
	for(auto& vec : m_vertices) { vec = glm::vec3(transform * glm::vec4(vec, 1.0f)); }
	Light::transform(transform);
}

float QuadLight::angularAttenuation(const glm::vec3 dir) const
{
	return glm::dot(glm::normalize(-dir), m_direction);
}

bool QuadLight::isIntersection(const Ray ray, float& t) const
{
	const std::vector<glm::vec3>& v = m_vertices;
	float u1, u2;
	return intersectTriangle(ray, v[0], v[1], v[3], t, u1, u2) || intersectTriangle(ray, v[2], v[3], v[1], t, u1, u2);
}

glm::vec3 QuadLight::sampleLightSource(const glm::vec3 dir, const glm::vec3 origin, const float t) const
{
	float dirDotDir = glm::dot(glm::normalize(-dir), m_direction) < 0 ? 0 : 1;

	return glm::min(/*glm::one_over_pi<float>() */ dirDotDir * glm::vec3(m_color) / (m_const > 0 ? m_const : 1), 1.f);
}

void QuadLight::onGui()
{
	ImGui::PushID(this->m_title.c_str());
	if (ImGui::CollapsingHeader(m_title.c_str())) {
		bool dirty = false;
		ImGui::ColorEdit4("Color", &m_color[0], false);
		dirty |= ImGui::DragFloat3("Position", &m_position[0], 0.01f, -100.0f, 100.0f);
		dirty |= ImGui::DragFloat3("Direction", &m_gui_direction[0], 0.01f, -10.0f, 10.0f);
		ImGui::TextColored({ 0.0f, 1.0f, 0.0f, 1.0f }, "Normalized: %.3f %.3f %.3f", m_direction[0], m_direction[1], m_direction[2]);
		dirty |= ImGui::DragFloat2("Size", &m_size[0], 0.01f, 0.01f, 200.0f);
		ImGui::DragFloat("Constant", &m_const, 0.01f, 0.0f, 100.0f, "%.2f");
		ImGui::DragFloat("Linear", &m_lin_att, 0.01f, 0.0f, 100.0f, "%.2f");
		ImGui::DragFloat("Quadratic", &m_quad_att, 0.00001f, 0.0f, 10.0f, "%.5f");
		if (dirty)
		{
			checkDirection();
			calcParams();
		}
	}
	ImGui::PopID();
}

void QuadLight::calcParams()
{
	glm::vec3 s, t;
	orthonormalBase(m_direction, s, t);

	m_vertices.resize(4);

	m_vertices[0] = m_position + s * -m_size.x / 2.0f + t * -m_size.y / 2.0f;
	m_vertices[1] = m_position + s * m_size.x / 2.0f + t * -m_size.y / 2.0f;
	m_vertices[2] = m_position + s * m_size.x / 2.0f + t * m_size.y / 2.0f;
	m_vertices[3] = m_position + s * -m_size.x / 2.0f + t * m_size.y / 2.0f;

	m_radius = glm::sqrt((m_size.x * m_size.y) / glm::pi<float>());
}

KIRK::Ray QuadLight::calcLightdir(const glm::vec3& samplePosition, float& attenuation, const bool randomize) const
{
	glm::vec3 lightDir = m_position - samplePosition;
	if(randomize)
	{
		float u = m_dist(m_gen);
		float v = m_dist(m_gen);
		// Bilinear Interpolation
		glm::vec3 x1 = m_vertices[0] + u * (m_vertices[1] - m_vertices[0]);
		glm::vec3 x2 = m_vertices[3] + u * (m_vertices[2] - m_vertices[3]);
		glm::vec3 interpol = x1 + v * (x2 - x1);

		lightDir = interpol - samplePosition;
	}
	float dirDotDir = glm::clamp(glm::dot(glm::normalize(-lightDir), m_direction), 0.0f, 1.0f);
	attenuation = dirDotDir * distanceAttenuation(glm::length(lightDir));

	return KIRK::Ray(samplePosition, lightDir);
}

KIRK::Ray QuadLight::calcRandomPhotonRay() const
{
	float u = m_dist(m_gen);
	float v = m_dist(m_gen);
	// Bilinear Interpolation
	glm::vec3 x1 = m_vertices[0] + u * (m_vertices[1] - m_vertices[0]);
	glm::vec3 x2 = m_vertices[3] + u * (m_vertices[2] - m_vertices[3]);
	glm::vec3 interpol = x1 + v * (x2 - x1);

	const glm::vec3& normal = m_direction;
	glm::vec3 s, t;
	orthonormalBase(normal, s, t);

	glm::vec3 hemi = cosineHemisphereSample();
	glm::vec3 direction = hemi.x * s + hemi.y * t + hemi.z * normal;

	return KIRK::Ray(interpol, direction);
}

KIRK::Ray QuadLight::calcLightBounce(glm::vec3& importance, float& att_lin, float& att_quad) const
{
	Ray ray = calcRandomPhotonRay();
	float dirDotDir = glm::clamp(glm::dot(glm::normalize(ray.m_direction), m_direction), 0.0f, 1.0f);
	importance = glm::vec3(m_color) * dirDotDir;
	att_lin = m_lin_att;
	att_quad = m_quad_att;
	return ray;
}

KIRK::Ray SpotLight::calcLightdir(const glm::vec3& samplePosition, float& attenuation, const bool randomize) const
{
	glm::vec3 lightDir = m_position - samplePosition;
	if(randomize)
	{
		glm::vec3 point = sampleDisk(m_direction, m_radius);

		lightDir = m_position + point - samplePosition;
	}

	float lightToSurfaceAngle = glm::degrees(glm::acos(glm::dot(glm::normalize(-lightDir), m_direction)));
	float delta = 1 - glm::clamp((lightToSurfaceAngle - m_inner_angle) / ((m_outer_angle - m_inner_angle)), 0.f, 1.f);
	delta = delta * delta * delta * delta;
	attenuation = delta * distanceAttenuation(glm::length(lightDir));

	return KIRK::Ray(samplePosition, lightDir);
}

KIRK::Ray SpotLight::calcRandomPhotonRay() const
{
	glm::vec3 point = sampleDisk(m_direction, m_radius);
	glm::vec3 direction = sampleAngle(glm::radians(m_outer_angle));
	glm::vec3 s, t;
	orthonormalBase(m_direction, s, t);
	glm::vec3 worldDirection = direction.x * s + direction.y * t + direction.z * m_direction;
	return KIRK::Ray(point + m_position, worldDirection);
}

KIRK::Ray SpotLight::calcLightBounce(glm::vec3& importance, float& att_lin, float& att_quad) const
{
	Ray ray = calcRandomPhotonRay();
	float lightToSurfaceAngle = glm::degrees(acos(glm::dot(glm::normalize(ray.m_direction), m_direction)));
	float delta = 1 - glm::clamp((lightToSurfaceAngle - m_inner_angle) / ((m_outer_angle - m_inner_angle)), 0.f, 1.f);
	delta = delta * delta * delta * delta;
	importance = delta * glm::vec3(m_color);
	att_lin = m_lin_att;
	att_quad = m_quad_att;
	return ray;
}

bool SpotLight::isIntersection(const Ray ray, float& t) const
{
	if(m_radius == 0)
		return false;
	// Create a Triangle that fits in the Spotlight's unit-disk
	glm::vec3 v1 = m_position;
	glm::vec3 normal = m_direction;
	// orthonormal basis for spotlight
	glm::vec3 x = glm::abs(normal.x) > glm::abs(normal.y) ?
		              glm::vec3(-normal.z, 0, normal.x) / glm::sqrt(normal.x * normal.x + normal.z * normal.z) :
		              glm::vec3(0, normal.z, -normal.y) / glm::sqrt(normal.y * normal.y + normal.z * normal.z);
	glm::vec3 y = glm::cross(normal, x);
	glm::vec3 v2 = m_position + x;
	glm::vec3 v3 = m_position + y;

	glm::vec3 e1, e2; //Edge 1, Edge 2
	glm::vec3 P, Q, T;
	float det, inv_det;

	//Find vectors for two edges sharing V1
	e1 = v2 - v1;
	e2 = v3 - v1;

	//Begin calculating determinant - also used to calculate u parameter
	P = glm::cross(ray.m_direction, e2);

	//if determinant is near zero, ray lies in plane of triangle
	det = glm::dot(e1, P);

	//NOT CULLING
	if(det > -FLT_EPSILON && det < FLT_EPSILON)
		return false;

	inv_det = 1.f / det;

	//calculate distance from V1 to ray origin
	T = ray.m_origin - v1;

	//Calculate u parameter
	float u = glm::dot(T, P) * inv_det;

	//Prepare to test v parameter
	Q = glm::cross(T, e1);

	//Calculate V parameter
	float v = glm::dot(ray.m_direction, Q) * inv_det;

	// u and v are in coordinates of the x and y axis in the spotlights orthonormal coordinate system
	// the length of the vector (u, v) has to be smaller or equal to the spotlights radius
	if(u * u + v * v > m_radius * m_radius)
		return false;

	t = glm::dot(e2, Q) * inv_det;

	if(t > FLT_EPSILON)
	{ //ray intersection
		return true;
	}

	// No hit, no win
	return false;
}

float SpotLight::angularAttenuation(const glm::vec3 dir) const
{
	float lightToSurfaceAngle = glm::degrees(glm::acos(glm::dot(glm::normalize(-dir), m_direction)));
	return 1.0f - glm::clamp((lightToSurfaceAngle - m_inner_angle) / ((m_outer_angle - m_inner_angle)), 0.0f, 1.0f);
}

glm::vec3 SpotLight::sampleLightSource(const glm::vec3 dir, const glm::vec3 origin, const float t) const
{
	float dirDotDir = glm::dot(glm::normalize(-dir), m_direction) < 0.0f ? 0.0f : 1.0f;
	return glm::min(glm::one_over_pi<float>() * dirDotDir * glm::vec3(m_color) / (m_const > 0 ? m_const : 1), 1.0f);
}

void SpotLight::onGui()
{
	ImGui::PushID(this->m_title.c_str());
	if (ImGui::CollapsingHeader(m_title.c_str())) {
		bool dirty = false;
		ImGui::ColorEdit4("Color", &m_color[0], false);
		ImGui::DragFloat3("Position", &m_position[0], 0.01f, -100.0f, 100.0f);
		dirty |= ImGui::DragFloat3("Direction", &m_gui_direction[0], 0.01f, -10.0f, 10.0f);
		ImGui::TextColored({ 0.0f, 1.0f, 0.0f, 1.0f }, "Normalized: %.3f %.3f %.3f", m_direction[0], m_direction[1], m_direction[2]);
		ImGui::DragFloat("Radius", &m_radius, 0.01f, 0.0f, 100.0f);
		ImGui::DragFloat("Outer", &m_outer_angle, 0.1f, m_inner_angle, 180.0f);
		dirty |= ImGui::DragFloat("Inner", &m_inner_angle, 0.1f, 0.0f, m_outer_angle);
		ImGui::DragFloat("Constant", &m_const, 0.01f, 0.0f, 100.0f, "%.2f");
		ImGui::DragFloat("Linear", &m_lin_att, 0.01f, 0.0f, 100.0f, "%.2f");
		ImGui::DragFloat("Quadratic", &m_quad_att, 0.00001f, 0.0f, 10.0f, "%.5f");
		if (dirty)
			checkDirection();
	}
	ImGui::PopID();
}

KIRK::Ray SunLight::calcLightdir(const glm::vec3& samplePosition, float& attenuation, const bool randomize) const
{
	glm::vec3 point = m_radius * uniformSampleSphere();

	point -= m_direction;

	glm::vec3 direction = glm::normalize(point);

	glm::vec3 position = 1e16*direction;

	attenuation = 1;
	return KIRK::Ray(samplePosition, position - samplePosition);
}

KIRK::Ray SunLight::calcRandomPhotonRay() const
{
	glm::vec3 point = m_radius * uniformSampleSphere();

	point -= m_direction;

	glm::vec3 direction = glm::normalize(point);

	glm::vec3 position = 1e16*direction;
	return KIRK::Ray(position, m_direction);
}

KIRK::Ray SunLight::calcLightBounce(glm::vec3& importance, float& att_lin, float& att_quad) const
{
	importance = glm::vec3(m_color);
	att_lin = 0;
	att_quad = 0;
	return calcRandomPhotonRay();
}

bool SunLight::isIntersection(const Ray ray, float& t) const
{
	//We do not test for intersection with the sun.
	return false;
}

float SunLight::angularAttenuation(const glm::vec3 dir) const
{
	return 1;
}

glm::vec3 SunLight::sampleLightSource(const glm::vec3 dir, const glm::vec3 origin, const float t) const
{
	return glm::vec3(m_color);
}

void SunLight::onGui()
{
	ImGui::PushID(this->m_title.c_str());
	if (ImGui::CollapsingHeader(m_title.c_str())) {
		bool dirty = false;
		ImGui::ColorEdit4("Color", &m_color[0], false);
		dirty |= ImGui::DragFloat3("Direction", &m_gui_direction[0], 0.01f, -10.0f, 10.0f);
		ImGui::TextColored({ 0.0f, 1.0f, 0.0f, 1.0f }, "Normalized: %.3f %.3f %.3f", m_direction[0], m_direction[1], m_direction[2]);
		m_direction = glm::normalize(m_direction);
		ImGui::DragFloat("Radius", &m_radius, 0.001f, 0.0f, 1.0f);
		if (dirty)
			checkDirection();
	}
	ImGui::PopID();
}
}
