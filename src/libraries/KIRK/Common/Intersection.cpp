#include "Intersection.h"

namespace KIRK
{

	Intersection::Intersection(const Ray ray) : m_ray(ray)
	{
		m_lambda = FLT_MAX;
		m_triangle = 0;
	}

	void Intersection::update(Triangle *obj, float lambda, bool enter)
	{
		m_triangle = obj;
		m_lambda = lambda;
		m_enter = enter;
		m_location = m_ray.followDistance(lambda);
		m_barycentric_coord = glm::vec3(0);
	}

	void Intersection::update(Triangle *obj, float lambda, bool enter, glm::vec3 bcc)
	{
		m_triangle = obj;
		m_lambda = lambda;
		m_enter = enter;
		m_location = m_ray.followDistance(lambda);
		m_barycentric_coord = bcc;
	}


}