#include "Intersection.h"

namespace KIRK
{

	Intersection::Intersection(const Ray ray) : m_ray(ray)
	{
		m_lambda = FLT_MAX;
		m_object = 0;
	}

	void Intersection::update(Object *obj, float lambda, bool enter)
	{
		m_object = obj;
		m_lambda = lambda;
		m_enter = enter;
		m_location = m_ray.followDistance(lambda);
		m_barycentric_coord = glm::vec3(0);
	}

	void Intersection::update(Object *obj, float lambda, bool enter, glm::vec3 bcc)
	{
		m_object = obj;
		m_lambda = lambda;
		m_enter = enter;
		m_location = m_ray.followDistance(lambda);
		m_barycentric_coord = bcc;
	}


}