#include "Ray.h"

namespace KIRK {

	Ray::Ray()
	{
		m_origin = glm::vec3(0, 0, 0);
		m_direction = glm::vec3(0, 0, 1);
	}

	Ray::Ray(const glm::vec3 origin, const glm::vec3 dir)
	{
		m_origin = origin;
		m_direction = dir;
	}

	Ray::Ray(glm::vec3 &&origin, glm::vec3 &&dir)
	{
		m_origin = origin;
		m_direction = dir;
	}

	glm::vec3 Ray::followDistance(const float distance) const
	{
		return (distance * m_direction + m_origin);
	}

	void Ray::jitterBy(float strength, int sharpness) {
		//Use the maximum of those two cross products for the rare case that one of them is of zero length
		glm::vec3 perpendicular = glm::normalize(glm::max(cross(m_direction, glm::vec3(0.0, 0.0, 1.0)), cross(m_direction, glm::vec3(0.0, 1.0, 0.0))));
		glm::vec3 secondary_perpendicular = glm::normalize(glm::cross(m_direction, perpendicular));

		//We can use those static variables (as they are only assigned once, not per ray)...
		//(They should be thread safe as they are thread_local)
		static thread_local std::mt19937_64 generator;
		static thread_local std::uniform_real_distribution<float> distribution(0.f, 1.f);
		//...to calculate those two random numbers...
		float rnd_a = distribution(generator);
		float rnd_b = distribution(generator);
		//...and swap them if necessary...
		if (rnd_b < rnd_a) {
			//Swap
			float tmp = rnd_b;
			rnd_b = rnd_a;
			rnd_a = tmp;
		}
		//...for creating a thread safe disk_rand functionality
		glm::vec2 rand_disk = glm::vec2(strength*rnd_b * glm::cos(2 * glm::pi<float>()*rnd_a / rnd_b), strength*rnd_b*glm::sin(2 * glm::pi<float>()*rnd_a / rnd_b));

		int sign_perpendicular = 1;
		int sign_perpendicular_sec = 1;

		if (sharpness % 2 == 0) {
			//For an even number, the sign gets lost, so we have to remember it.vg
			sign_perpendicular = glm::sign(rand_disk.x);
			sign_perpendicular_sec = glm::sign(rand_disk.y);
		}

		m_direction += sign_perpendicular * glm::pow(rand_disk.x, sharpness)*perpendicular + sign_perpendicular_sec * glm::pow(rand_disk.y, sharpness)*secondary_perpendicular;
		m_direction = glm::normalize(m_direction);

	}


}