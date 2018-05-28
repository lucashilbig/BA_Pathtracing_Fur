#include "Math.h"

namespace KIRK
{
	namespace Math
	{
		glm::vec3 worldToLocal(const glm::vec3 &vector, const glm::vec3 &local_x, const glm::vec3 &local_y, const glm::vec3 &local_z)
		{
			return glm::vec3(glm::dot(vector, local_x), glm::dot(vector, local_y), glm::dot(vector, local_z));
		}

		glm::vec3 localToWorld(const glm::vec3 &vector, const glm::vec3 &local_x, const glm::vec3 &local_y, const glm::vec3 &local_z)
		{
			return vector.x * local_x + vector.y * local_y + vector.z * local_z;
		}
		glm::vec3 localToWorldNormal(const glm::vec3 &vector, const glm::vec3 &normal)
		{
			glm::vec3 n = normal;
			glm::vec3 dx0 = glm::vec3(0, n.z, -n.y);
			glm::vec3 dx1 = glm::vec3(-n.z, 0, n.x);
			glm::vec3 s = glm::normalize(n.y*n.y > n.x*n.x ? dx0 : dx1);
			glm::vec3 t = glm::normalize(glm::cross(n, s));
			return localToWorld(vector, s, t, n);
		}
	}
}