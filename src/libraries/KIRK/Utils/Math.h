#ifndef __KIRK_MATH_H
#define __KIRK_MATH_H

#include <glm/glm.hpp>

namespace KIRK
{
	namespace Math
	{
		/**
		 * \brief Transforms a given direction vector into a local coordinate system determined by the x, y and z axises.
		 * \param vector World direction vector.
		 * \param local_x The local coordinate system's global x-axis-direction.
		 * \param local_y The local coordinate system's global y-axis-direction.
		 * \param local_z The local coordinate system's global z-axis-direction.
		 * \return The direction vector transformed into the local space.
		 */
		glm::vec3 worldToLocal(const glm::vec3 &vector, const glm::vec3 &local_x, const glm::vec3 &local_y, const glm::vec3 &local_z);

		/**
		* \brief Transforms a given direction vector from a local coordinate system determined 
		* by the x, y and z axises back to a global coordinate system.
		* \param vector Local direction vector.
		* \param local_x The local coordinate system's global x-axis-direction.
		* \param local_y The local coordinate system's global y-axis-direction.
		* \param local_z The local coordinate system's global z-axis-direction.
		* \return The direction vector transformed from the local space back to global space.
		*/
		glm::vec3 localToWorld(const glm::vec3 &vector, const glm::vec3 &local_x, const glm::vec3 &local_y, const glm::vec3 &local_z);
	/**
		 * \brief Transforms a given local vector into the coordinate system given by normal
		 * \param vector Local vector to transform.
		 * \param normal Normal that defines the coordinate system.
		 * \return The transformed vector.
		 */
		glm::vec3 localToWorldNormal(const glm::vec3 &vector, const glm::vec3 &normal);
	}
}

#endif // !__KIRK_MATH_H
