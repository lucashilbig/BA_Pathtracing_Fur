#ifndef __CVK_RT_RAY_H
#define __CVK_RT_RAY_H

#include <glm/ext.hpp>
#include <random>

namespace KIRK {

	constexpr float cRayEpsilon = 1e-4f;

	/**
	 * @class Ray
	 * @brief A class that represents a sngle Ray, its origin and direction.
	 *        since this is more a collection of data then a class attributes are declared public
	 *        for direct access.
	 */
	class Ray
	{
	public:
		/**
		 * Default Constructor.
		 * Initalises origin to (0,0,0) and direction to (0,0,1)
		 */
		Ray();

		/**
		 * Construct ray from origin and direction as vectors.
		 * @param origin The Origin of the Ray (3d vector).
		 * @param dir The Direction of the Ray (3d vector).
		 */
		Ray(const glm::vec3 origin, const glm::vec3 dir);

		/**
		 * Construct ray from origin and direction as vectors.
		 * @param origin The Origin of the Ray (3d vector).
		 * @param dir The Direction of the Ray (3d vector).
		 */
		Ray(glm::vec3 &&origin, glm::vec3 &&dir);

		/**
		 * Follow the ray the distance lambda and get the position.
		 * @param distance The distance to follow the ray (float).
		 * @return A 3d vector to the point you get when following the direction of the ray from the origin the distance of lambda.
		 */
		glm::vec3 followDistance(const float distance) const;

		/**
		Jitters a ray by moving it's distance-normalized target position by a random amount perpendicular to it's direction. 
		@param strength Amount that the ray target will be moved.
		@param sharpness The random amount will have an exponent determined by this value.
		*/
		void jitterBy(float strength, int sharpness = 1);

		// public attributes
		glm::vec3 m_origin; ///< the origin of the ray
		glm::vec3 m_direction; ///< the direction of the ray
	};

};

#endif /* __CVK_RT_RAY_H */

