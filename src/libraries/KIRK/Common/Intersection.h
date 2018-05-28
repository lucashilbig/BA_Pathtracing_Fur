#ifndef __CVK_RT_INTERSECION_H
#define __CVK_RT_INTERSECION_H

#include <glm/ext.hpp>
#include <random>

#include "Ray.h"

namespace KIRK {

	class Triangle;
	
	class Intersection
	{
	public:

		/**
		* Construct Intersection with ray
		* @param ray  to be put in the Intersection
		*/
		Intersection(const Ray ray);

		/**
		* Updates some attributes of the intersection.
		* Calculates position from lambda and the ray and sets barycentric coordinates to zero.
		* @param obj pointer to the object that is intersected YOU NEED TO ENSURE THE LIFETIME OF THE OBJECT EXEEDS THE LIFETIME OF THE INTERSECTION
		* @param lambda distance from ray origin to intersection
		* @param enter probably if the ray enters the object or leaves it, but not realy sure since it is never used (WTF)
		*/
		void update(Triangle *obj, const float lambda, const bool enter);

		/**
		* Updaes some attributes of the intersection.
		* Calculates position from lambda and the ray.
		* @param obj pointer to the object that is intersected YOU NEED TO ENSURE THE LIFETIME OF THE OBJECT EXEEDS THE LIFETIME OF THE INTERSECTION
		* @param lambda distance from ray origin to intersection
		* @param enter probably if the ray enters the object or leaves it, but not realy sure since it is never used (WTF)
		* @param bcc barycentric coordinates of the intersection point e.g. for triangles
		*/
		void update(Triangle *obj, const float lambda, const bool enter, const glm::vec3 bcc);

		// public attributes
		Ray m_ray; ///< the Ray object which intersects with.
		float m_lambda; ///< distance from ray origin to intersection
		Triangle* m_triangle; ///< pointer to the object that is intersected YOU NEED TO ENSURE THE LIFETIME OF THE OBJECT EXEEDS THE LIFETIME OF THE INTERSECTION
		bool m_enter; ///< probably if the ray enters the object or leaves it, but not realy sure since it is never used (WTF)
		glm::vec3 m_location; ///< position of the intersection
		glm::vec3 m_barycentric_coord; ///< barycentric coordinates of the intersection point e.g. for triangles
		glm::vec2 m_texcoord; ///< texture coordinates of the intersection point
		glm::vec3 m_normal; ///< normal vector of the intersection point
	};

};

#include "Triangle.h"

#endif /* __CVK_RT_RAY_H */

