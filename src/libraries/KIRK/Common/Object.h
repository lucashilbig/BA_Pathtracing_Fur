#ifndef KIRK_OBJECT_H
#define KIRK_OBJECT_H

#include "KIRK/Common/Intersection.h"

namespace KIRK {
	class Material;
	class Triangle;
	class Cylinder;
	/*
	 * Baseclass for Raytracing Objects which have to implement the Intersection Methods
	 */
	class Object
	{
	public:

		virtual bool closestIntersection(Intersection *hit, float tMin, float tMax) = 0;
		virtual bool isIntersection(Ray *ray, float tMax) = 0;
		virtual void calcTcoord(Intersection *hit) = 0;
		virtual void calcNormal(Intersection *hit) = 0;
		virtual bool isInAABB(glm::vec3 *bbox) = 0;

		int getLongestAxis();

	    void setMaterial(KIRK::Material *material);
		KIRK::Material *getMaterial();

		glm::vec3 *getBounds();
		glm::vec3 getMinBound() const;
		glm::vec3 getMaxBound() const;
		glm::vec3 getCentroid() const;

	protected:

		virtual void computeBounds() = 0;

		glm::vec3 m_bound[2];//Objects Axis-Aligned-Bounding-Box
		glm::vec3 m_sizeIndicator;//sizeIndicator of the BB
		glm::vec3 m_centroid;//Objects Centert point. Used for BVH Construction

		KIRK::Material *m_material;//Objects Material
		int m_lA; // longest Axis of the BB: 0 = x-axis, 1 = y-axis, 2 = z-axis
	};
};
#include "KIRK/Common/Material.h"
#include "KIRK/Common/Triangle.h"
#include "KIRK/Common/Cylinder.h"

#endif /* KIRK_OBJECT_H */