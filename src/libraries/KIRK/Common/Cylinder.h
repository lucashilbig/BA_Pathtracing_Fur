#ifndef KIRK_OBJCYLINDER_H
#define KIRK_OBJCYLINDER_H

#include "KIRK/Common/Intersection.h"

namespace KIRK {

	class Material;

	/* Class for Cylinder Intersection used for Raytracing. In this implementation the Cylinder is in a 
	*  local object space and the ray will be transformed from world space to this local object space.
	*  Because Kirk works in world space, every public variable that goes in or out of this class will
	*  internal be transformed to local space for calculations and than the return value will be transformed
	*  back to world space.
	*  This class is a modified version of the Cylinder.h class of pbrt. You can check the code here:
	*  https://www.csie.ntu.edu.tw/~cyy/courses/rendering/pbrt-2.00/html/cylinder_8h_source.html
	*  Therefor all Rights go to: pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.
	*/
	class Cylinder
	{
	public:
		// Cylinder Public Methods
		Cylinder(glm::vec3 basePoint, glm::vec3 apexPoint, float baseRadius, float apexRadius, glm::mat4 transform_objectToWolrd, glm::mat4 transform_worldToObject);
		Cylinder(glm::vec3 basePoint, glm::vec3 apexPoint, float baseRadius, float apexRadius);

		~Cylinder();

		bool closestIntersection(Intersection *hit, float tMin, float tMax);
		bool isIntersection(Ray *ray, float tMax);
		bool closestIntersectionAsPlane(Intersection *hit, float tMin, float tMax);
		void calcNormal(Intersection *hit);
		bool isInAABB(glm::vec3 *bbox);

		glm::vec3 *getBasePoint() { return &m_basepoint; }
		glm::vec3 *getApexPoint() { return &m_apexpoint; }
		glm::vec3 getVecBaseApex() { return (m_basepoint - m_apexpoint); }
		float getBaseRadius() { return m_baseradius; }
		float getApexRadius() { return m_apexradius; }
		int getLongestAxis() { return m_lA; }

		inline void setMaterial(KIRK::Material *material) {
			m_material = material;
		}
		inline KIRK::Material *getMaterial() {
			return m_material;
		}

		glm::vec3 *getBounds()  { return m_bound; }
		glm::vec3 getMinBound() const { return glm::vec3((trans_o2w * glm::vec4(m_bound[0], 1.0f))); }//transform bound from local to world space before returning
		glm::vec3 getMaxBound() const { return glm::vec3((trans_o2w * glm::vec4(m_bound[1], 1.0f))); }// "  "
		glm::vec3 getCentroid() const { return glm::vec3((trans_o2w * glm::vec4(m_centroid, 1.0f))); }// "  "

	private:
		// Cylinder Private Data

		void computeBounds();

		bool testAxis(const glm::vec3 e, const glm::vec3 v[], const glm::vec3 boxhalfsize, const int i, const int j, const int vi1, const int vi2, const int f, float &min, float &max, float &p0, float &p2, float &rad);

		glm::vec3 m_bound[2];//Objects BoundingBox in local space
		glm::vec3 m_sizeIndicator;
		glm::vec3 m_centroid;//Cylinders centroid in world space. Note: currently only correct if the base and apex radius are the same TODO: Calculate centroid with different apex radius

		glm::mat4 trans_w2o;//Transformation matrix world to object space
		glm::mat4 trans_o2w;//Transformation matrix object to world space

		KIRK::Material *m_material;

		glm::vec3 m_basepoint, m_apexpoint, m_Normal;//stored in world space
		float m_baseradius, m_apexradius, m_zmin, m_zmax, m_phiMax;//stored in local space
		int m_lA; // longest Axis: 0 = x-axis, 1 = y-axis, 2 = z-axis
		const float cCylinderEpsilon = 1e-7f;//Epsilon for Bounding Box

		bool planeBoxOverlap(const glm::vec3 &n, const glm::vec3 &v0, const glm::vec3 &halfsize) const;
	};

};

#include "KIRK/Common/Material.h"

#endif /* KIRK_OBJCYLINDER_H */