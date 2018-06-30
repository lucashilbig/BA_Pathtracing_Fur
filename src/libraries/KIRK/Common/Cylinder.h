#ifndef KIRK_OBJCYLINDER_H
#define KIRK_OBJCYLINDER_H

#include "KIRK/Common/Object.h"

namespace KIRK {

	/* Class for Cylinder Intersection used for Raytracing. In this implementation the Cylinder is in a 
	*  local object space and the ray will be transformed from world space to this local object space.
	*  Because Kirk works in world space, every public variable that goes in or out of this class will
	*  internal be transformed to local space for calculations and than the return value will be transformed
	*  back to world space.
	*  This class is a modified version of the CVK_RT_ObjCone.h class of CVK.
	*/
	class Cylinder : public KIRK::Object
	{
	public:
		Cylinder(glm::vec3 basePoint, glm::vec3 apexPoint, float baseRadius, float apexRadius, glm::mat4 *modelMatrix);

		~Cylinder();

		bool closestIntersection(Intersection *hit, float tMin, float tMax) override;
		bool isIntersection(Ray *ray, float tMax) override;
		void calcTcoord(Intersection *hit) override;
		void calcNormal(Intersection *hit) override;
		bool isInAABB(glm::vec3 *bbox) override;

		/*glm::vec3 *getBasePoint() { return &m_basepoint; }
		glm::vec3 *getApexPoint() { return &m_apexpoint; }
		glm::vec3 getU() { return &m_u; }
		glm::vec3 getV() { return &m_v; }
		glm::vec3 getW() { return &m_w; }
		float getBaseRadius() { return m_baseradius; }
		float getApexRadius() { return m_apexradius; }*/

	private:
		void computeBounds() override;//calculates Cylinders BoundingBox in WorldSpace
		
		glm::vec3 m_basepoint, m_apexpoint;//stored in world space
		glm::vec3 m_u, m_v, m_w;//The cones local space u,v,w axis
		float m_baseradius, m_apexradius, m_height;
		float m_slope, m_base_d, m_min_d, m_max_d;//cylinders slope and bottom and top caps
	};

};

#endif /* KIRK_OBJCYLINDER_H */