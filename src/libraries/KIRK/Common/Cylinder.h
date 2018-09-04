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
		
		float getBaseRadius() const { return m_baseradius; }
		float getApexRadius() const { return m_apexradius; }

		/*glm::vec3 *getBasePoint() { return &m_basepoint; }
		glm::vec3 *getApexPoint() { return &m_apexpoint; }		
		*/

	private:
		void computeBounds() override;//calculates Cylinders BoundingBox in WorldSpace
		
		glm::vec3 m_basepoint, m_apexpoint;//stored in world space
		float m_baseradius, m_apexradius, m_height;
		float m_slope, m_base_d, m_min_d, m_max_d;//cylinders slope and bottom and top caps
	};

};

#endif /* KIRK_OBJCYLINDER_H */