#ifndef __CVK_RT_OBJTRIANGLE_H
#define __CVK_RT_OBJTRIANGLE_H

#include "KIRK/Common/Object.h"

namespace KIRK {
	
/*
 * Every Tri consists of _a, _b, _c, where _a is always referring to the
 * vertex with smallest value on longest axis of triangle, _b refers to the
 * middle vertex and _c is always referring to the vertex with highest value
 * on longest axis. The longest axis is computed with the bb, which is stored
 * in bounds[]. Because of the ordering explained above, the .cpp file
 * is a mess... Furthermore m_ab referrs to the vector between _a and _b, m_bc
 * refers to the vector between _b and _c and m_ac refers to the vector between
 * _a and _c, which is the longest one. All these Attributes are needed for
 * faster voxelization.
*/
class Triangle : public KIRK::Object
{
public:
    Triangle();
    Triangle(glm::vec3 a, glm::vec3 b, glm::vec3 c, glm::vec3 na, glm::vec3 nb, glm::vec3 nc, glm::vec2 tca,
             glm::vec2 tcb, glm::vec2 tcc);
    ~Triangle();

    bool closestIntersection(Intersection *hit, float tMin, float tMax) override;
    bool isIntersection(Ray *ray, float tMax) override;
    void calcNormal(Intersection *hit) override;
    void calcTcoord(Intersection *hit) override;
    bool isInAABB(glm::vec3 *bbox) override;
	bool closestIntersectionAsPlane(Intersection *hit, float tMin, float tMax);

    glm::vec3 *getA() { return &m_A; }
    glm::vec3 *getB() { return &m_B; }
    glm::vec3 *getC() { return &m_C; }
    glm::vec3 *getVecAB() { return &m_ab; }//previously named u
    glm::vec3 *getVecAC() { return &m_ac; }//previously named v
    glm::vec3 *getVecBC() { return &m_bc; }//previously named uv
	
    int geti0() { return i0; }
    int geti1() { return i1; }
    int geti2() { return i2; }
	
private:
	const float cTriangleEpsilon = 1e-7f;

    void computeBounds() override;

	bool testAxis(const glm::vec3 e, const glm::vec3 v[], const glm::vec3 boxhalfsize, const int i, const int j, const int vi1, const int vi2, const int f, float &min, float &max, float &p0, float &p2, float &rad);
		
    glm::vec3 m_A, m_B, m_C, m_ab, m_ac, m_bc, m_Normal;
    glm::vec3 m_na, m_nb, m_nc;
    glm::vec2 m_tca, m_tcb, m_tcc;

    // Save original indices for triangle
    int i0 = 0;
    int i1 = 0;
    int i2 = 0;

    bool planeBoxOverlap(const glm::vec3 &n, const glm::vec3 &v0, const glm::vec3 &halfsize) const;
};

};

#endif /* __CVK_RT_OBJTRIANGLE_H */

