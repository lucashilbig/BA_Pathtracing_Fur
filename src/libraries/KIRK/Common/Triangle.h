#ifndef __CVK_RT_OBJTRIANGLE_H
#define __CVK_RT_OBJTRIANGLE_H

#include "KIRK/Common/Intersection.h"

namespace KIRK {

	class Material;

/*
 * Every Tri consists of _a, _b, _c, where _a is always referring to the
 * vertex with smallest value on longest axis of triangle, _b refers to the
 * middle vertex and _c is always referring to the vertex with highest value
 * on longest axis. The longest axis is computed with the bb, which is stored
 * in bounds[]. Because of the ordering explained above, the .cpp file
 * is a mess... Furthermore m_u referrs to the vector between _a and _b, m_uv
 * refers to the vector between _b and _c and m_v refers to the vector between
 * _a and _c, which is the longest one. All these Attributes are needed for
 * faster voxelization.
*/
class Triangle
{
public:
    Triangle();
    Triangle(glm::vec3 a, glm::vec3 b, glm::vec3 c, glm::vec3 na, glm::vec3 nb, glm::vec3 nc, glm::vec2 tca,
             glm::vec2 tcb, glm::vec2 tcc);
    ~Triangle();

    bool closestIntersection(Intersection *hit, float tMin, float tMax);
    bool isIntersection(Ray *ray, float tMax);
    bool closestIntersectionAsPlane(Intersection *hit, float tMin, float tMax);
    void calcNormal(Intersection *hit);
    void calcTcoord(Intersection *hit);
    bool isInAABB(glm::vec3 *bbox);

    glm::vec3 *getA() { return &m_A; }

    glm::vec3 *getB() { return &m_B; }

    glm::vec3 *getC() { return &m_C; }

    glm::vec3 *getU() { return &m_u; }

    glm::vec3 *getV() { return &m_v; }

    glm::vec3 *getUV() { return &m_uv; }

    glm::vec3 *getVecAB() { return &m_u; }

    glm::vec3 *getVecAC() { return &m_v; }

    glm::vec3 *getVecBC() { return &m_uv; }

    int getLongestAxis() { return m_lA; }

    int geti0() { return i0; }

    int geti1() { return i1; }

    int geti2() { return i2; }

	inline void setMaterial(KIRK::Material *material){
	    m_material = material;
	}
	inline KIRK::Material *getMaterial(){
	    return m_material;
	}
	
	glm::vec3 *getBounds();
	
	glm::vec3 getMinBound() const { return m_bound[0]; }
	
	glm::vec3 getMaxBound() const { return m_bound[1]; }
	
	glm::vec3 getCentroid() const { return m_centroid; }

private:
	const float cTriangleEpsilon = 1e-7f;

    void computeBounds();

	bool testAxis(const glm::vec3 e, const glm::vec3 v[], const glm::vec3 boxhalfsize, const int i, const int j, const int vi1, const int vi2, const int f, float &min, float &max, float &p0, float &p2, float &rad);

	glm::vec3 m_bound[2];
	glm::vec3 m_sizeIndicator;
	glm::vec3 m_centroid;

	KIRK::Material *m_material;

    glm::vec3 m_A, m_B, m_C, m_u, m_v, m_uv, m_Normal;
    glm::vec3 m_na, m_nb, m_nc;
    glm::vec2 m_tca, m_tcb, m_tcc;
    int m_lA; // longest Axis

    // Save original indices for triangle
    int i0 = 0;
    int i1 = 0;
    int i2 = 0;

    bool planeBoxOverlap(const glm::vec3 &n, const glm::vec3 &v0, const glm::vec3 &halfsize) const;
};

};

#include "KIRK/Common/Material.h"

#endif /* __CVK_RT_OBJTRIANGLE_H */

