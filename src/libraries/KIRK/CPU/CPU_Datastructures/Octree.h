#ifndef __CVK_DS_OCTREE_H__
#define __CVK_DS_OCTREE_H__

#include "CPU_DataStructure.h"
#include "Tree.h"
#include "Voxel.h"

namespace KIRK {
namespace CPU {

class Octree : public CPU_DataStructure, public Tree<Octree>
{
public:
    Octree(glm::vec3 minBound, glm::vec3 maxBound, int maxRecursionDepth);
    virtual ~Octree();

    void addBaseDataStructure(KIRK::CPU::Scene *scene);

    bool closestIntersection(KIRK::Intersection *hit);
    bool isIntersection(KIRK::Ray *ray, float tMax);

    virtual void printDebugInfo();
    virtual int getSizeInBytes();

private:
    Octree(int depth, glm::vec3 mid);
    void subdivide();
    void addObject(KIRK::Object *obj);

    bool traverseNode(glm::vec3 *tmin, glm::vec3 *tmax, unsigned char directionBits, KIRK::Intersection *hit);
    bool traverseNode(glm::vec3 *tmin, glm::vec3 *tmax, unsigned char directionBits, float tMax, KIRK::Ray *ray);

    unsigned char getFirstNode(glm::vec3 *tmin, glm::vec3 *tmid);
    unsigned char getNextNode(glm::vec3 *tmax, unsigned char nodeX, unsigned char nodeY, unsigned char nodeZ);

    glm::vec3 m_mid, m_halfSize;
    static int s_maxRecursionDepth;
    static glm::vec3 *s_NodeRadius;
    //glm::vec3 *m_bound;
};
}}

#endif //__CVK_DS_OCTREE_H__
