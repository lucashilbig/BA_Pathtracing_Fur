#ifndef __CVK_DS_UNIFORM_GRID_H__
#define __CVK_DS_UNIFORM_GRID_H__

#include <vector>

#include "CPU_DataStructure.h"

namespace KIRK {
namespace CPU {

/*** WARNING: THIS MIGHT NOT WORK. The UniformGrid was not yet modified to work with our new Pathtracer / contains known bugs. Use at own risk. */
class UniformGrid : public CPU_DataStructure
{
public:
    //UniformGrid(KIRK::Scene *scene, float sizeOfVoxel);
    //UniformGrid(KIRK::Scene *scene, glm::vec3 sizeOfVoxel);
    UniformGrid(glm::vec3 minBound, glm::vec3 maxBound, int numVoxelsPerAxis);
    virtual ~UniformGrid();

    void addBaseDataStructure(KIRK::CPU::Scene *scene);

    bool closestIntersection(KIRK::Intersection *hit);
    bool isIntersection(KIRK::Ray *ray, float tMax);

    inline unsigned int getIndex(glm::ivec3 v)
    {
        return v.z * m_numVoxels.x * m_numVoxels.y + v.y * m_numVoxels.x + v.x;
    }

    inline unsigned int getIndex(int x, int y, int z)
    {
        return z * m_numVoxels.x * m_numVoxels.y + y * m_numVoxels.x + x;
    }

    virtual void printDebugInfo();
    virtual int getSizeInBytes();

private:
    void init();

    glm::vec3 m_voxelSize, m_invVoxelSize;
    glm::ivec3 m_numVoxels, m_maxIndex;
    unsigned int m_totalVoxels;

    std::vector<std::vector<KIRK::Object *>> m_voxelCandidates;
};

}}

#endif //__CVK_DS_UNIFORM_GRID_H__
