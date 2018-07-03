#ifndef __CVK_DS_CPU_DATASTRUCTURE_H
#define __CVK_DS_CPU_DATASTRUCTURE_H

#include "glm/gtc/type_precision.hpp"
#include "glm/gtx/component_wise.hpp"

#include "KIRK/CPU/CPU_Raytracer/CPU_Scene.h"

#include "Voxel.h"
#include "KIRK/Utils/Log.h"

namespace KIRK {
namespace CPU {

class CPU_DataStructure : public Voxel
{

public:
    CPU_DataStructure(glm::vec3 minBound, glm::vec3 maxBound) : Voxel(minBound, maxBound) {};

    CPU_DataStructure() : Voxel() {};

    virtual ~CPU_DataStructure() {};

    virtual void addBaseDataStructure(KIRK::CPU::Scene *scene) = 0;

    virtual bool closestIntersection(KIRK::Intersection *hit) = 0;
    virtual bool isIntersection(KIRK::Ray *ray, float tMax) = 0;

    virtual void printDebugInfo() {};

    virtual int getSizeInBytes()
    {
        return -1;
    };

protected:
    bool testClosestIntersectionWithCandidates(std::vector<KIRK::Object *> *candidates, KIRK::Intersection *hit, float tMin = 0.0f,
                                               float tMax = FLT_MAX);
    bool testIsIntersectionWithCandidates(std::vector<KIRK::Object *> *candidates, KIRK::Ray *ray, float tMax = FLT_MAX);

};
}}

#endif
