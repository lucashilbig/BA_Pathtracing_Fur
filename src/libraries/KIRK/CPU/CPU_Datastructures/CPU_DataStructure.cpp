//
// Created by Maximilian Luzius on 11/05/16.
//

#include "CPU_DataStructure.h"


bool KIRK::CPU::CPU_DataStructure::testClosestIntersectionWithCandidates(std::vector<KIRK::Triangle *> *candidates,
                                                                      KIRK::Intersection *hit, float tMin, float tMax)
{
    bool intersectionFoundInLeaf = false;

    for(KIRK::Triangle *candidate : *candidates)
        if(candidate->closestIntersection(hit, tMin, tMax))
        {
            tMax = hit->m_lambda;
            intersectionFoundInLeaf |= true;
        };

    return intersectionFoundInLeaf;
}

bool KIRK::CPU::CPU_DataStructure::testIsIntersectionWithCandidates(std::vector < KIRK::Triangle * > *candidates, KIRK::Ray *ray,
                                                                 float tMax)
{
    for(KIRK::Triangle *candidate : *candidates)
        if(candidate->isIntersection(ray, tMax))
            return true;

    return false;
}
