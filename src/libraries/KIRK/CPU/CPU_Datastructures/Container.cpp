#include "Container.h"

KIRK::CPU::Container::Container()
{
    m_hasCandidates = false;
}

KIRK::CPU::Container::~Container()
{

}

bool KIRK::CPU::Container::closestIntersectionWithCandidates(KIRK::Intersection *hit, float tMin, float tMax)
{
    bool intersectionFoundInLeaf = false;

    for(KIRK::Triangle *candidate : m_candidateList)
        if(candidate->closestIntersection(hit, tMin, tMax))
        {
            tMax = hit->m_lambda;
            intersectionFoundInLeaf |= true;
        };

    return intersectionFoundInLeaf;
}

bool KIRK::CPU::Container::isIntersectionWithCandidates(KIRK::Ray *ray, float tMax)
{
    for(KIRK::Triangle *candidate : m_candidateList)
        if(candidate->isIntersection(ray, tMax))
            return true;

    return false;
}
