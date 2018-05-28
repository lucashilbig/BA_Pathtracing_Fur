#include "CPU_NoDataStructure.h"

KIRK::CPU::NoDataStructure::NoDataStructure() : CPU_DataStructure()
{
}

KIRK::CPU::NoDataStructure::~NoDataStructure()
{
}

void KIRK::CPU::NoDataStructure::addBaseDataStructure(Scene *scene)
{
    if(m_hasCandidates)
    {
        m_minBound = glm::min(m_minBound, scene->getBounds()[0]);
        m_maxBound = glm::max(m_maxBound, scene->getBounds()[1]);
    } else
    {
        m_minBound = scene->getBounds()[0];
        m_maxBound = scene->getBounds()[1];
    }

    for(KIRK::Triangle *o : scene->getSceneObjects())
        addCandidate(o);
}

bool KIRK::CPU::NoDataStructure::closestIntersection(KIRK::Intersection *hit)
{
    return closestIntersectionWithCandidates(hit);
}

bool KIRK::CPU::NoDataStructure::isIntersection(KIRK::Ray *ray, float tMax)
{
    return isIntersectionWithCandidates(ray, tMax);
}
