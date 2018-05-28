#include "Voxel.h"

KIRK::CPU::Voxel::Voxel(glm::vec3 minBound, glm::vec3 maxBound) : Container()
{
    m_minBound = minBound;
    m_maxBound = maxBound;
    m_size = m_maxBound - m_minBound;
}

KIRK::CPU::Voxel::~Voxel()
{

}
