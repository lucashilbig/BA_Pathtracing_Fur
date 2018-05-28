#ifndef __CVK_DS_VOXEL_H__
#define __CVK_DS_VOXEL_H__

#include <glm/glm.hpp>

#include "Container.h"

namespace KIRK {
namespace CPU {

class Voxel : public Container
{
public:
    Voxel(glm::vec3 minBound = glm::vec3(0.0f), glm::vec3 maxBound = glm::vec3(0.0f));
    virtual ~Voxel();

    inline const glm::vec3 *getMinBound() { return &m_minBound; }

    inline const glm::vec3 *getMaxBound() { return &m_maxBound; }

    inline void setMinBound(glm::vec3 bound)
    {
        m_minBound = bound;
        m_size = m_maxBound - m_minBound;
    }

    inline void setMaxBound(glm::vec3 bound)
    {
        m_maxBound = bound;
        m_size = m_maxBound - m_minBound;
    }

    inline void setBounds(const glm::vec3 bounds[2])
    {
        m_minBound = bounds[0];
        m_maxBound = bounds[1];
        m_size = m_maxBound - m_minBound;
    }

protected:
    glm::vec3 m_minBound, m_maxBound, m_size;
};

}}

#endif //__CVK_DS_VOXEL_H__
