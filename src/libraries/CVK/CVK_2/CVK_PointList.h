#ifndef __CVK_POINTLIST_H
#define __CVK_POINTLIST_H

#include "CVK_Defs.h"
#include "CVK_Geometry.h"

namespace CVK {
class PointList : public CVK::Geometry
{
public:
    PointList();
    ~PointList();

    void init_PointList();
    void add_Point(glm::vec3 p, glm::vec3 size, glm::vec3 col);
    void finish_PointList();
    void render();

private:
    GLuint m_sizebuffer;
    GLuint m_colorbuffer;
    std::vector<glm::vec3> m_size;
    std::vector<glm::vec3> m_color;
};

};

#endif /* __CVK_POINTLIST_H */

