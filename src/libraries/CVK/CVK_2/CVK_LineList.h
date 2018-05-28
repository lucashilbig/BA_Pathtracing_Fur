#ifndef __CVK_LINELIST_H
#define __CVK_LINELIST_H

#include "CVK_Defs.h"
#include "CVK_Geometry.h"

namespace CVK {
class LineList : public CVK::Geometry
{
public:
    LineList();
    ~LineList();

    void init_LineList();
    void add_Line(glm::vec3 s, glm::vec3 e, glm::vec3 col);
    void finish_LineList();
    void render();

private:
    GLuint m_colorbuffer;
    std::vector<glm::vec3> m_color;
};

};

#endif /* __CVK_LINELIST_H */

