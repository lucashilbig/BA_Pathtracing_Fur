#include "CVK_PointList.h"

CVK::PointList::PointList()
{
    m_vertexbuffer = INVALID_OGL_VALUE;
    m_sizebuffer = INVALID_OGL_VALUE;
    m_colorbuffer = INVALID_OGL_VALUE;

    m_vao = INVALID_OGL_VALUE;
    m_geotype = CVK_POINTS;
}

CVK::PointList::~PointList()
{
}

void CVK::PointList::init_PointList()
{
    m_index.clear();
    m_vertices.clear();
    m_size.clear();
    m_color.clear();
}

void CVK::PointList::add_Point(glm::vec3 p, glm::vec3 size, glm::vec3 col)
{
    m_index.push_back(m_vertices.size());
    m_vertices.push_back(glm::vec4(p, 1.f));
    m_size.push_back(size);
    m_color.push_back(col);
}

void CVK::PointList::finish_PointList()
{
    m_indices = m_index.size();
    m_points = m_vertices.size();

    if(m_vao == INVALID_OGL_VALUE)
        glGenVertexArrays(1, &m_vao);
    glBindVertexArray(m_vao);

    if(m_vertexbuffer == INVALID_OGL_VALUE)
    {
        glGenBuffers(1, &m_vertexbuffer);
        glBindBuffer(GL_ARRAY_BUFFER, m_vertexbuffer);
        glEnableVertexAttribArray(VERTICES);
        glVertexAttribPointer(VERTICES, 4, GL_FLOAT, GL_FALSE, 0, 0);
    }
    if(m_sizebuffer == INVALID_OGL_VALUE)
    {
        glGenBuffers(1, &m_sizebuffer);
        glBindBuffer(GL_ARRAY_BUFFER, m_sizebuffer);
        glEnableVertexAttribArray(ATTRIBUTE_0);
        glVertexAttribPointer(ATTRIBUTE_0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    }
    if(m_colorbuffer == INVALID_OGL_VALUE)
    {
        glGenBuffers(1, &m_colorbuffer);
        glBindBuffer(GL_ARRAY_BUFFER, m_colorbuffer);
        glEnableVertexAttribArray(COLOUR);
        glVertexAttribPointer(COLOUR, 3, GL_FLOAT, GL_FALSE, 0, 0);
    }
    if(m_indexlist == INVALID_OGL_VALUE)
    {
        glGenBuffers(1, &m_indexlist);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_indexlist);
    }

    glBindBuffer(GL_ARRAY_BUFFER, m_vertexbuffer);
    glBufferData(GL_ARRAY_BUFFER, m_points * sizeof(glm::vec4), &m_vertices[0], GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, m_sizebuffer);
    glBufferData(GL_ARRAY_BUFFER, m_points * sizeof(glm::vec3), &m_size[0], GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, m_colorbuffer);
    glBufferData(GL_ARRAY_BUFFER, m_points * sizeof(glm::vec3), &m_color[0], GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_indexlist);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_indices * sizeof(unsigned int), &m_index[0], GL_STATIC_DRAW);

    glBindVertexArray(0);
    glDisableVertexAttribArray(VERTICES);
    glDisableVertexAttribArray(ATTRIBUTE_0);
    glDisableVertexAttribArray(COLOUR);
}

void CVK::PointList::render()
{
    glBindVertexArray(m_vao);
    glDrawElements(GL_POINTS, m_indices, GL_UNSIGNED_INT, 0);
}
