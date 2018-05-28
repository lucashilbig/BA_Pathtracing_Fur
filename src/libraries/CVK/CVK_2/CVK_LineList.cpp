#include "CVK_LineList.h"

CVK::LineList::LineList()
{
    m_vertexbuffer = INVALID_OGL_VALUE;
    m_colorbuffer = INVALID_OGL_VALUE;

    m_vao = INVALID_OGL_VALUE;
    m_geotype = CVK_LINES;
}

CVK::LineList::~LineList()
{
}

void CVK::LineList::init_LineList()
{
    m_index.clear();
    m_vertices.clear();
    m_color.clear();
}

void CVK::LineList::add_Line(glm::vec3 s, glm::vec3 e, glm::vec3 col)
{
    m_index.push_back(m_vertices.size());
    m_vertices.push_back(glm::vec4(s, 1.f));
    m_index.push_back(m_vertices.size());
    m_vertices.push_back(glm::vec4(e, 1.f));
    m_color.push_back(col);
    m_color.push_back(col);
}

void CVK::LineList::finish_LineList()
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

    if(m_points > 0)
    {
        glBindBuffer(GL_ARRAY_BUFFER, m_vertexbuffer);
        glBufferData(GL_ARRAY_BUFFER, m_points * sizeof(glm::vec4), &m_vertices[0], GL_STATIC_DRAW);

        glBindBuffer(GL_ARRAY_BUFFER, m_colorbuffer);
        glBufferData(GL_ARRAY_BUFFER, m_points * sizeof(glm::vec3), &m_color[0], GL_STATIC_DRAW);
    }

    if(m_indices > 0)
    {
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_indexlist);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_indices * sizeof(unsigned int), &m_index[0], GL_STATIC_DRAW);
    }

    glBindVertexArray(0);
    glDisableVertexAttribArray(VERTICES);
    glDisableVertexAttribArray(ATTRIBUTE_0);
    glDisableVertexAttribArray(COLOUR);
}

void CVK::LineList::render()
{
    glBindVertexArray(m_vao);
    glDrawElements(GL_LINES, m_indices, GL_UNSIGNED_INT, 0);
}
