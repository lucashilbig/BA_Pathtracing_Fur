#include "CVK_Triangle.h"

glm::vec3 def_A(0.f, 1.f, 0.f);
glm::vec3 def_B(-1.f, -1.f, 0.f);
glm::vec3 def_C(1.f, -1.f, 0.f);

glm::vec3 def_n(0.f, 0.f, 1.f);

glm::vec2 def_ta(0.f, 1.f);
glm::vec2 def_tb(0.f, 0.f);
glm::vec2 def_tc(1.f, 0.f);

CVK::Triangle::Triangle()
{
    create(def_A, def_B, def_C, def_n, def_n, def_n, def_ta, def_tb, def_tc);
    m_geotype = CVK_TRIANGLE;
}

CVK::Triangle::Triangle(glm::vec3 a, glm::vec3 b, glm::vec3 c)
{
    create(a, b, c, def_ta, def_tb, def_tc);
    m_geotype = CVK_TRIANGLE;
}

CVK::Triangle::Triangle(glm::vec3 a, glm::vec3 b, glm::vec3 c, glm::vec3 na, glm::vec3 nb, glm::vec3 nc)
{
    create(a, b, c, na, nb, nc, def_ta, def_tb, def_tc);
    m_geotype = CVK_TRIANGLE;
}

CVK::Triangle::Triangle(glm::vec3 a, glm::vec3 b, glm::vec3 c, glm::vec2 tca, glm::vec2 tcb, glm::vec2 tcc)
{
    create(a, b, c, tca, tcb, tcc);
    m_geotype = CVK_TRIANGLE;
}


CVK::Triangle::Triangle(glm::vec3 a, glm::vec3 b, glm::vec3 c, glm::vec3 na, glm::vec3 nb, glm::vec3 nc, glm::vec2 tca,
                        glm::vec2 tcb, glm::vec2 tcc)
{
    create(a, b, c, na, nb, nc, tca, tcb, tcc);
    m_geotype = CVK_TRIANGLE;
}

CVK::Triangle::~Triangle()
{
}

void CVK::Triangle::create(glm::vec3 a, glm::vec3 b, glm::vec3 c, glm::vec2 tca, glm::vec2 tcb, glm::vec2 tcc)
{
    glm::vec3 n = glm::normalize(glm::cross(b - a, c - a));
    create(a, b, c, n, n, n, tca, tcb, tcc);
}

void
CVK::Triangle::create(glm::vec3 a, glm::vec3 b, glm::vec3 c, glm::vec3 na, glm::vec3 nb, glm::vec3 nc, glm::vec2 tca,
                      glm::vec2 tcb, glm::vec2 tcc)
{
    m_vertices.push_back(glm::vec4(a, 1.0f));
    m_vertices.push_back(glm::vec4(b, 1.0f));
    m_vertices.push_back(glm::vec4(c, 1.0f));

    m_normals.push_back(na);
    m_normals.push_back(nb);
    m_normals.push_back(nc);

    m_uvs.push_back(tca);
    m_uvs.push_back(tcb);
    m_uvs.push_back(tcc);

    m_points = 3;

    m_index.push_back(0);
    m_index.push_back(1);
    m_index.push_back(2);

    m_indices = 3;

    createBuffers();
}

void CVK::Triangle::set_Points(glm::vec3 a, glm::vec3 b, glm::vec3 c)
{
    m_vertices.clear();
    m_vertices.push_back(glm::vec4(a, 1.0f));
    m_vertices.push_back(glm::vec4(b, 1.0f));
    m_vertices.push_back(glm::vec4(c, 1.0f));

    glBindBuffer(GL_ARRAY_BUFFER, m_vertexbuffer);
    glBufferData(GL_ARRAY_BUFFER, m_points * sizeof(glm::vec4), &m_vertices[0], GL_STATIC_DRAW);
}

void CVK::Triangle::get_Points(glm::vec3 *a, glm::vec3 *b, glm::vec3 *c)
{
    *a = glm::vec3(m_vertices[0]);
    *b = glm::vec3(m_vertices[1]);
    *c = glm::vec3(m_vertices[2]);
}

void CVK::Triangle::set_Normal(glm::vec3 n)
{
    m_normals.clear();
    m_normals.push_back(n);
    m_normals.push_back(n);
    m_normals.push_back(n);

    glBindBuffer(GL_ARRAY_BUFFER, m_normalbuffer);
    glBufferData(GL_ARRAY_BUFFER, m_points * sizeof(glm::vec3), &m_normals[0], GL_STATIC_DRAW);
}

void CVK::Triangle::get_Normal(glm::vec3 *n)
{
    *n = m_normals[0];
}

void CVK::Triangle::set_Normals(glm::vec3 na, glm::vec3 nb, glm::vec3 nc)
{
    m_normals.clear();
    m_normals.push_back(na);
    m_normals.push_back(nb);
    m_normals.push_back(nc);

    glBindBuffer(GL_ARRAY_BUFFER, m_normalbuffer);
    glBufferData(GL_ARRAY_BUFFER, m_points * sizeof(glm::vec3), &m_normals[0], GL_STATIC_DRAW);
}

void CVK::Triangle::get_Normals(glm::vec3 *na, glm::vec3 *nb, glm::vec3 *nc)
{
    *na = m_normals[0];
    *nb = m_normals[1];
    *nc = m_normals[2];
}

void CVK::Triangle::set_Tcoords(glm::vec2 tca, glm::vec2 tcb, glm::vec2 tcc)
{
    m_uvs.clear();
    m_uvs.push_back(tca);
    m_uvs.push_back(tcb);
    m_uvs.push_back(tcc);

    glBindBuffer(GL_ARRAY_BUFFER, m_uvbuffer);
    glBufferData(GL_ARRAY_BUFFER, m_points * sizeof(glm::vec2), &m_uvs[0], GL_STATIC_DRAW);
}

void CVK::Triangle::get_Tcoords(glm::vec2 *tca, glm::vec2 *tcb, glm::vec2 *tcc)
{
    *tca = m_uvs[0];
    *tcb = m_uvs[1];
    *tcc = m_uvs[2];
}
