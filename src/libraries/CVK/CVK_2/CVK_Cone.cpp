#include "CVK_Cone.h"
#include "KIRK/Common/Ray.h"

CVK::Cone::Cone()
{
    m_baseradius = 1.0f;
    m_apexradius = 1.0f;
    m_basepoint = glm::vec3(0);
    m_apexpoint = glm::vec3(0.f, 1.f, 0.f);

    m_resolution = 20;
    m_geotype = CVK_CONE;
    create();
}

CVK::Cone::Cone(float baseradius, float apexradius, float height, int resolution)
{
    m_baseradius = baseradius;
    m_apexradius = apexradius;
    m_basepoint = glm::vec3(0);
    m_apexpoint = glm::vec3(0.f, height, 0.f);

    m_resolution = resolution;
    m_geotype = CVK_CONE;
    create();
}

CVK::Cone::Cone(glm::vec3 basepoint, glm::vec3 apexpoint, float baseradius, float apexradius, int resolution)
{
    m_baseradius = baseradius;
    m_apexradius = apexradius;
    m_basepoint = basepoint;
    m_apexpoint = apexpoint;

    m_resolution = resolution;
    m_geotype = CVK_CONE;
    create();
}

CVK::Cone::~Cone()
{
}

glm::vec3 *CVK::Cone::getBasepoint()
{
    return &m_basepoint;
}

glm::vec3 *CVK::Cone::getApexpoint()
{
    return &m_apexpoint;
}

glm::vec3 *CVK::Cone::get_u()
{
    return &m_u;
}

glm::vec3 *CVK::Cone::get_v()
{
    return &m_v;
}

glm::vec3 *CVK::Cone::get_w()
{
    return &m_w;
}

float CVK::Cone::getBaseradius()
{
    return m_baseradius;
}

float CVK::Cone::getSlope()
{
    return m_slope;
}

float CVK::Cone::getApexradius()
{
    return m_apexradius;
}

void CVK::Cone::create()
{
    // iniatialize the variable we are going to use
    float u, v;
    float radius, phi;
    glm::vec3 q;
    int offset = 0, i, j;
    glm::vec3 n1, n2, n;

    m_v = m_apexpoint - m_basepoint;
    m_height = glm::length(m_v);
    m_v = glm::normalize(m_v);

    /* find two axes which are at right angles to cone_v */
    glm::vec3 tmp(0.f, 1.f, 0.f);
    if(1.f - fabs(glm::dot(tmp, m_v)) < KIRK::cRayEpsilon)
        tmp = glm::vec3(0.f, 0.f, 1.f);

    m_u = glm::normalize(glm::cross(m_v, tmp));
    m_w = glm::normalize(glm::cross(m_u, m_v));

    m_slope = (m_baseradius - m_apexradius) / m_height;

    // Envelope
    for(j = 0; j <= m_resolution; j++)  //radius
        for(i = 0; i <= m_resolution; i++) //phi
        {
            u = i / (float)m_resolution;
            phi = 2 * glm::pi<float>() * u;
            v = j / (float)m_resolution;
            v = m_height * v;

            radius = m_baseradius - m_slope * v;
            q = m_basepoint + radius * sinf(phi) * m_u + v * m_v + radius * cosf(phi) * m_w;

            float t = glm::dot(q, m_v) - glm::dot(m_basepoint, m_v);
            glm::vec3 q_1 = q - t * m_v;
            glm::vec3 n = glm::normalize(q_1 - m_basepoint);
            n = glm::normalize(n + m_slope * m_v);

            m_vertices.push_back(glm::vec4(q, 1.0f));
            m_normals.push_back(n);
            m_uvs.push_back(glm::vec2(u, v / m_height));
        }

    m_points = m_vertices.size();

    // create index list
    for(j = 0; j < m_resolution; j++)
    {
        for(i = 0; i < m_resolution; i++)
        {
            // 1. Triangle
            m_index.push_back(offset + i + m_resolution + 1);
            m_index.push_back(offset + i);
            m_index.push_back(offset + i + 1);

            // 2. Triangle
            m_index.push_back(offset + i + 1);
            m_index.push_back(offset + i + m_resolution + 1 + 1);
            m_index.push_back(offset + i + m_resolution + 1);
        }
        offset += m_resolution + 1;
    }
    m_indices = m_index.size();

    createBuffers();
}
