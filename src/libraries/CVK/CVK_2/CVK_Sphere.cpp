#include "CVK_Sphere.h"

CVK::Sphere::Sphere()
{
    m_radius = 1.0f;
    m_resolution = 20;
    m_Center = glm::vec3(0);
    create(m_Center, m_radius, m_resolution);
    m_geotype = CVK_SPHERE;
}

CVK::Sphere::Sphere(float radius)
{
    m_radius = radius;
    m_resolution = 20;
    m_Center = glm::vec3(0);
    create(m_Center, m_radius, m_resolution);
    m_geotype = CVK_SPHERE;
}

CVK::Sphere::Sphere(glm::vec3 center, float radius)
{
    m_radius = radius;
    m_resolution = 20;
    m_Center = center;
    create(m_Center, m_radius, m_resolution);
    m_geotype = CVK_SPHERE;
}

CVK::Sphere::Sphere(float radius, int resolution)
{
    m_radius = radius;
    m_resolution = resolution;
    m_Center = glm::vec3(0);
    create(m_Center, m_radius, m_resolution);
    m_geotype = CVK_SPHERE;
}

CVK::Sphere::Sphere(glm::vec3 center, float radius, int resolution)
{
    m_radius = radius;
    m_resolution = resolution;
    m_Center = center;
    create(m_Center, m_radius, m_resolution);
    m_geotype = CVK_SPHERE;
}

CVK::Sphere::~Sphere()
{
}

glm::vec3 *CVK::Sphere::get_center()
{
    return &m_Center;
}

float CVK::Sphere::get_radius()
{
    return m_radius;
}

void CVK::Sphere::create(glm::vec3 center, float radius, int resolution)
{
    // iniatialize the variable we are going to use
    float u, v;
    float phi, theta;
    float x, y, z;
    int offset = 0, i, j;

    // create points
    for(j = 0; j <= resolution; j++)  //theta
        for(i = 0; i <= resolution; i++) //phi
        {
            u = i / (float)resolution;
            phi = 2 * glm::pi<float>() * u;
            v = j / (float)resolution;
            theta = glm::pi<float>() * v;

            x = m_Center.x + radius * sin(theta) * sin(phi);
            y = m_Center.y + radius * cos(theta);
            z = m_Center.z + radius * sin(theta) * cos(phi);

            m_vertices.push_back(glm::vec4(x, y, z, 1.0f));
            m_normals.push_back(glm::vec3(x, y, z) / radius);
            m_uvs.push_back(glm::vec2(u, 1 - v));
        }
    m_points = m_vertices.size();

    // create index list
    for(j = 0; j < resolution; j++)
    {
        for(i = 0; i < resolution; i++)
        {
            // 1. Triangle
            m_index.push_back(offset + i);
            m_index.push_back(offset + i + resolution + 1);
            m_index.push_back(offset + i + resolution + 1 + 1);

            // 2. Triangle
            m_index.push_back(offset + i + resolution + 1 + 1);
            m_index.push_back(offset + i + 1);
            m_index.push_back(offset + i);
        }
        offset += resolution + 1;
    }
    m_indices = m_index.size();

    createBuffers();
}
