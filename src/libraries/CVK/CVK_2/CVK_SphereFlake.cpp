#include "CVK_SphereFlake.h"

#include "KIRK/Utils/Log.h"

std::shared_ptr<CVK::Geometry> CVK::SphereFlake::s_sphere = NULL;

CVK::SphereFlake::SphereFlake(int recursionDepth, std::shared_ptr<Material> mat, std::shared_ptr<Node> scene)
{
    glm::mat4 modelmatrix = glm::mat4(1.f);
    m_mat = mat;
    m_scene = scene;
    float radius = 0.5f;

    if(s_sphere == NULL)
        s_sphere = std::make_shared<Sphere>(radius, 20);

    generateSphereFlake(0, recursionDepth, modelmatrix, radius);

    int n = 0;
    for(float i = 0.f; i <= recursionDepth; i++)
        n += (int)glm::pow(9.f, i);
    LOG_INFO("Sphere flake with % spheres", n);
}

CVK::SphereFlake::SphereFlake(int recursionDepth, std::shared_ptr<Material> mat, std::shared_ptr<Node> scene, glm::mat4 ModelMatrix, float radius)
{
    m_mat = mat;
    m_scene = scene;

    if(s_sphere == NULL)
        s_sphere = std::make_shared<Sphere>(radius, 20);

    generateSphereFlake(0, recursionDepth, ModelMatrix, radius);

    int n = 0;
    for(float i = 0.f; i <= recursionDepth; i++)
        n += (int)glm::pow(9.f, i);
    LOG_INFO("Sphere flake with % spheres", n);
}

CVK::SphereFlake::~SphereFlake()
{
}

void CVK::SphereFlake::generateSphereFlake(int ctr, int recursionDepth, glm::mat4 ModelMatrix, float radius)
{
    //CVK::Sphere* sphere = new CVK::Sphere( radius, glm::max( 5, 20-ctr*5));
    char string[100];
    sprintf(string, "Sphere %d", ctr);
	std::shared_ptr<CVK::Node> sphere_node = std::make_shared<CVK::Node>(string);
    sphere_node->setMaterial(m_mat);
    sphere_node->setGeometry(s_sphere);
    sphere_node->setModelMatrix(glm::scale(ModelMatrix, glm::pow(glm::vec3(1 / 3.0f), glm::vec3((float)ctr))));
    m_scene->addChild(sphere_node);

    if(recursionDepth-- > 0)
    {
        float newRadius = radius / 3.f;

        glm::vec3 xVec = glm::vec3(glm::column(ModelMatrix, 0));
        glm::vec3 yVec = glm::vec3(glm::column(ModelMatrix, 1));
        glm::vec3 zVec = glm::vec3(glm::column(ModelMatrix, 2));
        glm::vec3 center = glm::vec3(glm::column(ModelMatrix, 3));

        float rho = 2.f * glm::pi<float>() / 6.f;

        for(int i = 0; i < 6; i++)
        {
            float phi = i * rho;

            glm::vec3 newy = ((float)cos(phi)) * xVec + ((float)sin(phi)) * zVec;
            glm::vec3 newx = glm::cross(newy, yVec);
            glm::vec3 newz = glm::cross(newx, newy);
            glm::vec3 newcenter = center + newy * (radius + newRadius);

            glm::mat4 newmat = glm::mat4(glm::vec4(newx, 0), glm::vec4(newy, 0), glm::vec4(newz, 0),
                                         glm::vec4(newcenter, 1.f));

            generateSphereFlake(ctr + 1, recursionDepth, newmat, newRadius);
        }

        for(int i = 0; i < 3; i++)
        {
            float phi = (i + 0.5f) * 2.f * rho;

            glm::vec3 tmpy = ((float)cos(phi)) * xVec + ((float)sin(phi)) * zVec;
            glm::vec3 newy = ((float)cos(rho)) * tmpy + ((float)sin(rho)) * yVec;
            glm::vec3 newx = glm::cross(newy, yVec);
            newx = glm::normalize(newx);
            glm::vec3 newz = glm::cross(newx, newy);
            glm::vec3 newcenter = center + newy * (radius + newRadius);

            glm::mat4 newmat = glm::mat4(glm::vec4(newx, 0), glm::vec4(newy, 0), glm::vec4(newz, 0),
                                         glm::vec4(newcenter, 1.f));

            generateSphereFlake(ctr + 1, recursionDepth, newmat, newRadius);
        }
    }
}
