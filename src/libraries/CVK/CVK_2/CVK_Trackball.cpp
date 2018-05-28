#include "CVK_Trackball.h"

CVK::Trackball::Trackball(int width, int height, CVK::Projection *projection)
{
    m_cameraPos = glm::vec3(0.0f, 0.0f, 5.0);
    m_center = glm::vec3(0.0f, 0.0f, 0.0f);
    m_up = glm::vec3(0.0f, 1.0f, 0.0f);

    m_width = width;
    m_height = height;

    m_sensitivity = 0.010f;
    m_theta = glm::pi<float>() / 2.0f;
    m_phi = 0.f;
    m_radius = 1.5;

    m_viewmatrix = glm::lookAt(m_center + m_cameraPos, m_center, m_up);

    m_oldX = width / 2.f;
    m_oldY = height / 2.f;

    setProjection(projection);
}

CVK::Trackball::Trackball(int width, int height)
{
    m_cameraPos = glm::vec3(0.0f, 0.0f, 5.0);
    m_center = glm::vec3(0.0f, 0.0f, 0.0f);
    m_up = glm::vec3(0.0f, 1.0f, 0.0f);

    m_width = width;
    m_height = height;

    m_sensitivity = 0.010f;
    m_theta = glm::pi<float>() / 2.0f;
    m_phi = 0.f;
    m_radius = 1.5;

    m_viewmatrix = glm::lookAt(m_center + m_cameraPos, m_center, m_up);

    m_oldX = width / 2.f;
    m_oldY = height / 2.f;

    CVK::Perspective *projection = new CVK::Perspective(width / (float)height);

    setProjection(projection);
}

CVK::Trackball::~Trackball()
{
}

void CVK::Trackball::update(GLFWwindow *window)
{
    double x, y;

    glfwGetCursorPos(window, &x, &y);

    if(glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS)
    {
        float changeX = ((float)x - m_oldX) * m_sensitivity;
        float changeY = ((float)y - m_oldY) * m_sensitivity;

        m_theta -= changeY;
        if(m_theta < 0.01f)
            m_theta = 0.01f;
        else if(m_theta > glm::pi<float>() - 0.01f)
            m_theta = glm::pi<float>() - 0.01f;

        m_phi -= changeX;
        if(m_phi < 0)
            m_phi += 2 * glm::pi<float>();
        else if(m_phi > 2 * glm::pi<float>())
            m_phi -= 2 * glm::pi<float>();
    }

    m_oldX = (float)x;
    m_oldY = (float)y;

    if(glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
        m_radius -= 0.1f;
    if(glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
        m_radius += 0.1f;
    if(m_radius < 0.1f)
        m_radius = 0.1f;

    m_cameraPos.x = m_center.x + m_radius * sin(m_theta) * sin(m_phi);
    m_cameraPos.y = m_center.y + m_radius * cos(m_theta);
    m_cameraPos.z = m_center.z + m_radius * sin(m_theta) * cos(m_phi);

    m_viewmatrix = glm::lookAt(m_cameraPos, m_center, m_up);
}

void CVK::Trackball::setCenter(glm::vec3 *center)
{
    m_center = *center;
}

void CVK::Trackball::setRadius(float radius)
{
    m_radius = radius;
}

void CVK::Trackball::setUpvector(glm::vec3 up)
{
    m_up = up;
}

void CVK::Trackball::setPosition(glm::vec3 pos)
{
    m_theta = acosf(glm::dot(glm::normalize(pos), m_up));
    m_phi = acosf(glm::dot(glm::normalize(glm::vec3(pos.x, 0.0001f, pos.z)), glm::vec3(0.0f, 0.0f, 1.0f)));
    if(pos.x < 0)
        m_phi *= -1.0f;

    if(m_theta < 0.01f)
        m_theta = 0.01f;
    else if(m_theta > glm::pi<float>() - 0.01f)
        m_theta = glm::pi<float>() - 0.01f;

    if(m_phi < 0)
        m_phi += 2 * glm::pi<float>();
    else if(m_phi > 2 * glm::pi<float>())
        m_phi -= 2 * glm::pi<float>();

    m_radius = glm::length(pos - m_center);
}
