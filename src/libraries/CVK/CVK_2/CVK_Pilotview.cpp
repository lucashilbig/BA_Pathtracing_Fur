#include "CVK_Pilotview.h"

CVK::Pilotview::Pilotview(int width, int height, CVK::Projection *projection)
{
    m_cameraPos = glm::vec3(0.0f, 0.0f, 5.0f);
    m_direction = glm::vec3(0.0f, 0.0f, -1.0f);
    m_up = glm::vec3(0.0f, 1.0f, 0.0f);

    m_width = width;
    m_height = height;

    m_step = 0.0f;
    m_speed = 0.0f;
    m_sensitivity = 0.01f;
    m_theta = glm::pi<float>() / 2.0f;
    m_phi = glm::pi<float>();

    m_viewmatrix = glm::lookAt(m_cameraPos, m_cameraPos + m_direction, m_up);

    m_oldX = width / 2.f;
    m_oldY = height / 2.f;

    setProjection(projection);
}

CVK::Pilotview::Pilotview(int width, int height)
{
    m_cameraPos = glm::vec3(0.0f, 0.0f, 5.0f);
    m_direction = glm::vec3(0.0f, 0.0f, -1.0f);
    m_up = glm::vec3(0.0f, 1.0f, 0.0f);

    m_width = width;
    m_height = height;

    m_step = 0.0f;
    m_speed = 0.0f;
    m_sensitivity = 0.01f;
    m_theta = glm::pi<float>() / 2.0f;
    m_phi = glm::pi<float>();

    m_viewmatrix = glm::lookAt(m_cameraPos, m_cameraPos + m_direction, m_up);

    m_oldX = width / 2.f;
    m_oldY = height / 2.f;

    CVK::Perspective *projection = new CVK::Perspective(width / (float)height);

    setProjection(projection);
}

CVK::Pilotview::~Pilotview()
{
}

void CVK::Pilotview::update(GLFWwindow *window)
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

        m_oldX = (float)x;
        m_oldY = (float)y;

        m_direction.x = sin(m_theta) * sin(m_phi);
        m_direction.y = cos(m_theta);
        m_direction.z = sin(m_theta) * cos(m_phi);

    } else
    {
        m_oldX = (float)x;
        m_oldY = (float)y;
    }

    if(glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
        m_speed += 0.00001f;
    if(glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
        m_speed -= 0.00001f;
    if(glfwGetKey(window, ' ') == GLFW_PRESS)
        m_speed = 0.0f;

    m_cameraPos += m_speed * m_direction;
    m_viewmatrix = glm::lookAt(m_cameraPos, m_cameraPos + m_direction, m_up);
}

void CVK::Pilotview::setPosition(glm::vec3 pos)
{
    m_cameraPos = pos;
}

void CVK::Pilotview::setLookAt(glm::vec3 *position, glm::vec3 *center)
{
    m_cameraPos = *position;
    glm::vec3 diff = glm::normalize(*center - *position);
    m_direction = diff;

    m_theta = acos(diff.y);
    if(m_theta < 0.01f)
        m_theta = 0.01f;
    m_phi = acos(diff.z / (sin(m_theta) + 0.0001f));
    if(diff.x > 0)
        m_phi = 2 * glm::pi<float>() - m_phi;
}

void CVK::Pilotview::setUpvector(glm::vec3 *up)
{
    m_up = *up;
}