#include "CVK_Camera.h"

CVK::Camera::Camera()
{
    m_focusDistance = 4.5f;
    m_blurRadius = 1.0f;
}

CVK::Camera::~Camera()
{
}

glm::mat4 CVK::Camera::getView()
{
    return m_viewmatrix;
}

void CVK::Camera::getView(glm::vec3 *x, glm::vec3 *y, glm::vec3 *z, glm::vec3 *pos)
{
    *x = glm::vec3(glm::row(m_viewmatrix, 0));
    *y = glm::vec3(glm::row(m_viewmatrix, 1));
    *z = glm::vec3(glm::row(m_viewmatrix, 2));
    *pos = glm::vec3(glm::column(m_viewmatrix, 3));
    glm::mat3 mat_inv = glm::inverse(glm::mat3(m_viewmatrix));
    *pos = -mat_inv * *pos;
}

void CVK::Camera::setView(glm::mat4 view)
{
    m_viewmatrix = view;
}

void CVK::Camera::setWidthHeight(int width, int height)
{
	m_width = width;
	m_height = height;
	m_projection->updateRatio((float)width / (float)height);
}

void CVK::Camera::getWidthHeight(int *width, int *height)
{
    *width = m_width;
    *height = m_height;
}

void CVK::Camera::lookAt(glm::vec3 position, glm::vec3 center, glm::vec3 up)
{
    m_viewmatrix = glm::lookAt(position, center, up);
}

void CVK::Camera::setProjection(CVK::Projection *projection)
{
    m_projection = projection;
}

CVK::Projection *CVK::Camera::getProjection()
{
    return m_projection;
}

float CVK::Camera::getFocus()
{
    return m_focusDistance;
}

float CVK::Camera::getBlurRadius()
{
    return m_blurRadius;
}

void CVK::Camera::setFocus(float input)
{
    m_focusDistance = input;
}

void CVK::Camera::setBlurRadius(float input)
{
    m_blurRadius = input;
}
