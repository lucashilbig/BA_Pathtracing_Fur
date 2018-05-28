#include "CVK_Perspective.h"

CVK::Perspective::Perspective(float ratio)
{
    m_fov = glm::radians(60.f);
    m_ratio = ratio;

    m_znear = 0.001f;
    m_zfar = 10.f;
    setPerspective(m_fov, m_ratio, m_znear, m_zfar);
}

CVK::Perspective::Perspective(float fov, float ratio, float near, float far)
{
    m_fov = fov;
    m_ratio = ratio;

    m_znear = near;
    m_zfar = far;
    setPerspective(m_fov, m_ratio, m_znear, m_zfar);
}

CVK::Perspective::~Perspective()
{
}

void CVK::Perspective::setPerspective(float fov, float ratio, float near, float far)
{
    m_projection = glm::perspective(m_fov, ratio, m_znear, m_zfar);
}

void CVK::Perspective::setNearFar(float near, float far)
{
    m_znear = near;
    m_zfar = far;
    setPerspective(m_fov, m_ratio, m_znear, m_zfar);
}


void CVK::Perspective::setFov(float fov)
{
    m_fov = fov;
    setPerspective(m_fov, m_ratio, m_znear, m_zfar);
}

float CVK::Perspective::getFov()
{
    return m_fov;
}

void CVK::Perspective::updateRatio(float ratio)
{
    m_ratio = ratio;
    setPerspective(m_fov, m_ratio, m_znear, m_zfar);
}