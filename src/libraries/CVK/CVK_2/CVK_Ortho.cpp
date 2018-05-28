#include "CVK_Ortho.h"

CVK::Ortho::Ortho(float ratio)
{
    m_left = -ratio;
    m_right = ratio;
    m_bottom = -1.0f;
    m_top = 1.0f;
    m_znear = 0.f;
    m_zfar = 10.f;
    setOrtho(m_left, m_right, m_bottom, m_top, m_znear, m_zfar);
}

CVK::Ortho::Ortho(float left, float right, float bottom, float top, float near, float far)
{
    m_left = left;
    m_right = right;
    m_bottom = bottom;
    m_top = top;
    m_znear = near;
    m_zfar = far;
    setOrtho(m_left, m_right, m_bottom, m_top, m_znear, m_zfar);
}

CVK::Ortho::~Ortho()
{
}

void CVK::Ortho::setOrtho(float left, float right, float bottom, float top, float near, float far)
{
    m_projection = glm::ortho(left, right, bottom, top, near, far);
}

void CVK::Ortho::setNearFar(float near, float far)
{
    m_znear = near;
    m_zfar = far;
    setOrtho(m_left, m_right, m_bottom, m_top, m_znear, m_zfar);
}

void CVK::Ortho::updateRatio(float ratio)
{
    m_left = -ratio / 2.f;
    m_right = ratio / 2.f;
    setOrtho(m_left, m_right, m_bottom, m_top, m_znear, m_zfar);
}
