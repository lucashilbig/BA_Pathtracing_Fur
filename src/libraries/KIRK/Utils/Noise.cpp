#include "Noise.h"

KIRK::Noise::Noise()
{
    m_fbo = std::unique_ptr<CVK::FBO>(new CVK::FBO(512, 512, 4, true, false));

    const char *shadernames[2] = {SHADERS_PATH "/Pathtracer/pathtracer_gbuffer_noise.vert",
                                  SHADERS_PATH "/Pathtracer/pathtracer_gbuffer_noise.frag"};
    m_shader_gbuffer = std::unique_ptr<CVK::ShaderMinimal>(
            new CVK::ShaderMinimal(VERTEX_SHADER_BIT | FRAGMENT_SHADER_BIT, shadernames));

    initPlane(512, 512);
}

KIRK::Noise::Noise(int width, int height)
{
    m_fbo = std::unique_ptr<CVK::FBO>(new CVK::FBO(width, height, 4, true, false));

    const char *shadernames[2] = {SHADERS_PATH "/Pathtracer/pathtracer_gbuffer_noise.vert",
                                  SHADERS_PATH "/Pathtracer/pathtracer_gbuffer_noise.frag"};
    m_shader_gbuffer = std::unique_ptr<CVK::ShaderMinimal>(
            new CVK::ShaderMinimal(VERTEX_SHADER_BIT | FRAGMENT_SHADER_BIT, shadernames));

    initPlane(width, height);

}


KIRK::Noise::~Noise()
{

}

void KIRK::Noise::generate_noise()
{
    CVK::State::getInstance()->getCamera()->setPosition(glm::vec3(0.0, 0.0, -49.f));

    m_fbo->bind();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    CVK::State::getInstance()->setShader(m_shader_gbuffer.get());

    m_shader_gbuffer->update();

    m_plane_node->render();

    m_fbo->unbind();

    CVK::State::getInstance()->getCamera()->setPosition(glm::vec3(0.0, 0.0, 0.0f));
}

GLuint KIRK::Noise::get_noise_texture()
{
    return m_fbo->getColorTexture(3);
}


void KIRK::Noise::initPlane(int width, int height)
{
    m_plane_node = std::unique_ptr<CVK::Node>(new CVK::Node("PlaneNoise"));
    m_plane = std::make_shared<CVK::Plane>();

    m_plane_node->setGeometry(m_plane);
    m_plane_node->setModelMatrix(glm::scale(*m_plane_node->getModelMatrix(), glm::vec3(width, height, 1.0f)));
    m_plane_node->setModelMatrix(glm::rotate(*m_plane_node->getModelMatrix(), 55.f, glm::vec3(1.0f, 0.0f, 0.f)));

}

