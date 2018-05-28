
#ifndef __SHADER_UTILS_H__
#define __SHADER_UTILS_H__

#include <memory>

#include <CVK/CVK_2/CVK_Framework.h>


namespace CVK {
namespace ShaderUtils {
void init(std::string pathToShaders = SHADERS_PATH);
std::string toShaderPath(std::string path);

class ScreenFillShader
{
public:
    ScreenFillShader();

    virtual ~ScreenFillShader() {}

    void setRenderTexture(unsigned int texture);
    void render(unsigned int texture = 0);
private:
    CVK::ShaderSimpleTexture *m_shader;
};

class PhongToScreenShader
{
public:
    // shaderName f.e. "/Phong/Phong"
    PhongToScreenShader(std::shared_ptr<CVK::Node> scene = 0, std::string shaderName = "/Phong/Phong");

    virtual ~PhongToScreenShader() {}

    void setRenderScene(std::shared_ptr<CVK::Node> scene) { m_scene = scene; }

    void render(std::shared_ptr<CVK::Node> scene = 0);

    CVK::ShaderPhong *getShader() { return m_shader; }

private:
    CVK::ShaderPhong *m_shader;
	std::shared_ptr<CVK::Node> m_scene;
};

class GBufferShader : public CVK::ShaderMinimal
{
public:
    GBufferShader() {}

    GBufferShader(CVK::FBO *fbo = NULL, CVK::Node *scene = NULL, std::string shader_name = "/Minimal/simple_gbuffer");

    virtual ~GBufferShader() {}

    /** scene to render */
    void setRenderScene(CVK::Node *scene) { m_scene = scene; }

    /** FBO render target */
    void setRenderFBO(CVK::FBO *fbo) { m_fbo = fbo; }

    void render(CVK::Node *scene = NULL);

    void update() { CVK::ShaderMinimal::update(); }

    void update(CVK::Node *node);
private:
    CVK::Node *m_scene;
    CVK::FBO *m_fbo;
};

};

};


#endif //__SHADER_UTILS_H__
