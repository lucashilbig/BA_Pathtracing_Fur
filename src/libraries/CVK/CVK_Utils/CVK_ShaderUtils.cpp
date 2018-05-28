#include "CVK_ShaderUtils.h"

std::string CVK_UtilsShadersPath = SHADERS_PATH;

void CVK::ShaderUtils::init(std::string pathToShaders)
{
    CVK_UtilsShadersPath = pathToShaders;
}

std::string CVK::ShaderUtils::toShaderPath(std::string path)
{
    return (CVK_UtilsShadersPath + path).c_str();
}


CVK::ShaderUtils::ScreenFillShader::ScreenFillShader()
{
    std::string v = toShaderPath("/Examples/screenFill.vert");
    std::string f = toShaderPath("/Examples/simpleTexture.frag");
    const char *shadernamesSimpleTexture[2] = {v.c_str(), f.c_str()};
    m_shader =
            new ShaderSimpleTexture(VERTEX_SHADER_BIT | FRAGMENT_SHADER_BIT, shadernamesSimpleTexture);
}

void CVK::ShaderUtils::ScreenFillShader::setRenderTexture(unsigned int texture)
{
    if(texture != 0)
        m_shader->setTextureInput(0, texture);
}

void CVK::ShaderUtils::ScreenFillShader::render(unsigned int texture)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    if(texture != 0)
        m_shader->setTextureInput(0, texture);
    m_shader->useProgram();
    m_shader->update();
    glDisable(GL_DEPTH_TEST);
    m_shader->render();
    glEnable(GL_DEPTH_TEST);
}


CVK::ShaderUtils::PhongToScreenShader::PhongToScreenShader(std::shared_ptr<CVK::Node> scene, std::string shaderName)
{
    std::string v = toShaderPath(shaderName + ".vert");
    std::string f = toShaderPath(shaderName + ".frag");
    const char *shadernames[2] = {v.c_str(), f.c_str()};
    m_shader =
            new CVK::ShaderPhong(VERTEX_SHADER_BIT | FRAGMENT_SHADER_BIT, shadernames);
    m_scene = scene;
}

void CVK::ShaderUtils::PhongToScreenShader::render(std::shared_ptr<CVK::Node> scene)
{
    if(scene == 0)
        scene = m_scene;
    if(scene == 0)
        return;
    CVK::State::getInstance()->setShader(m_shader);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    m_shader->useProgram();
    m_shader->update();
    scene->render();
}


CVK::ShaderUtils::GBufferShader::GBufferShader(CVK::FBO *fbo, CVK::Node *scene, std::string shader_name)
{
    std::string v = toShaderPath(shader_name + ".vert");
    std::string f = toShaderPath(shader_name + ".frag");
    const char *shadernames[2] = {v.c_str(), f.c_str()};

    setShaderPaths(VERTEX_SHADER_BIT | FRAGMENT_SHADER_BIT, shadernames);

    m_scene = scene;
    m_fbo = fbo;
}

void CVK::ShaderUtils::GBufferShader::render(CVK::Node *scene)
{
    if(scene != NULL)
        m_scene = scene;

    m_fbo->bind();
    glClearColor(0.f, 0.f, 0.f, 0.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    CVK::State::getInstance()->setShader(this);

    this->useProgram();
    this->update();

    m_scene->render();

    m_fbo->unbind();
}

void CVK::ShaderUtils::GBufferShader::update(CVK::Node *node)
{
    if(node->hasMaterial())
    {
        auto mat = node->getMaterial();
        CVK::Texture *color_texture;

        setUniform("mat.kd", mat->getKdVecP()->x);
        setUniform("mat.ks", mat->getKsVecP()->x);
        setUniform("mat.kt", mat->getKt());
        setUniform("mat.diffColor", *mat->getdiffColor());
        setUniform("mat.specColor", *mat->getspecColor());
        setUniform("mat.shininess", mat->getShininess());

        bool has_col_tex = mat->hasTexture(COLOR_TEXTURE);
        setUniform("useColorTexture", has_col_tex);

        if(has_col_tex)
        {
            setUniform("colorTexture", 0);

            glActiveTexture(COLOR_TEXTURE_UNIT);
            color_texture = mat->getTexture(COLOR_TEXTURE);
            color_texture->bind();
        }
    }
}
