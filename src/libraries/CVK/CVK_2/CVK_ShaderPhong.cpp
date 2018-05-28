#include "CVK_ShaderPhong.h"


CVK::ShaderPhong::ShaderPhong(GLuint shader_mask, const char **shaderPaths) : CVK::ShaderMinimal(shader_mask,
                                                                                                 shaderPaths)
{
    //Material
    m_kdID = glGetUniformLocation(m_ProgramID, "mat.kd");
    m_ksID = glGetUniformLocation(m_ProgramID, "mat.ks");
    m_ktID = glGetUniformLocation(m_ProgramID, "mat.kt");

    m_diffuseID = glGetUniformLocation(m_ProgramID, "mat.diffColor");
    m_specularID = glGetUniformLocation(m_ProgramID, "mat.specColor");
    m_shininessID = glGetUniformLocation(m_ProgramID, "mat.shininess");

    //Light
    m_numLightsID = glGetUniformLocation(m_ProgramID, "numLights");
    for(unsigned int i = 0; i < MAX_LIGHTS; i++)
    {
        char string[100];
        sprintf(string, "light[%d].pos", i);
        m_lightposID[i] = glGetUniformLocation(m_ProgramID, string);
        sprintf(string, "light[%d].col", i);
        m_lightcolID[i] = glGetUniformLocation(m_ProgramID, string);
        sprintf(string, "light[%d].spot_direction", i);
        m_lightsdirID[i] = glGetUniformLocation(m_ProgramID, string);
        sprintf(string, "light[%d].spot_exponent", i);
        m_lightsexpID[i] = glGetUniformLocation(m_ProgramID, string);
        sprintf(string, "light[%d].spot_cutoff", i);
        m_lightscutID[i] = glGetUniformLocation(m_ProgramID, string);
    }

    //Textures
    m_useColorTexture = glGetUniformLocation(m_ProgramID, "useColorTexture");
    m_colorTextureID = glGetUniformLocation(m_ProgramID, "colorTexture");

    //Fog
    m_lightambID = glGetUniformLocation(m_ProgramID, "lightAmbient");
    m_fogcolID = glGetUniformLocation(m_ProgramID, "fog.col");
    m_fogstartID = glGetUniformLocation(m_ProgramID, "fog.start");
    m_fogendID = glGetUniformLocation(m_ProgramID, "fog.end");
    m_fogdensID = glGetUniformLocation(m_ProgramID, "fog.density");
    m_fogmodeID = glGetUniformLocation(m_ProgramID, "fog.mode");
}

void CVK::ShaderPhong::update()
{
    int numLights = CVK::State::getInstance()->getLights()->size();
    CVK::ShaderMinimal::update();

    glUniform1i(m_numLightsID, numLights);
    for(int i = 0; i < numLights; i++)
    {
        CVK::Light *light = &CVK::State::getInstance()->getLights()->at(i);
        glUniform4fv(m_lightposID[i], 1, glm::value_ptr(*light->getPosition()));
        glUniform3fv(m_lightcolID[i], 1, glm::value_ptr(*light->getColor()));
        glUniform3fv(m_lightsdirID[i], 1, glm::value_ptr(*light->getSpotDirection()));
        glUniform1f(m_lightsexpID[i], light->getSpotExponent());
        glUniform1f(m_lightscutID[i], light->getSpotCutoff());
    }

    glUniform3fv(m_lightambID, 1, glm::value_ptr(CVK::State::getInstance()->getLightAmbient()));

    glUniform3fv(m_fogcolID, 1, glm::value_ptr(CVK::State::getInstance()->getFogCol()));
    glUniform1i(m_fogmodeID, CVK::State::getInstance()->getFogMode());
    glUniform1f(m_fogstartID, CVK::State::getInstance()->getFogStart());
    glUniform1f(m_fogendID, CVK::State::getInstance()->getFogEnd());
    glUniform1f(m_fogdensID, CVK::State::getInstance()->getFogDens());
}

void CVK::ShaderPhong::update(CVK::Node *node)
{

    if(node->hasMaterial())
    {
        auto mat = node->getMaterial();
        CVK::Texture *color_texture;

        glUniform3fv(m_kdID, 1, glm::value_ptr(*mat->getKdVecP()));
        glUniform3fv(m_ksID, 1, glm::value_ptr(*mat->getKsVecP()));
        glUniform1f(m_ktID, mat->getKt());
        glUniform3fv(m_diffuseID, 1, glm::value_ptr(*mat->getdiffColor()));
        glUniform3fv(m_specularID, 1, glm::value_ptr(*mat->getspecColor()));
        glUniform1f(m_shininessID, mat->getShininess());

        bool colorTexture = mat->hasTexture(COLOR_TEXTURE);
        glUniform1i(m_useColorTexture, colorTexture);

        if(colorTexture)
        {
            //TODO: COLOR_TEXTURE_UNIT
            glUniform1i(m_colorTextureID, 0);

            glActiveTexture(COLOR_TEXTURE_UNIT);
            color_texture = mat->getTexture(COLOR_TEXTURE);
            color_texture->bind();
        }
    }
}
