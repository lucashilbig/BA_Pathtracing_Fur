#ifndef __CVK_SHADERPHONG_H
#define __CVK_SHADERPHONG_H

#include "CVK_Defs.h"
#include "CVK_Camera.h"
#include "CVK_Light.h"
#include "CVK_Material.h"
#include "CVK_Geometry.h"
#include "CVK_ShaderMinimal.h"

namespace CVK {

/**
* Phong shader class impementation using the ShaderSet and ShaderMinimal. The model, view and projection
* matrices are set using ShaderMinimal. Sets additional variables (f.e. light and fog informations).
* Uses values set in State. Has to collect uniform locations for variables from shader first.
* @brief Phong shader that sets light and fog informations
* @see State
*/
class ShaderPhong : public CVK::ShaderMinimal
{
public:
    /**
    * Constructor for ShaderPhong with given parameters. Collects uniform locations for
    * all used variables from Shader Program.
    * @param shader_mask Describes which shader files are used
    * @param shaderPaths Array of paths to shader files
    */
    ShaderPhong(GLuint shader_mask, const char **shaderPaths);

    /**
    * Sets scene dependent variables in Shader. Namely light and fog informations set in State.
    * @brief Sets scene variables
    * @see State
    * @see Light
    */
    void update();
    /**
    * Sets node dependent variables in Shader, like Material information.
    * @brief Sets node variables
    * @see Material
    */
    void update(CVK::Node *node);

private:
    GLuint m_kdID, m_ksID, m_ktID;
    GLuint m_diffuseID, m_specularID, m_shininessID;
    GLuint m_lightambID;

    GLuint m_numLightsID;
    GLuint m_lightposID[MAX_LIGHTS], m_lightcolID[MAX_LIGHTS], m_lightsdirID[MAX_LIGHTS], m_lightsexpID[MAX_LIGHTS], m_lightscutID[MAX_LIGHTS];

    GLuint m_useColorTexture, m_colorTextureID;

    GLuint m_fogcolID, m_fogstartID, m_fogendID, m_fogdensID, m_fogmodeID;
};

};

#endif /*__CVK_SHADERPHONG_H*/
