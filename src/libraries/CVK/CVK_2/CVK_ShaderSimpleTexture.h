#ifndef __SHADER_SIMPLE_TEXTURE_H
#define __SHADER_SIMPLE_TEXTURE_H

#include "CVK_Defs.h"
#include "CVK_ShaderPostProcessing.h"
#include "CVK_State.h"

namespace CVK {

/**
* Shader impementation for Post Processing using the ShaderPostProcessing. Simply sets the first
* texture from ShaderSet and renders screen filling quad.
* @brief Shader setting one texture for post processing
*/
class ShaderSimpleTexture : public CVK::ShaderPostProcessing
{
public:

    /**
    * Constructor for ShaderSimpleTexture with given parameters. Collects uniform locations for
    * all used variables from Shader Program.
    * @param shader_mask Describes which shader files are used
    * @param shaderPaths Array of paths to shader files
    */
    ShaderSimpleTexture(GLuint shader_mask, const char **shaderPaths);

    /**
    * Sets rendering dependent variables in Shader. Namely the first texture from ShaderSet texture list.
    * @brief Sets first texture as rendering variable
    */
    virtual void update();

private:
    GLuint m_colorTextureID;

};

};

#endif /*__SHADER_SIMPLE_TEXTURE_H*/