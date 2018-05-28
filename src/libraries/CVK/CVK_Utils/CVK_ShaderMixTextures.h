#ifndef __SHADERMIXTEXTURES_H__
#define __SHADERMIXTEXTURES_H__

#include "CVK/CVK_2/CVK_ShaderPostProcessing.h"

namespace CVK {

class ShaderMixTextures : public CVK::ShaderPostProcessing
{
public:

    ShaderMixTextures(GLuint shader_mask, const char **shader_paths);
    void update();

private:
    GLuint m_colorTextureID1, m_colorTextureID2;
};
}

#endif

