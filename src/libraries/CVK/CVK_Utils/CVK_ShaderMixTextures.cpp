#include "CVK_ShaderMixTextures.h"

CVK::ShaderMixTextures::ShaderMixTextures(GLuint shader_mask, const char **shader_paths) :
        CVK::ShaderPostProcessing(shader_mask, shader_paths),
        m_colorTextureID1(glGetUniformLocation(m_ProgramID, "colorTexture1")),
        m_colorTextureID2(glGetUniformLocation(m_ProgramID, "colorTexture2"))
{
}

void CVK::ShaderMixTextures::update()
{
    if(m_textures.size() > 0)
    {
        glUniform1i(m_colorTextureID1, 0);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, m_textures[0]);
    }

    if(m_textures.size() > 1)
    {
        glUniform1i(m_colorTextureID2, 1);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, m_textures[1]);
    }
}
