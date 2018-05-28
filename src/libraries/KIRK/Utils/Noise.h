#ifndef __CVK_RT_NOISE_H
#define __CVK_RT_NOISE_H

#include <memory>
#include "CVK/CVK_Utils/CVK_OpenGLUtils.h"

namespace KIRK {
/*@class: This class provides functions to generate a noise texture*/
class Noise
{
public:
    /*@brief: The constructors which uses a widht and height value. If not given by the user the standard value will be 512x512*/
    Noise();
    Noise(int width, int height);
    ~Noise();

    /*@brief: This method generates the noise textur. It will use a fbo and gbuffer-shader to render noise onto a plane. Beware: The camera will be set to (0, 0 , 0) afterwards.*/
    void generate_noise();

    /* @brief: Returns the texture of the fbo-object which contains the noise texture */
    GLuint get_noise_texture();


private:
    /* @brief: Initialize a plane and the necessary node to render it later */
    void initPlane(int width, int height);

    std::unique_ptr<CVK::FBO> m_fbo;
    std::unique_ptr<CVK::ShaderMinimal> m_shader_gbuffer;

    std::unique_ptr<CVK::Node> m_plane_node;
    std::shared_ptr<CVK::Plane> m_plane;

};

};

#endif /* __CVK_RT_NOISE_H */

