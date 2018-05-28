#ifndef __CVK_FBO_H
#define __CVK_FBO_H

#include <cstdio>

#include "CVK_Defs.h"
#include "CVK_State.h"

namespace CVK {

/**
 * @brief Class for OpenGL frame buffer objects (FBOs).
 *
 * A FBO has a given width and height. All textures used by it
 * or created by it have these values for width and height. It
 * is possible to create a fbo with multiple color buffers, each
 * with rgba values per pixel. A fbo can use a depth buffer and
 * a stencil buffer, which are not necesserily needed.
 *
 * To create a fbo, the constructor has to be called with the
 * preferred width and height. Additionally the number of color
 * textures and the use of depth or stencil buffer can be assigned.
 * An example creation of a fbo with 1 color buffer and a depth
 * buffer is:
 * CVK::FBO fbo( 400, 300, 1, true);
 *
 * To use a fbo, it has to be bound to OpenGL. If the standard
 * OpenGL frame buffer should be used, the fbo has to be unbound.
 * fbo.bind();
 * ...
 * fbo.unbind();
 *
 * To use the result of the fbo, getters for the colorTextures
 * and the depth texture exist.
 * GLuint colorTexture = fbo.getColorTexture(0);
 * GLuint depthTexture = fbo.getDepthTexture();
 */
class FBO
{
public:
    /**
    * Constructor for FBO with given parameters
    * @param width The width for each texture of the FBO
    * @param height The height for each texture of the FBO
    * @param numColorTextures The number of color textures for this FBO
    * @param depthTexture Declares if the FBO has a depth texture for depth comparison
    * @param stencilTexture Declares if the FBO has a stencil texture for stencil comparison
    */
    FBO(int width, int height, int numColorTextures = 1, bool depthTexture = false, bool stencilTexture = false);
    /**
     * Standard Destructor for FBO
     */
    ~FBO();

    // called by constructors
    /**
     * Creates a new FBO with given parameters
     * @brief Creates FBO
     * @param width The width for each texture of the FBO
     * @param height The height for each texture of the FBO
     * @param numColorTextures The number of color textures for this FBO
     * @param depthTexture Declares if the FBO has a depth texture for depth comparison
     * @param stencilTexture Declares if the FBO has a stencil texture for stencil comparison
     */
    void create(int width, int height, int numColorTextures, bool depthTexture, bool stencilTexture);
    /**
     * Resets FBO and all corresponding textures and sets them to the invalid value
     * @brief Resets FBO and all textures
     */
    void reset();
    /**
     * Resizes all textures in this FBO to the given parameters
     * @brief Resizes FBO textures
     * @param width the new width of all textures
     * @param height the new height of all textures
     */
    void resize(int width, int height);

    /**
     * Binds the FBO, so that all OpenGL draw commands execute in this Frame Buffer.
     * @brief Bind FBO for usage with OpenGL
     */
    void bind();
    /**
    * Unbinds the FBO, so that all OpenGL draw commands execute in the default Frame Buffer.
    * @brief Unbind FBO
    */
    void unbind();

    /**
     * @brief Standard Getter for n-th color texture object
     * @param index the index of the color texture object of this FBO
     * @return the color texture object at given index
     */
    GLuint getColorTexture(unsigned int index);
    /**
     * @brief Standard Getter for depth texture object
     * @return the depth texture object of this FBO
     */
    GLuint getDepthTexture();

private:
    GLuint createTexture();

    int m_width, m_height; //!< the size of each texture object
    int m_screen_x, m_screen_y; //!< the screen resolution (temp for unbind() )
    std::vector<GLuint> m_colorTextures; //!< List containing all OpenGL color textures
    GLuint m_depthTexture, m_stencilTexture;
    GLuint m_frameBufferHandle; //!< The actual OpenGL FBO id
};

};

#endif /* __CVK_FBO_H */

