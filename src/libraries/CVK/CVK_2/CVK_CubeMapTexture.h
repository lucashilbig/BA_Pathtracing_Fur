#ifndef __CVK_CUBE_MAP_TEXTURE_H
#define __CVK_CUBE_MAP_TEXTURE_H

#include "CVK_Defs.h"

namespace CVK {
/**
 * CubeMapTexture is a Texture class for loading environment maps as cube maps.
 * This class loads and stores cube maps and creates the corresponding OpenGL
 * texture objects.
 * @brief Load cube maps and create OpenGL Textures
 */
class CubeMapTexture
{
public:
    /**
     * Constructor for CubeMapTexture with given parameters. Loads the 6 image files.
     * The order for the image objects is: +x, -x, +y, -y, +z, -z.
     * @param fileNames An array containing 6 paths to image files
     * @see load
     */
    CubeMapTexture(const char **fileNames);
    /**
     * Constructor for CubeMapTexture with given parameters
     * @param texture The OpenGL Texture object to use
     */
    CubeMapTexture(GLuint texture);
    /**
     * Standard Destructor for CubeMapTexture
     */
    ~CubeMapTexture();

    /**
     * Loads 6 image objects and creates one OpenGL Cube Map Texture with them.
     * The order for the image objects is: +x, -x, +y, -y, +z, -z.
     * @brief Loads images and creates Texture Object
     * @param fileNames An array containing 6 paths to image files
     * @return true, if everything worked, false otherwise
     */
    bool load(const char **fileNames);
    /**
     * Binds the OpenGL Texture Object, so that GLSL shaders are able to use it.
     * @brief Binds the Texture Object
     */
    void bind();

    /**
     * @brief Standard Setter for texture object
     * @param texture the new texture object of this object
     */
    void setCubeMapTexture(GLuint texture);
    /**
     * @brief Standard Getter for texture object
     * @return the texture object of this object
     */
    unsigned int getCubeMapTexture();

private:
    void createCubeMapTexture();

    unsigned int m_textureID; //!< OpenGL texture object
    int m_width, m_height; //!< size of each image object
};

};

#endif /*__CVK_CUBE_MAP_TEXTURE_H*/
