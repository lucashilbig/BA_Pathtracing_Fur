#ifndef __CVK_TEXTURE_H
#define __CVK_TEXTURE_H

#include "CVK_Defs.h"

namespace CVK {

/**
 * Class for abstraction of OpenGL Texture Objects. Image objects are loaded and
 * OpenGL texture objects are created and stored.
 * @brief Class for creating and storing OpenGL Texture Objects
 */
class Texture
{
public:
    /**
     * Constructor for Texture with given parameters. Loads Image File first.
     * @param fileName The path to a file which needs to be loaded
     */
    Texture(const char *fileName);
    /**
    * Constructor for Texture with given parameters
    * @param width The width of the texture
    * @param height The height of the texture
    * @param bytesPerPixel The number of channels per Pixel
    * @param data The data of an already loaded Image File
    */
    Texture(int width, int height, int bytesPerPixel, unsigned char *data);
    /**
    * Constructor for Texture with given parameters
    * @param texture An already created OpenGL Texture Objects
    */
    Texture(GLuint texture);
    /**
     * Standard Destructor for Texture
     */
    ~Texture();

    /**
     * Loads an image file at the given path and creates an OpenGL texture object with it
     * @brief Loads image and creates texture
     * @param fileName The path to a file which needs to be loaded
     * @return false, if an error occured, true otherwise
     */
    bool load(const char *fileName);
    /**
     * Binds the OpenGL texture object to the current GL_TEXTURE_2D, so that it can be used in the shader
     * @brief Binds texture for usage
     */
    void bind();

    /**
    * @brief Setter for texture
    * @param texture An already created OpenGL Texture Objects
    */
    void setTexture(GLuint texture);
    /**
     * @brief Setter for texture
     * @param width The width of the texture
     * @param height The height of the texture
     * @param bytesPerPixel The number of channels per Pixel
     * @param data The data of an already loaded Image File
     */
    void setTexture(int width, int height, int bytesPerPixel, unsigned char *data);
    /**
     * @brief Standard Getter for OpenGL texture object
     * @return The OpenGL texture object of this object
     */
    unsigned int getTexture();
    /**
     * Returns the color of the data at the given texture coordinate
     * @brief Returns color at texture coordinate
     * @param tcoord The texture coordinate
     * @return the color at the given texture coordinate
     */
    glm::vec3 getValue(glm::vec2 tcoord);

private:
    void createTexture();

    unsigned int m_textureID;
    int m_width, m_height;

    unsigned char *m_data; //keep copy e.g. for CPU ray tracing
    int m_bytesPerPixel;
};

};

#endif /*__CVK_TEXTURE_H*/
