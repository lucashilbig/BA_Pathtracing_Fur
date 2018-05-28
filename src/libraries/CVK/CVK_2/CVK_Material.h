#ifndef __CVK_MATERIAL_H
#define __CVK_MATERIAL_H

#include <string>

#include "CVK_Defs.h"
#include "CVK_Texture.h"

namespace CVK {
/**
 * Enum for different texture types. New types might follow.
 * @brief Different texture types
 */
enum TextureType
{
    COLOR_TEXTURE,
    NORMAL_TEXTURE
};

/**
 * Material are used for lighting within a shader. For this the objects of this class store
 * the needed values. Furthermore they store all necessary OpenGL texture objects for an
 * object of the scene.
 * @brief Used for storing information for lighting
 */
class Material
{
public:
    /**
     * Constructor for Material with given parameters
     * @param diffuse The color for direct illumination only dependent on angle between normal and light position
     * @param specular The color for simple reflective illumination dependent on light, normal and camera
     * @param shininess The exponent for the angle for reflective specular illumination
     */
    Material(glm::vec3 diffuse, glm::vec3 specular, float shininess, float ior);
    /**
    * Constructor for Material with given parameters
    * @param kd The factor of diffuse illumination
    * @param diffuse The color for direct illumination only dependent on angle between normal and light position
    */
    Material(float kd, glm::vec3 diffuse, float ior);
    /**
    * Constructor for Material with given parameters
    * @param kd The factor of diffuse illumination
    * @param diffuse The color for direct illumination only dependent on angle between normal and light position
    * @param ks The factor of specular illumination
    * @param specular The color for simple reflective illumination dependent on light, normal and camera
    * @param shininess The exponent for the angle for reflective specular illumination
    */

    Material(float kd, glm::vec3 diffuse, float ks, glm::vec3 specular, float shininess, float ior = 1.0f,
             float alpha = 1.0f, int illum = 0);
    Material(glm::vec3 kd, glm::vec3 diffuse, glm::vec3 ks, glm::vec3 specular, float shininess, float ior = 1.0f,
             float alpha = 1.0f, int illum = 0);
    /**
    * Constructor for Material with given parameters
    * @param colorTexturePath The path to the color image for this material, used instead of diffuse color. Image is loaded first.
    * @param kd The factor of diffuse illumination
    */
    Material(const char *colorTexturePath, float kd, float ior);
    /**
    * Constructor for Material with given parameters
    * @param colorTextureID The OpenGL id for the already created OpenGL texture object
    * @param kd The factor of diffuse illumination
    */
    Material(GLuint colorTextureID, float kd, float ior);
    /**
    * Constructor for Material with given parameters
    * @param colorTexturePath The path to the color image for this material, used instead of diffuse color. Image is loaded first.
    * @param kd The factor of diffuse illumination
    * @param ks The factor of specular illumination
    * @param specular The color for simple reflective illumination dependent on light, normal and camera
    * @param shininess The exponent for the angle for reflective specular illumination
    */
    Material(const char *colorTexturePath, float kd, float ks, glm::vec3 specular, float shininess, float ior = 1.0f,
             float alpha = 1.0f, int illum = 0);
    /**
    * Constructor for Material with given parameters
    * @param colorTextureID The OpenGL id for the already created OpenGL texture object
    * @param kd The factor of diffuse illumination
    * @param ks The factor of specular illumination
    * @param specular The color for simple reflective illumination dependent on light, normal and camera
    * @param shininess The exponent for the angle for reflective specular illumination
    */
    Material(GLuint colorTextureID, float kd, float ks, glm::vec3 specular, float shininess, float ior);
    /**
     * Standard Destructor for Material
     */
    ~Material();

    /**
    * This function is called by the constructors. Normally it is not necessary to call it manually.
    * @brief Initializes Phong Material parameters
    * @param kd The factor of diffuse illumination
    * @param diffuse The color for direct illumination only dependent on angle between normal and light position
    * @param ks The factor of specular illumination
    * @param specular The color for simple reflective illumination dependent on light, normal and camera
    * @param shininess The exponent for the angle for reflective specular illumination
    */
    void init(float kd, glm::vec3 diffuse, float ks, glm::vec3 specular, float shininess, float ior = 1.0f,
              float alpha = 1.0f, int illum = 0);
    void init(glm::vec3 kd, glm::vec3 diffuse, glm::vec3 ks, glm::vec3 specular, float shininess, float ior = 1.0f,
              float alpha = 1.0f, int illum = 0);
    /**
     * @brief Standard Setter for factor of diffuse illumination
     * @param kd the new factor of diffuse illumination of this object
     */
    void setKd(float kd);
    /**
     * @brief Standard Getter for factor of diffuse illumination
     * @return the factor of diffuse illumination of this object
     */
    float getKd();
    glm::vec3 getKdVec();
    glm::vec3 *getKdVecP();
    /**
    * @brief Standard Setter for factor of specular illumination
    * @param ks the new factor of specular illumination of this object
    */
    void setKs(float ks);
    /**
    * @brief Standard Getter for factor of specular illumination
    * @return the factor of specular illumination of this object
    */
    float getKs();
    glm::vec3 getKsVec();
    glm::vec3 *getKsVecP();
    /**
    * @brief Standard Setter for factor of transparency
    * @param kt the new factor of transparency of this object
    */
    void setKt(float kt);
    /**
    * @brief Standard Getter for factor of transparency
    * @return the factor of transparency of this object
    */
    float getKt();
    glm::vec3 getKtVec();
    glm::vec3 *getKtVecP(); //TODO

    /**
     * @brief Standard Setter for diffuse color
     * @param col the new diffuse color of this object
     */
    void setdiffColor(glm::vec3 col);
    /**
     * @brief Standard Getter for diffuse color
     * @return the diffuse color of this object
     */
    glm::vec3 *getdiffColor();

    /**
    * @brief Standard Setter for specular color
    * @param col the new specular color of this object
    */
    void setspecColor(glm::vec3 col);
    /**
    * @brief Standard Getter for specular color
    * @return the specular color of this object
    */
    glm::vec3 *getspecColor();
    /**
    * @brief Standard Setter for emitted color
    * @param col the new emitted color of this object
    */
    void setEmittedColor(const glm::vec3 &col);
    /**
    * @brief Standard Getter for emitted color
    * @return the emitted color of this object
    */
    const glm::vec3 &getEmittedColor() const;
    /**
     * @brief Standard Setter for shininess
     * @param shininess the new shininess of this object
     */
    void setShininess(float shininess);
    /**
     * @brief Standard Getter for shininess
     * @return the shininess of this object
     */
    float getShininess();
    /**
     * @brief Standard Setter for factor of refraction
     * @param ior the new factor of refraction of this object
     */
    void setIor(float ior);
    /**
     * @brief Standard Getter for factor of refraction
     * @return the factor of refraction of this object
     */
    float getIor();

    /**
     * @brief Standard Setter for illumination model
     * @param illum the new illumination model for this object
     */
    void setIllum(int illum);
    /**
     * @brief Standard Getter ffor illumination model
     * @return the illumnation model
     */
    int getIllum();

    /**
     * @brief Standard Setter for factor of transparancy
     * @param alpha the new transparency value of this object
     */
    void setAlpha(float alpha);
    /**
     * @brief Standard Getter for factor of transparancy
     * @return the transparency of this object
     */
    float getAlpha();

    /**
     * @brief Setter for texture of given texture type. Image is loaded first.
     * @param type the type of the texture
     * @param fileName the type of the texture
     */
    void setTexture(TextureType type, const char *fileName);
    /**
    * @brief Setter for texture of given texture type.
    * @param type the type of the texture
    * @param textureID the OpenGL ID of the OpenGL texture object
    */
    void setTexture(TextureType type, GLuint textureID);
    /**
     * @brief Getter for determining, if texture of given type is set
     * @return true, if such a texture is set, false otherwise
     */
    bool hasTexture(TextureType type);
    /**
     * @brief Standard Getter for texture of given type
     * @param type the type of the texture
     * @return The texture of given type of this object
     */
    CVK::Texture *getTexture(TextureType type);

    /**
     * @brief Sets the index for the Material in the GPU Material Buffer
     * @param index
     */
    void setMaterialsIndex(int index) { materials_index = index; }

    /**
     * @brief Return the index of the Material in the GPU Material Buffer
     * @return index
     */
    int getMaterialsIndex() { return materials_index; }

    /**
     * @brief Standard Getter for id
     * @return the id
     */
    int getId() { return m_id; }

    inline std::string getColorTexturePath() { return m_colorTexturePath; }

    inline std::string getNormalTexturePath() { return m_normalTexturePath; }


private:
    // good reference page for mtl files
    // http://www.fileformat.info/format/material/
    // in the MTL file

    int m_id;
    int materials_index;
    static int count;

    float m_kd, m_ks, m_kt;                                                        // Kd Ks Kt
    glm::vec3 m_kdVec;                                                            // Kd just as a vec if they differ
    glm::vec3 m_ksVec;                                                            // Ks just as a vec if they differ
    glm::vec3 m_diffColor; //!< diffuse color of material
    glm::vec3 m_specColor; //!< specular color of material
    glm::vec3 m_emitColor; //!< emitted color of material
    float m_shininess; //!< shininess exponent for specular illumination 		// Ns
    float m_ior; //!< index of refraction										// Ni
    int m_illum; // illumination model as written in the comments at Node.cpp	// illum
    float m_alpha; // transparencie of the object								// d or Tr = 1 - d

    std::string m_colorTexturePath;
    std::string m_normalTexturePath;

    CVK::Texture *m_colorTexture; //!< color texture for this material
    CVK::Texture *m_normalTexture; //!< normal texture for this material
};

};

#endif /* __CVK_MATERIAL_H */

