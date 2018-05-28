#ifndef __CVK_STATE_H
#define __CVK_STATE_H

#include "CVK_Defs.h"
#include "CVK_Light.h"
#include "CVK_Camera.h"
#include "CVK_CubeMapTexture.h"

namespace CVK {

// forward declaration
class ShaderMinimal;

/**
 * Class implemented as singleton to store all variables related to rendering
 * in a similar way to global variables
 * @brief Singleton class for the rendering state
 */
class State
{
public:
    /**
     * @brief Getter for instance
     * @return The instance of this class as pointer
     */
    static State *getInstance();

    /**
     * @brief Standard Setter for shader
     * @param shader The new shader of this object as pointer
     */
    void setShader(CVK::ShaderMinimal *shader);
    /**
     * @brief Standard Getter for shader
     * @return The shader of this object as pointer
     */
    CVK::ShaderMinimal *getShader();

    /**
     * @brief Standard Setter for cube map texture
     * @param cubeMap The new cube map texture of this object as pointer
     */
    void setCubeMapTexture(CVK::CubeMapTexture *cubeMap);
    /**
     * @brief Standard Getter for cube map texture
     * @return The cube map texture of this object as pointer
     */
    CVK::CubeMapTexture *getCubeMapTexture();

    /**
     * Adds a Light to the list of all light sources
     * @brief Adds a light
     * @param light The light to add as pointer
     */
    void addLight(CVK::Light *light);
    /**
    * Sets a Light in the list of all light sources
    * @brief Sets a light
    * @param index The index of the light in the list of lights which is to change
    * @param light The light to set as pointer
    */
    void setLight(int index, CVK::Light *light);
    /**
     * Removes the light at the given index from the list of lights
     * @brief Remove light at index
     * @param indexToRemove The index of the light which is to remove
     */
    void removeLight(int indexToRemove);
    /**
     * @brief Standard Getter for list of lights
     * @return The list of lights of this object as pointer
     */
    std::vector<CVK::Light> *getLights();

    /**
     * @brief Standard Setter for camera
     * @param camera The new camera of this object as pointer
     */
    void setCamera(CVK::Camera *camera);
    /**
     * @brief Standard Getter for camera
     * @return The camera of this object as pointer
     */
    CVK::Camera *getCamera();

    /**
    * @brief Standard Setter for glfw window
    * @param window The new glfw window of this object as pointer
    */
    void setWindow(GLFWwindow *window);
    /**
    * @brief Standard Getter for window
    * @return The glfw window of this object as pointer
    */
    GLFWwindow *getWindow();

    /**
     * Setter for given scene Settings
     * @brief Sets all attributes accordingly
     * @param lightAmbient The ambient light color
     * @param fogMode The fog mode to be used by rendering (the shader has to implement the fog)
     * @param fogCol The color of the fog
     * @param fogStart Describes at which point the fog starts
     * @param fogEnd Describes at which point the fog ends
     * @param fogDens Describes strength of fog
     */
    void updateSceneSettings(glm::vec3 lightAmbient, int fogMode, glm::vec3 fogCol, float fogStart, float fogEnd,
                             float fogDens);
	/**
	* @brief Standard Setter for ambient light color
	*/
	void setLightAmbient(glm::vec3 lightAmbient);
    /**
     * @brief Standard Getter for ambient light color
     * @return The ambient light color of this object
     */
    glm::vec3 getLightAmbient();
    /**
     * @brief Standard Setter for background color
     * @param color The new background color of this object
     */
    void setBackgroundColor(glm::vec3 color);
    /**
     * @brief Standard Getter for background color
     * @return The background color of this object
     */
    glm::vec3 getBackgroundColor();
    /**
     * @brief Standard Getter for the specified fog mode
     * @return The specified fog mode of this object
     */
    int getFogMode();
    /**
     * @brief Standard Getter for fog color
     * @return The fog color of this object
     */
    glm::vec3 getFogCol();
    /**
     * @brief Standard Getter for start of the fog
     * @return The start of the fog of this object
     */
    float getFogStart();
    /**
     * @brief Standard Getter for end of the fog
     * @return The end of the fog of this object
     */
    float getFogEnd();
    /**
     * @brief Standard Getter for strenth / density of the fog
     * @return The strenth / density of the fog of this object
     */
    float getFogDens();

private:
    State();
    ~State();

    GLFWwindow *m_window;

    CVK::ShaderMinimal *m_shader;
    CVK::CubeMapTexture *m_cubeMap;
    std::vector<CVK::Light> m_lights;
    CVK::Camera *m_camera;

    glm::vec3 m_lightAmbient;
    glm::vec3 m_BackgroundColor;
    int m_fogMode;
    glm::vec3 m_fogCol;
    float m_fogStart;
    float m_fogEnd;
    float m_fogDens;

    static State *g_instance;
};

};

#endif /* __CVK_STATE_H */

