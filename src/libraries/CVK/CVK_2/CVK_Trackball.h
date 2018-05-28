#ifndef __CVK_TRACKBALL_H
#define __CVK_TRACKBALL_H

#include "CVK_Camera.h"
#include "CVK_Perspective.h"

namespace CVK {

/**
* Implementation of the Camera. Movement is done with constant radius around a given center point.
* @brief Camera movement around center
*/
class Trackball : public CVK::Camera
{
public:
    /**
    * Constructor for Trackball with given parameters
    * @param width The width of the camera, used for projection
    * @param height The height of the camera, used for projection
    * @param projection The corresponding projecion matrix
    */
    Trackball(int width, int height, CVK::Projection *projection);
    /**
    * Constructor for Trackball with given parameters
    * @param width The width of the camera, used for projection
    * @param height The height of the camera, used for projection
    */
    Trackball(int width, int height);
    /**
     * Standard Destructor for Trackball
     */
    ~Trackball();

    /**
    * Update Function to move camera according to the mouse position and the key controls
    * @brief Update position and look
    * @param window The window where OpenGL is running
    */
    void update(GLFWwindow *window);
    /**
     * @brief Standard Setter for center point
     * @param center The new center point of this object as pointer
     */
    void setCenter(glm::vec3 *center);
    /**
     * @brief Standard Setter for radius
     * @param radius The new radius of this object
     */
    void setRadius(float radius);
    /**
    * @brief Standard Setter for up vector
    * @param up The new up vector of this object
    */
    void setUpvector(glm::vec3 up);
    void setPosition(glm::vec3 pos);

private:
    float m_sensitivity; //!< the sensitivity of the mouse movement
    float m_theta, m_phi, m_radius;
    glm::vec3 m_center;
};

};

#endif /* __CVK_TRACKBALL_H */
