#ifndef __CVK_FREEFLIGHT_H
#define __CVK_FREEFLIGHT_H

#include <chrono>
#include "CVK_Camera.h"
#include "CVK_Perspective.h"

namespace CVK {

/**
* Implementation of the Camera. Movement is done with constant radius around a given center point.
* @brief Camera movement around center
*/
class FreeFlight : public CVK::Camera
{
public:
	static constexpr float epsilon = 0.99f;
    /**
    * Constructor for Trackball with given parameters
    * @param width The width of the camera, used for projection
    * @param height The height of the camera, used for projection
    * @param projection The corresponding projecion matrix
    */
	FreeFlight(int width, int height, CVK::Projection *projection = nullptr);

    /**
     * Standard Destructor for Trackball
     */
    ~FreeFlight();

    /**
    * Update Function to move camera according to the mouse position and the key controls
    * @brief Update position and look
    * @param window The window where OpenGL is running
    */
    void update(GLFWwindow *window) override;

    /**
    * @brief Standard Setter for up vector
    * @param up The new up vector of this object
    */
    void setUpvector(glm::vec3 up);
    void setPosition(glm::vec3 pos) override;

private:
	std::chrono::time_point<std::chrono::steady_clock> m_clock;
	GLFWwindow* m_window;

    float m_sensitivity; //!< the sensitivity of the mouse movement
	float m_flyspeed;
	float m_altitudespeed;
    float m_theta, m_phi;
	glm::vec3 m_dirdefault, m_positiondefault;
};

};

#endif /* __CVK_TRACKBALL_H */
