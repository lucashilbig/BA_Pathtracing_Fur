#ifndef __CVK_CAMERA_H
#define __CVK_CAMERA_H

#include "CVK_Defs.h"
#include "CVK_Projection.h"

namespace CVK {
/**
 * Abstract Class for storing the Camera values with movement. It represents
 * a camera in 3D space and therefore contains a pointer to the corresponding
 * projection matrix.
 * @brief Class for Camera with movement
 */
class Camera
{
public:
    /**
     * Standard constructor
     */
    Camera();
    /**
    * Standard destructor
    */
	virtual ~Camera();

    /**
     * Virtual function which needs to be implemented in the subclasses. Is used for all movements.
     * @brief Moves the Camera.
     * @param window The window of the application. Is used to get the mouse movements.
     */
    virtual void update(GLFWwindow *window) = 0;

    /**
     * Returns the view matrix corresponding to this Camera as pointer.
     * @brief Getter for the view matrix
     * @return a pointer to the view matrix
     */
    glm::mat4 getView();
    /**
    * Returns the view matrix corresponding to this Camera as pointer arguments.
    * @brief Getter for the view matrix
    * @param x pointer to the x Axis
    * @param y pointer to the y Axis
    * @param z pointer to the z Axis
    * @param pos pointer to the position of the camera
    */
    void getView(glm::vec3 *x, glm::vec3 *y, glm::vec3 *z, glm::vec3 *pos);
    /**
    * Sets the view matrix corresponding to this Camera.
    * @brief Setter for the view matrix
    * @param view the new view matrix as pointer
    */
    void setView(glm::mat4 view);

    /**
     * deprecated.
     * Used in Raytracer. TODO: Delete... Put it in the Projection class
     */
	void setWidthHeight(int width, int height);
    /**
    * deprecated.
    * Used in Raytracer. TODO: Delete... Put it in the Projection class
    */
    void getWidthHeight(int *width, int *height);

    /**
     * Standard lookAt Function for the camera.
     * @brief Moves Camera and look at the center
     * @param position the new Position for the camera.
     * @param center the center point to look at.
     * @param up the tilt of the camera.
     */
    void lookAt(glm::vec3 position, glm::vec3 center, glm::vec3 up);

    /**
     * Sets the projection matrix.
     * @brief Setter for the projection matrix
     * @param projection the new projection matrix as pointer
     */
    void setProjection(CVK::Projection *projection);

    /**
    * Returns the projection matrix.
    * @brief Getter for the projection matrix
    * @return the projection matrix as pointer
    */
    CVK::Projection *getProjection();

    /**
    * @brief Standard Setter for position
    * @param pos The new position of this object
    */
    virtual void setPosition(glm::vec3 pos) = 0;

    /**
    * Returns the focus distance
    * @brief Getter for the focus distance
    * @return the focus distance value
    */
    float getFocus();

    /**
    * Sets the focus distance
    * @brief Setter for the focus distance
    * @param input the focus distance value
    */
    void setFocus(float input);

    /**
    * Returns the blur radius
    * @brief Getter for the blur radius
    * @return the blur radius value
    */
    float getBlurRadius();

    /**
    * Sets the blur radius
    * @brief Setter for the blur radius
    * @param input new blur radius value
    */
    void setBlurRadius(float input);

    glm::vec3  m_cameraPos, m_direction, m_up; //!< camera attributes

protected:
    int m_width, m_height; //!< deprecated width and height for the projection. TODO remove...
    glm::mat4 m_viewmatrix; //!< the view matrix for this object
    CVK::Projection *m_projection; //!< a pointer to the corresponding projection matrix.
    float m_oldX, m_oldY; //!< Old values for movement

    float m_focusDistance;        //!< camera focus distance
    float m_blurRadius;            //!< camera blur radius
};
};

#endif /* __CVK_CAMERA_H */
