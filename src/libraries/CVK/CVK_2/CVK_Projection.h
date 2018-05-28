#ifndef __CVK_PROJECTION_H
#define __CVK_PROJECTION_H

#include "CVK_Defs.h"

namespace CVK {
/**
 * Abstract Class to store and change the projection matrix.
 * @brief Abstract Class for projection matrix.
 */
class Projection
{
public:
    /**
     * @brief Standard Getter for projection matrix
     * @return The projection matrix of this object as pointer
     */
    glm::mat4 getProjMatrix();
    /**
     * @brief Standard Setter for projection matrix
     * @param projection The new projection matrix of this object as pointer
     */
    void setProjMatrix(glm::mat4 projection);

    /**
     * @brief Getter for near and far value
     * @param near The near value of this object as pointer
     * @param far The far value of this object as pointer
     */
    void getNearFar(float *near, float *far);
    /**
     * @brief Standard Getter for near value
     * @return The near value of this object
     */
    float getNear();
    /**
     * Updates the Ratio of the projection matrix. Dependant on subclass implementation and therefore abstract
     * @brief Update needs implementation in subclass
     */
    virtual void updateRatio(float ratio) = 0;

protected:
    float m_znear, m_zfar; //!< near and far value of the projection matrix

    glm::mat4 m_projection; //!< projection matrix
};
};

#endif /* __CVK_PROJECTION_H */
