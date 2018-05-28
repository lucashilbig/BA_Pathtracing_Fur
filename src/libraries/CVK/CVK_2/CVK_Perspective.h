#ifndef __CVK_PERSP_H
#define __CVK_PERSP_H

#include "CVK_Defs.h"
#include "CVK_Projection.h"

namespace CVK //symmetic frustum
{
/**
 * A Projection is perspective when the size on the screen of an object depends on its distance.
 * @brief Class for perspective projection matrix
 */
class Perspective : public CVK::Projection
{
public:
    /**
     * Constructor for Perspective with given parameters and standard parameters (z in [0.001, 10.0])
     * @param ratio The ratio of the projection
     */
    Perspective(float ratio);
    /**
    * Constructor for Perspective with given parameters
    * @param fov The field of view for this projection
    * @param ratio The ratio of the projection
    * @param near The nearest depth visible
    * @param far The farthest depth visible
    */
    Perspective(float fov, float ratio, float near, float far);
    /**
     * Standard Destructor for Perspective
     */
    ~Perspective();

    /**
    * @brief Standard Setter for all Perspective attributes
    * @param fov The field of view for this projection
    * @param ratio The ratio of the projection
    * @param near The nearest depth visible
    * @param far The farthest depth visible
    */
    void setPerspective(float fov, float ratio, float near, float far);

    /**
    * @brief Setter for near and far attributes and rebuilds a perspective projection
    * @param near The nearest depth visible
    * @param far The farthest depth visible
    */
    void setNearFar(float near, float far);
    /**
     * @brief Standard Setter for field of view
     * @param fov The new field of view of this object
     */
    void setFov(float fov);
    /**
     * @brief Standard Getter for field of view
     * @return The field of view of this object
     */
    float getFov();

    /**
     * Updates the ratio of the projection specific to a perspective projection
     * @param ratio The ratio of the projection
     */
    void updateRatio(float ratio);

protected:
    float m_fov, m_ratio;
};
};

#endif /* __CVK_PERSP_H */
