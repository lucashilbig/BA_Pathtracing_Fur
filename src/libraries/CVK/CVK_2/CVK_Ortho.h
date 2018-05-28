#ifndef __CVK_ORTHO_H
#define __CVK_ORTHO_H

#include "CVK_Defs.h"
#include "CVK_Projection.h"

namespace CVK {
/**
* A Projection is orthographic when the size on the screen of an object is independent to its distance.
* @brief Class for orthographic projection matrix
*/
class Ortho : public CVK::Projection
{
public:
    /**
    * Constructor for Ortho with given parameters and standard parameters (z in [0.0, 10.0])
    * @param ratio The ratio of the projection
    */
    Ortho(float ratio);
    /**
     * Constructor for Ortho with given parameters
     * @param left The left value of the projection cuboid (-x)
     * @param right The right value of the projection cuboid (+x)
     * @param bottom The bottom value of the projection cuboid (-y)
     * @param top The top value of the projection cuboid (+y)
     * @param near The near value of the projection cuboid (-z)
     * @param far The far value of the projection cuboid (+z)
     */
    Ortho(float left, float right, float bottom, float top, float near, float far);
    /**
     * Standard Destructor for Ortho
     */
    ~Ortho();

    /**
    * Sets the orthographic projection to a cuboid with given values.
    * @brief Sets orthographic values
    * @param left The left value of the projection cuboid (-x)
    * @param right The right value of the projection cuboid (+x)
    * @param bottom The bottom value of the projection cuboid (-y)
    * @param top The top value of the projection cuboid (+y)
    * @param near The near value of the projection cuboid (-z)
    * @param far The far value of the projection cuboid (+z)
    */
    void setOrtho(float left, float right, float bottom, float top, float near, float far);
    /**
    * @brief Setter for near and far attributes and rebuilds an orthographic projection
    * @param near The nearest depth visible
    * @param far The farthest depth visible
    */
    void setNearFar(float near, float far);
    /**
    * Updates the ratio of the projection specific to an orthographic projection
    * @param ratio The ratio of the projection
    */
    void updateRatio(float ratio);

protected:
    float m_left, m_right, m_bottom, m_top;
};
};

#endif /* __CVK_ORTHO_H */
