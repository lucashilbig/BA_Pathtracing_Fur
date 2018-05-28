// Cone along y-axis, origin in base midpoint

#ifndef __CVK_CONE_H
#define __CVK_CONE_H

#include "CVK_Defs.h"
#include "CVK_Geometry.h"

namespace CVK {
/**
 * The cone is one of the example geometry classes in the CVK.
 * @brief Cone class as Geometry
 */
class Cone : public CVK::Geometry
{
public:
    /**
     * Standard Constructor for Cone with 1 unit range
     */
    Cone();
    /**
    * Constructor for Cone using given parameters
    * @param baseradius The radius of the base of the Cone
    * @param apexradius The radius of the top of the Cone
    * @param height The height from base to top
    * @param resolution The resolution of the mantle
    */
    Cone(float baseradius, float apexradius, float height, int resolution);
    /**
    * Constructor for Cone with 1 unit range
    * @param basepoint The position of the base
    * @param apexpoint The position of the top
    * @param baseradius The radius of the base of the Cone
    * @param apexradius The radius of the top of the Cone
    * @param resolution The resolution of the mantle
    */
    Cone(glm::vec3 basepoint, glm::vec3 apexpoint, float baseradius, float apexradius, int resolution);
    /**
     * Standard Destructor for Cone
     */
    ~Cone();

    /**
     * @brief Standard Getter for position of the base
     * @return the position of the base of this object as pointer
     */
    glm::vec3 *getBasepoint();
    /**
     * @brief Standard Getter for position of the top
     * @return the position of the top of this object as pointer
     */
    glm::vec3 *getApexpoint();
    /**
     * @brief Standard Getter for u Axis of the local coordinate system
     * @return the u Axis of the local coordinate system of this object as pointer
     */
    glm::vec3 *get_u();
    /**
    * @brief Standard Getter for v Axis of the local coordinate system
    * @return the v Axis of the local coordinate system of this object as pointer
    */
    glm::vec3 *get_v();
    /**
    * @brief Standard Getter for w Axis of the local coordinate system
    * @return the w Axis of the local coordinate system of this object as pointer
    */
    glm::vec3 *get_w();

    /**
     * @brief Standard Getter for radius at base
     * @return the radius at base of this object
     */
    float getBaseradius();
    /**
     * @brief Standard Getter for radius at top
     * @return the radius at top of this object
     */
    float getApexradius();
    /**
     * @brief Standard Getter for the slope
     * @return the the slope of this object
     */
    float getSlope();

private:
    /**
     * Create the cone and the buffers with the given attributes
     * @brief Create the cone and the buffers
     */
    void create();

    glm::vec3 m_basepoint, m_apexpoint;
    glm::vec3 m_u, m_v, m_w;
    float m_baseradius, m_apexradius, m_height, m_slope;
    int m_resolution;
};

};

#endif /* __CVK_CONE_H */

