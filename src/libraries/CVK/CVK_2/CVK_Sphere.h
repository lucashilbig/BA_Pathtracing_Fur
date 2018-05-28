#ifndef __CVK_SPHERE_H
#define __CVK_SPHERE_H

#include "CVK_Defs.h"
#include "CVK_Geometry.h"

namespace CVK {

/**
* The sphere is one of the example geometry classes in the CVK.
* @brief Sphere class as Geometry
*/
class Sphere : public CVK::Geometry
{
public:
    /**
     * Standard Constructor for Sphere
     */
    Sphere();
    /**
     * Constructor for Sphere with given parameters around (0/0/0)
     * @param radius the radius of the sphere
     */
    Sphere(float radius);
    /**
    * Constructor for Sphere with given parameters around center
    * @param center the center of the sphere
    * @param radius the radius of the sphere
    */
    Sphere(glm::vec3 center, float radius);
    /**
    * Constructor for Sphere with given parameters around (0/0/0)
    * @param radius the radius of the sphere
    * @param resolution The resolution for triangulation
    */
    Sphere(float radius, int resolution);
    /**
    * Constructor for Sphere with given parameters around center
    * @param center the center of the sphere
    * @param radius the radius of the sphere
    * @param resolution The resolution for triangulation
    */
    Sphere(glm::vec3 center, float radius, int resolution);
    /**
     * Standard Destructor for Sphere
     */
    ~Sphere();

    /**
     * @brief Standard Getter for center position
     * @return the center position of this object as pointer
     */
    glm::vec3 *get_center();
    /**
     * @brief Standard Getter for radius
     * @return the radius of this object
     */
    float get_radius();

private:
    /**
    * Create the Sphere and the buffers with the given attributes
    * @brief Create the Sphere and the buffers
    * @param center The center point of the sphere
    * @param radius The radius of the sphere
    * @param resolution The resolution for triangulation
    */
    void create(glm::vec3 center, float radius, int resolution);

    glm::vec3 m_Center;
    float m_radius;
    int m_resolution;
};

};

#endif /* __CVK_SPHERE_H */

