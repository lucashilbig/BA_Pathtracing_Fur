#ifndef __CVK_PLANE_H
#define __CVK_PLANE_H

#include "CVK_Defs.h"
#include "CVK_Geometry.h"

namespace CVK {
/**
* The plane is one of the example geometry classes in the CVK. The normal is othogonal to the plane.
* @brief Plane class as Geometry
*/
class Plane : public CVK::Geometry
{
public:
    /**
     * Standard Constructor for Plane
     */
    Plane();
    /**
     * Constructor for Plane with given Parameters
     * @param a 1. vertex position
     * @param b 2. vertex position
     * @param c 3. vertex position
     * @param d 4. vertex position
     */
    Plane(glm::vec3 a, glm::vec3 b, glm::vec3 c, glm::vec3 d);
    /**
     * Constructor for Plane with given parameters
     * @param a 1. vertex position
     * @param b 2. vertex position
     * @param c 3. vertex position
     * @param d 4. vertex position
     * @param tca 1. vertex uv coordinate
     * @param tcb 2. vertex uv coordinate
     * @param tcc 3. vertex uv coordinate
     * @param tcd 4. vertex uv coordinate
     */
    Plane(glm::vec3 a, glm::vec3 b, glm::vec3 c, glm::vec3 d, glm::vec2 tca, glm::vec2 tcb, glm::vec2 tcc,
          glm::vec2 tcd);
    /**
     * Standard Destructor for Plane
     */
    ~Plane();

    /**
     * @brief Setter for vertex positions
     * @param a 1. vertex position
     * @param b 2. vertex position
     * @param c 3. vertex position
     * @param d 4. vertex position
     */
    void set_Points(glm::vec3 a, glm::vec3 b, glm::vec3 c, glm::vec3 d);
    /**
     * @brief Setter for vertex uv coordinates
     * @param tca 1. vertex uv coordinate
     * @param tcb 2. vertex uv coordinate
     * @param tcc 3. vertex uv coordinate
     * @param tcd 4. vertex uv coordinate
     */
    void set_Tcoords(glm::vec2 tca, glm::vec2 tcb, glm::vec2 tcc, glm::vec2 tcd);

private:
    /**
     * Create the Plane and the buffers with the given attributes. The normal is othogonal to the plane.
     * @brief Create the Plane and the buffers
     * @param a 1. vertex position
     * @param b 2. vertex position
     * @param c 3. vertex position
     * @param d 4. vertex position
     * @param tca 1. vertex uv coordinate
     * @param tcb 2. vertex uv coordinate
     * @param tcc 3. vertex uv coordinate
     * @param tcd 4. vertex uv coordinate
     */
    void create(glm::vec3 a, glm::vec3 b, glm::vec3 c, glm::vec3 d, glm::vec2 tca, glm::vec2 tcb, glm::vec2 tcc,
                glm::vec2 tcd);
};

};

#endif /* __CVK_PLANE_H */

