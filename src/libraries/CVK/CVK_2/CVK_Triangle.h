#ifndef __CVK_TRIANGLE_H
#define __CVK_TRIANGLE_H

#include "CVK_Defs.h"
#include "CVK_Geometry.h"

namespace CVK {

/**
 * The Triangle is one of the example geometry classes in the CVK. The normals of
 * the vertices don't need to be orthogonal to the plane of the Triangle.
 * @brief Triangle class as Geometry
 */
class Triangle : public CVK::Geometry
{
public:
    /**
     * Standard Constructor for Triangle
     */
    Triangle();
    /**
    * Constructor for Triangle with given parameters
    * @param a 1. vertex position
    * @param b 2. vertex position
    * @param c 3. vertex position
    */
    Triangle(glm::vec3 a, glm::vec3 b, glm::vec3 c);
    /**
    * Constructor for Triangle with given parameters
    * @param a 1. vertex position
    * @param b 2. vertex position
    * @param c 3. vertex position
    * @param na 1. vertex normal
    * @param nb 2. vertex normal
    * @param nc 3. vertex normal
    */
    Triangle(glm::vec3 a, glm::vec3 b, glm::vec3 c, glm::vec3 na, glm::vec3 nb, glm::vec3 nc);
    /**
     * Constructor for Triangle with given parameters
     * @param a 1. vertex position
     * @param b 2. vertex position
     * @param c 3. vertex position
     * @param tca 1. vertex uv coordinate
     * @param tcb 2. vertex uv coordinate
     * @param tcc 3. vertex uv coordinate
     */
    Triangle(glm::vec3 a, glm::vec3 b, glm::vec3 c, glm::vec2 tca, glm::vec2 tcb, glm::vec2 tcc);
    /**
     * Constructor for Triangle with given parameters
     * @param a 1. vertex position
     * @param b 2. vertex position
     * @param c 3. vertex position
     * @param na 1. vertex normal
     * @param nb 2. vertex normal
     * @param nc 3. vertex normal
     * @param tca 1. vertex uv coordinate
     * @param tcb 2. vertex uv coordinate
     * @param tcc 3. vertex uv coordinate
     */
    Triangle(glm::vec3 a, glm::vec3 b, glm::vec3 c, glm::vec3 na, glm::vec3 nb, glm::vec3 nc, glm::vec2 tca,
             glm::vec2 tcb, glm::vec2 tcc);

    /**
     * Standard Destructor for Triangle
     */
    ~Triangle();

    /**
     * @brief Setter for vertex positions
     * @param a 1. vertex position
     * @param b 2. vertex position
     * @param c 3. vertex position
     */
    void set_Points(glm::vec3 a, glm::vec3 b, glm::vec3 c);
    /**
    * @brief Getter for vertex positions as pointers
    * @param a 1. vertex position as pointer
    * @param b 2. vertex position as pointer
    * @param c 3. vertex position as pointer
    */
    void get_Points(glm::vec3 *a, glm::vec3 *b, glm::vec3 *c);

    /**
     * @brief Standard Setter for plane normal
     * @param n the new normal of this object
     */
    void set_Normal(glm::vec3 n);
    /**
    * @brief Standard Getter for plane normal
    * @param n the normal of this object as pointer
    */
    void get_Normal(glm::vec3 *n);

    /**
     * @brief Setter for vertex normals
     * @param na 1. vertex normal
     * @param nb 2. vertex normal
     * @param nc 3. vertex normal
     */
    void set_Normals(glm::vec3 na, glm::vec3 nb, glm::vec3 nc);
    /**
    * @brief Getter for vertex normals as pointers
     * @param na 1. vertex normal as pointer
     * @param nb 2. vertex normal as pointer
     * @param nc 3. vertex normal as pointer
    */
    void get_Normals(glm::vec3 *na, glm::vec3 *nb, glm::vec3 *nc);

    /**
    * @brief Setter for vertex uv coordinates
    * @param tca 1. vertex uv coordinate
    * @param tcb 2. vertex uv coordinate
    * @param tcc 3. vertex uv coordinate
    */
    void set_Tcoords(glm::vec2 tca, glm::vec2 tcb, glm::vec2 tcc);
    /**
    * @brief Getter for vertex uv coordinates as pointers
    * @param tca 1. vertex uv coordinate as pointer
    * @param tcb 2. vertex uv coordinate as pointer
    * @param tcc 3. vertex uv coordinate as pointer
    */
    void get_Tcoords(glm::vec2 *tca, glm::vec2 *tcb, glm::vec2 *tcc);

private:
    /**
    * Create the Triangle and the buffers with the given attributes
    * @brief Create the Triangle and the buffers
    * @param a 1. vertex position
    * @param b 2. vertex position
    * @param c 3. vertex position
    * @param tca 1. vertex uv coordinate
    * @param tcb 2. vertex uv coordinate
    * @param tcc 3. vertex uv coordinate
    */
    void create(glm::vec3 a, glm::vec3 b, glm::vec3 c, glm::vec2 tca, glm::vec2 tcb, glm::vec2 tcc);
    /**
    * Create the Triangle and the buffers with the given attributes
    * @brief Create the Triangle and the buffers
    * @param a 1. vertex position
    * @param b 2. vertex position
    * @param c 3. vertex position
    * @param na 1. vertex normal
    * @param nb 2. vertex normal
    * @param nc 3. vertex normal
    * @param tca 1. vertex uv coordinate
    * @param tcb 2. vertex uv coordinate
    * @param tcc 3. vertex uv coordinate
    */
    void create(glm::vec3 a, glm::vec3 b, glm::vec3 c, glm::vec3 na, glm::vec3 nb, glm::vec3 nc, glm::vec2 tca,
                glm::vec2 tcb, glm::vec2 tcc);

};

};

#endif /* __CVK_TRIANGLE_H */

