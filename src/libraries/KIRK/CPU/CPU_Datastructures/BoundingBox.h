#ifndef __CVK_DS_BB_H__
#define __CVK_DS_BB_H__

#include <utility>
#include <glm/glm.hpp>
#include <vector>
#include <algorithm>
#include "KIRK/Common/Triangle.h"
#include "KIRK/Common/Ray.h"

namespace KIRK {
namespace CPU {

typedef std::pair<unsigned int, unsigned int> IdRange;

/*
 * @brief A simple bounding box, marking its boundaries with minimum and maximum coordinates (two vectors)
 */
class BoundingBox
{
public:

	/*
	* @brief Default constructor; Make sure std::min and std::max don't return gibberish upon first expansion of the bounds!
	* Can't use initializer list because Visual Studio says no (C2536).
	*/
    BoundingBox();

	/*
	* @brief Alternative constructor; construct from another bounding box
	* @param other the other bounding box to "copy"
	*/
    BoundingBox(const BoundingBox &other);

	/*
	* @brief Alternative constructor.
	* @param min initial minimum bounding coordinates
	* @param max initial maximum bounding coordinates
	*/
    BoundingBox(const glm::vec3 &min, const glm::vec3 &max);

	/*
	* @brief This can be thought of as a "grow to include"-function.
	* @param rhs bounding box used to grow this bounding box
	* @return this bounding box, but it's (probably) bigger and now contains rhs
	*/
    void expandToInclude(const BoundingBox &rhs);

	/*
	* @brief Given a point in 3D-space, grow this bounding box to fully contain the point.
	* @param point the point to be contained in this bounding box
	*/
    void expandToInclude(const glm::vec3 &point);

	/*
	* @brief Another overload for growing the box so that it includes another box
	* @param bbox a bounding box as an array of two vectors (min,max coords)
	*/
    void expandToInclude(glm::vec3 *bbox);

	/*
	* @brief Check whether a bounding box is contained within this bounding box
	* @param other the other bounding box to test
	* @returns true if the other bounding box lies fully within this bounding box
	*/
    bool contains(const BoundingBox &other) const;

	/*
	* @brief Same as above
	* @param other the other bounding box as an array of glm::vec3s
	* @return same as above
	*/
    bool contains(const glm::vec3 *other) const;

	/*
	* @return the minimum bounds of this bounding box
	*/
    const glm::vec3 &min() const;

	/*
	* @return the maximum bounds of this bounding box
	*/
    const glm::vec3 &max() const;

	/*
	* @return the extend, i.e. the size of this bounding box.
	*/
    glm::vec3 extend() const;

	/*
	* @return the surface area of this bounding box
	*
	* This is used in the computation of the SAH-cost for a node split
	*/
    float surfaceArea() const;

	/*
	* @return True if the bounding box is big enough to be worth splitting into smaller ones during tree construction
	*/
    bool worthSplitting() const;
protected:
    glm::vec3 m_bounds[2]; // m_bounds[ 0 ] = minimum, m_bounds [ 1 ] = maximum
};

/*
 * @brief A specific bounding box marking the boundaries of a node in a Bounding Volume Hierarchy
 */
class BoundingVolume : public BoundingBox
{
public:

	/*
	* @brief Default BoundingVolume class constructor
	*/
    BoundingVolume();

	/*
	* @brief BoundingVolume class constructor
	* @param min initial minimum bounding coordinates
	* @param max initial maximum bounding coordinates
	*/
    BoundingVolume(const glm::vec3 &min, const glm::vec3 &max);

	/*
	* @brief Compute the boundaries of this bounding volume
	* @param objects all objects
	* @param object_ids id_range refers to this array holding indices to objects, but not in ascending order
	* @param id_range range of objects contained within this bounding volume (indices only, no pointers)
	*
	* This bounding volume has eight vertices marking its boundaries. We only need to store two vertices (min and max ones) in order to be able to compute all the others.
	*/
    void computeBoundaries(const std::vector<KIRK::Object *> &objects, const std::vector<unsigned int> &object_ids,
                           const IdRange &id_range);

	/*
	* @return true if a given ray intersects this axis-aligned bounding volume
	* @param ray the given ray
	* @param inv_dir precomputed inverse direction of ray
	* @param ray_dir_sign ray direction signs on each axis
	* @param tmin reference; overwrite with ray parameter at entry intersection point
	* @param tmax reference; overwrite with ray parameter at exit intersection point
	*/
    bool intersects(const KIRK::Ray &ray, const glm::vec3 &inv_dir, int ray_dir_sign[3], float &tmin, float &tmax) const;
};
}}

#endif
