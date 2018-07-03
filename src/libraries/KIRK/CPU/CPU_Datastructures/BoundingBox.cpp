#include "BoundingBox.h"
#include <sstream>

#define MIN m_bounds[0]
#define MAX m_bounds[1]

namespace KIRK {
CPU::BoundingBox::BoundingBox()
{
    MIN = glm::vec3(std::numeric_limits<float>::max());
    MAX = glm::vec3(std::numeric_limits<float>::lowest());
}

CPU::BoundingBox::BoundingBox(const CPU::BoundingBox &other)
{
    MIN = other.MIN;
    MAX = other.MAX;
}

CPU::BoundingBox::BoundingBox(const glm::vec3 &min, const glm::vec3 &max)
{
    MIN = min;
    MAX = max;
}

void CPU::BoundingBox::expandToInclude(const BoundingBox &rhs)
{
    // Need to go through all dimensions (axis) and see if we might have to expand in that direction
    for(int d = 0; d < 3; ++d)
    {
        MIN[d] = std::min(MIN[d], rhs.MIN[d]);
        MAX[d] = std::max(MAX[d], rhs.MAX[d]);
    }
}

void CPU::BoundingBox::expandToInclude(const glm::vec3 &point)
{
    for(int d = 0; d < 3; ++d)
    {
        MIN[d] = std::min(MIN[d], point[d]);
        MAX[d] = std::max(MAX[d], point[d]);
    }
}

void CPU::BoundingBox::expandToInclude(glm::vec3 *bbox)
{
    for(int d = 0; d < 3; ++d)
    {
        MIN[d] = std::min(MIN[d], bbox[0][d]);
        MAX[d] = std::max(MAX[d], bbox[1][d]);
    }
}

bool CPU::BoundingBox::contains(const CPU::BoundingBox &other) const
{
    for(int d = 0; d < 3; ++d)
    {
        if(other.MAX[d] < MIN[d] || other.MIN[d] > MAX[d])
        {
            return false;
        }
    }

    return true;
}

bool CPU::BoundingBox::contains(const glm::vec3 *other) const
{
    for(int d = 0; d < 3; ++d)
    {
        if(other[1][d] < MIN[d] || other[0][d] > MAX[d])
        {
            return false;
        }
    }

    return true;
}

const glm::vec3 &CPU::BoundingBox::min() const
{
    return MIN;
}

const glm::vec3 &CPU::BoundingBox::max() const
{
    return MAX;
}

glm::vec3 CPU::BoundingBox::extend() const
{
    return MAX - MIN;
}

float CPU::BoundingBox::surfaceArea() const
{
    glm::vec3 size = MAX - MIN;

    return 2.f * (size.x * size.y + size.x * size.z + size.y * size.z);
}

bool CPU::BoundingBox::worthSplitting() const
{
    glm::vec3 diff = MAX - MIN;

    return diff.x > 0.f && diff.y > 0.f && diff.z > 0.f;
}

CPU::BoundingVolume::BoundingVolume() : CPU::BoundingBox()
{
}

CPU::BoundingVolume::BoundingVolume(const glm::vec3 &min, const glm::vec3 &max) : BoundingBox(min, max)
{
}

void CPU::BoundingVolume::computeBoundaries(const std::vector<KIRK::Object *> &objects,
                                       const std::vector<unsigned int> &object_ids, const IdRange &id_range)
{
	
	if (objects.size() == 0)
		throw std::invalid_argument("Your scene is empty!");
	

    // Get all vertices from objects' bounding boxes and push them into a single vector of vertices
    for(unsigned int id = id_range.first; id <= id_range.second; ++id)
    {
        glm::vec3 *bbox = objects[object_ids[id]]->getBounds();

        // bbox[ 0 ] are the lowest coordinates of the object's bounding box
        // bbox[ 1 ] are the largest coordinates of the object's bounding box
        for(int d = 0; d < 3; ++d)
        {
            MIN[d] = std::min(MIN[d], bbox[0][d]);
            MAX[d] = std::max(MAX[d], bbox[1][d]);
        }

        // MIN, MAX are set - we're done. These are two opposite corners of the bounding box out of which you can compute all other vertices of the bounding box much the same way as we did with the objects' bounding boxes.
    }
}

bool CPU::BoundingVolume::intersects(const KIRK::Ray &ray, const glm::vec3 &inv_dir, int ray_dir_sign[3], float &tmin,
                                float &tmax) const
{
    const glm::vec3 origin = ray.m_origin;
    // Usually, this would be something like:
    // double t1 = ( MIN[ 0 ] - origin[ 0 ])/dir[ 0 ]; )
    // Spock says divisions are unwise, hence we use the precomputed inverse direction to get rid of it.
    // Also, using ray direction signs for each coordinate can help speeding things up.
    // Go over all coordinates of our minimum and maximum boundaries to choose the point for which tmin is the greatest (vice versa for tmax) since not all intersecting points lie on the bounding volume (remember we're testing against planes). Cyrus-Beck clipping does a very similar thing to ensure it's using the right clipping points.

    // tmax must be bigger than tmin for the ray to intersect. There are some nice illustrations on the internet showing why this works. This is also similar to Cyrus-Beck clipping.

    float tymin, tymax, tzmin, tzmax;
    tmin = (m_bounds[ray_dir_sign[0]].x - origin.x) * inv_dir.x;
    tmax = (m_bounds[1 - ray_dir_sign[0]].x - origin.x) * inv_dir.x;
    tymin = (m_bounds[ray_dir_sign[1]].y - origin.y) * inv_dir.y;
    tymax = (m_bounds[1 - ray_dir_sign[1]].y - origin.y) * inv_dir.y;

    if((tmin > tymax) || (tymin > tmax))
    {
        return false;
    }

    if(tymin > tmin)
    {
        tmin = tymin;
    }

    if(tymax < tmax)
    {
        tmax = tymax;
    }

    tzmin = (m_bounds[ray_dir_sign[2]].z - origin.z) * inv_dir.z;
    tzmax = (m_bounds[1 - ray_dir_sign[2]].z - origin.z) * inv_dir.z;

    if((tmin > tzmax) || (tzmin > tmax))
    {
        return false;
    }

    if(tzmin > tmin)
    {
        tmin = tzmin;
    }

    if(tzmax < tmax)
    {
        tmax = tzmax;
    }

    return true;
}
}
