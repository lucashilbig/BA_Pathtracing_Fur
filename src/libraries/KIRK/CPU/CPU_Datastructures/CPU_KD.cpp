#include "KIRK/Utils/Log.h"
#include "CPU_KD.h"

namespace KIRK {

    const std::bitset<2> CPU::Plane::cLeft("01");
const std::bitset<2> CPU::Plane::cRight("10");
const float CPU::KDNode::cKt = 15.f; // cost for traversal step
const float CPU::KDNode::cKi = 20.f; // cost for triangle intersection

/*
 * @brief A KD-Tree used for raytracing acceleration
 * @param leaf_threshold amount of objects required to stop splitting nodes
 * @param max_depth maximum allowed total depth of tree
 * @param split_method how to split the KD-Tree; Binned or Exact.
 */
CPU::KDTree::KDTree(const int leaf_threshold, const int max_depth, const SplitMethod split_method) :
        TreeAccelBase(leaf_threshold, max_depth),
        m_split_method(split_method)
{
}

/*
 * @brief Simply returns "KD-Tree" to be used in TreeAccelBase::printDebugInfo ()
 * @return "KD-Tree"
 */
const char *CPU::KDTree::printName() const
{
    return "KD-Tree";
}

/*
 * @brief Print debug information specific to KD-Trees.
 */
void CPU::KDTree::printAdditionalDebugInfo() const
{
    LOG_DEBUG("Split method: %", ((m_split_method == SplitMethod::Binned) ? "Binned Surface Area Heuristic"
                                                                          : "Exact Surface Area Heuristic"));
}

/*
 * @brief Create the KD-Tree based on the objects provided by a given scene
 * @param scene pointer to the provided scene; scene objects are needed for tree construction
 */
void CPU::KDTree::addBaseDataStructure(CPU::Scene *scene)
{
    std::vector<KIRK::Triangle *> objects = scene->getSceneObjects();
    m_bvol = BoundingVolume((scene->getBounds())[0], (scene->getBounds())[1]);
    m_minBound = m_bvol.min();
    m_maxBound = m_bvol.max();

    // Exact SAH-splitting uses objects' bounds as splitting plane candidates
    if(m_split_method == SplitMethod::Exact)
    {
        std::vector<Plane> planes[3];

        // Instead of copying object pointers, refer to objects with their index in the original pointer array.
        // The ID of an object is its index in the object array
        // A "Plane" object stores this id, along with its position. For each object, we construct two potential splitting planes - its bounding box's min and max coordinates on a given axis.
        for(int axis = 0; axis < 3; ++axis)
        {
            planes[axis].reserve(2 * objects.size()); // two planes per object per axis

            for(unsigned int o = 0; o < objects.size(); ++o)
            {
                planes[axis].emplace_back((objects[o]->getBounds())[0][axis], o, Plane::Type::Min);
                planes[axis].emplace_back((objects[o]->getBounds())[1][axis], o, Plane::Type::Max);
            }

            // Sorting is necessary for the SAH computation later on so that we can compare plane positions (bounding box min, max coordinates) to splitting planes to determine which side an object overlaps
            std::sort(planes[axis].begin(), planes[axis].end());
        }

        m_root.split(objects, planes, BoundingBox(m_bvol.min(), m_bvol.max()), m_prop, 1);
    }
        // Use binned SAH.
    else
    {
        m_root.split(objects, BoundingBox(m_bvol.min(), m_bvol.max()), m_prop, 1);
    }
}

/*
 * @brief Find an intersection between a ray and this node
 * @param hit an object that includes a ray pointer, but will also store the intersection point if there is an intersection
 * @return true if an intersection can be found in this node or anywhere further down the tree
 */
bool CPU::KDTree::closestIntersection(KIRK::Intersection *hit)
{
    const glm::vec3 dir = hit->m_ray.m_direction;
    glm::vec3 inv_dir = 1.f / dir;
    int dir_sign[3] = {dir[0] < 0.0f, dir[1] < 0.0f, dir[2] < 0.0f};
    float tmin, tmax;

    // The only time we check intersection with a bounding box is for the root node; if this test fails, no scene intersections can be found
    if(!m_bvol.intersects(hit->m_ray, inv_dir, dir_sign, tmin, tmax))
    {
        return false;
    }

    return m_root.traverse(hit, inv_dir, tmin, tmax);
}

/*
 * @brief Only check whether there is an intersection between a ray and this node
 * @param ray a given ray to test for intersection with the tree and (possibly) contained primitives
 * @param tMax maximum allowed ray parameter
 * @return true if an intersection can be found in this node or anywhere further down the tree
 */
bool CPU::KDTree::isIntersection(KIRK::Ray *ray, float tMax)
{
    const glm::vec3 &dir = ray->m_direction;
    glm::vec3 inv_dir = 1.f / dir;
    int dir_sign[3] = {dir[0] < 0.0f, dir[1] < 0.0f, dir[2] < 0.0f};
    float tmin, tmax;

    // The only time we check intersection with a bounding box is for the root node; if this test fails, no scene intersections can be found

    if(!m_bvol.intersects(*ray, inv_dir, dir_sign, tmin, tmax) || tmin > tMax)
    {
        return false;
    }

    return m_root.traverse(ray, inv_dir, tmin, tmax, tMax);
}

/*
 * @brief A KD-Tree node; initial split axis is -1 as an indication that it has not been assigned yet
 */
CPU::KDNode::KDNode() :
        m_axis(-1)
{
}

/*
 * @brief Recursively split the KD-Tree into binary partitions
 * @param objects contains all objects in the entire scene
 * @param box bounding box of this node; in this implementation, it might not contain all objects' bounding boxes
 * @param p KD-Tree properties such as depth or amount of nodes
 * @param depth current recursion depth
 *
 * This uses binned splitting: We set a final amount of bins and sort objects into them, using bin bounds as splitting planes.
 */
void CPU::KDNode::split(std::vector<KIRK::Triangle *> &objects, const BoundingBox &box, TreeAccelProperties &p, int depth)
{
    ++(p.amnt_nodes);

    // Keep splitting nodes until maximum depth or minimum amount of objects per leaf is reached
    if(depth < p.max_depth && objects.size() > p.leaf_threshold)
    {
        BoundingBox left_box, right_box; // to be filled
        std::vector<KIRK::Triangle *> left_range, right_range; // to be filled

        // Compute splitting plane and decide which objects should be on either side of it
        // This will set m_axis and m_split_pos
        partition(left_range, right_range, objects, box, left_box, right_box);

        // Might "fail" if node is not worth splitting; make this a leaf, then.
        if(m_axis == -1)
        {
            createLeaf(objects, p, depth);

            return;
        }

        // Our KD-Tree is a binary tree - each node has two children which might need to split recursively.
        LEFT_CHILD = new KDNode;
        LEFT_CHILD->split(left_range, left_box, p, depth + 1);
        RIGHT_CHILD = new KDNode;
        RIGHT_CHILD->split(right_range, right_box, p, depth + 1);
    }
        // Leaf node; add objects here, no further splitting
        // A little bit of explanation here because the hierarchy is confusing: KDNode is derived from Container. This Container class implements the method "addCandidate" we are using here. Our tree itself is not concerned with storing and testing ray intersection with primitives, we delegate this responsibility to the underlying Container class
    else
    {
        createLeaf(objects, p, depth);
    }
}

/*
 * @brief Recursively split the KD-Tree into binary partitions (this is an alternative to the above depending on the split method used)
 * @param objects contains all objects in the entire scene
 * @param planes Basically just objects represented by their bounds as splitting planes
 * @param box bounding box of this node; in this implementation, it might not contain all objects' bounding boxes
 * @param p KD-Tree properties such as depth or amount of nodes
 * @param depth current recursion depth
 *
 * Exact SAH-splitting will consider all object bounds (min/max) as potential splitting planes.
 * Any documentation in the function body assumes that you are familiar with the Surface Area Heuristic and what is needed to compute the SAH cost for a given split plane candidate.
 */
void
CPU::KDNode::split(const std::vector<KIRK::Triangle *> &objects, std::vector<Plane> ( &planes )[3], const BoundingBox &box,
              TreeAccelProperties &p, int depth)
{
    ++(p.amnt_nodes);
    unsigned int num_triangles = planes[0].size() / 2; // two planes per object

    if(num_triangles < p.leaf_threshold || depth >= p.max_depth)
    {
        createLeaf(objects, planes, p, depth);

        return;
    }
    // make this a branch node unless we fail to find a suitable splitting plane

    float best_cost = num_triangles * cKi; // this is from the papers...
    std::vector<Plane>::iterator best_plane;

    // Try to find the best splitting plane across all axis.
    for(int axis = 0; axis < 3; ++axis)
    {
        unsigned int amnt_left = 0; // for the first axis, there are 0 objects to the left...
        unsigned int amnt_right = num_triangles; // ... and all objects are to the right of the leftmost splitting plane

        // Think of this not as going through a number of planes and choosing the best one, but as going through the objects in this voxel and testing their bounds as splitting planes.
        for(auto plane = planes[axis].begin(); plane != planes[axis].end(); ++plane)
        {
            // If this plane was the max bounds of an object, we can be sure this object does not lie to the right of this plane candidate
            if(plane->type == Plane::Type::Max)
            {
                --amnt_right;
            }

            // If split plane is out of bounds, do not consider it
            if(plane->coord >= box.min()[axis] && plane->coord <= box.max()[axis])
            {
                glm::vec3 extend = box.extend();
                float inverted_surface_area = 0.5f / (extend.x * extend.y + extend.x * extend.z + extend.y * extend.z);
                float other_extends[2] = {
                        extend[(axis + 1) % 3],
                        extend[(axis + 2) % 3]
                };

                float other_axis_square_extend = other_extends[0] * other_extends[1];
                float left_size_on_this_axis = plane->coord - box.min()[axis];
                float right_size_on_this_axis = box.max()[axis] - plane->coord;
                float left_bounding_box_area =
                        2 * (other_axis_square_extend + left_size_on_this_axis * (other_extends[0] + other_extends[1]));
                float right_bounding_box_area = 2 * (other_axis_square_extend +
                                                     right_size_on_this_axis * (other_extends[0] + other_extends[1]));
                float p_left = left_bounding_box_area * inverted_surface_area;
                float p_right = right_bounding_box_area * inverted_surface_area;
                float cost = cKt + cKi * (p_left * amnt_left + p_right * amnt_right);

                if(cost < best_cost)
                {
                    best_cost = cost;
                    best_plane = plane;
                    m_axis = axis;
                }
            }

            // If this plane was the min bounds of an object, we can now be sure the corresponding object lies to the left of the next plane candidate
            if(plane->type == Plane::Type::Min)
            {
                ++amnt_left;
            }
        }
    }

    // Eh, not worth splitting after all. Create a leaf!
    if(m_axis == -1)
    {
        createLeaf(objects, planes, p, depth);

        return;
    }

    m_split_pos = best_plane->coord;
    std::vector<std::bitset<2>> object_sides(objects.size());
    auto plane = planes[m_axis].begin();

    // If an object has its min bounds to the left of the splitting plane, we can add it to the left side.
    for(; plane != best_plane; ++plane)
    {
        if(plane->type == Plane::Type::Min)
        {
            object_sides[plane->obj_index] |= Plane::cLeft;
        }
    }

    // If an object has its max bounds to the right of the splitting plane, we can add it to the right side.
    // It might also have been added to the left side if it overlaps (meaning its min bounds were left)!
    for(++plane; plane != planes[m_axis].end(); ++plane)
    {
        if(plane->type == Plane::Type::Max)
        {
            object_sides[plane->obj_index] |= Plane::cRight;
        }
    }

    // Now we can fill new vectors with the "objects" (represented by their planes) on either side based on our calculated splitting plane
    std::vector<Plane> left[3];
    std::vector<Plane> right[3];

    // Do that for each axis as we don't know which axis will be next
    for(int axis = 0; axis < 3; ++axis)
    {
        for(const auto &plane : planes[axis])
        {
            if(Plane::objectOnSide(object_sides[plane.obj_index], Plane::cLeft))
            {
                left[axis].push_back(plane);
            }

            if(Plane::objectOnSide(object_sides[plane.obj_index], Plane::cRight))
            {
                right[axis].push_back(plane);
            }
        }
    }

    // Create the new bounds for both sides.
    // Unfortunately, these are simply clipped at the split coordinates and thus might not fully contain their objects' bounding boxes.
    glm::vec3 left_clamp = box.max();
    glm::vec3 right_clamp = box.min();
    left_clamp[m_axis] = best_plane->coord;
    right_clamp[m_axis] = best_plane->coord;
    BoundingBox left_box(box.min(), left_clamp);
    BoundingBox right_box(right_clamp, box.max());

    LEFT_CHILD = new KDNode;
    LEFT_CHILD->split(objects, left, left_box, p, depth + 1);
    RIGHT_CHILD = new KDNode;
    RIGHT_CHILD->split(objects, right, right_box, p, depth + 1);
}

/*
 * @brief Traverse the kd-tree in order to find intersections with primitives stored in leaf nodes
 * @param hit an objects that includes a ray pointer, but will also store the intersection point if there is an intersection
 * @param inv_dir the inverted direction of the input ray
 * @param tmin intersection entry point between the ray and this node
 * @param tmax intersection exit point between the ray and this node
 * @return true if an intersection can be found in this node or anywhere further down the tree
 */
bool CPU::KDNode::traverse(KIRK::Intersection *hit, const glm::vec3 &inv_dir, float tmin, float tmax)
{
    if(tmax < 0.f)
    {
        return false;
    }

    // Leaf node
    // Since bounding boxes are not allowed to overlap in a kd-tree, there is no need to test for intersections in any other node
    if(m_axis == -1)
    {
        if(!hasCandidates())
        {
            return false;
        }

        return closestIntersectionWithCandidates(hit, 0, tmax);
    }

    // Branch node
    const glm::vec3 r_origin = hit->m_ray.m_origin;
    KDNode *near;
    KDNode *far;

    // Check which node is closer to the ray origin
    if(m_split_pos > r_origin[m_axis])
    {
        near = LEFT_CHILD;
        far = RIGHT_CHILD;
    } else
    {
        near = RIGHT_CHILD;
        far = LEFT_CHILD;
    }

    // Ray intersection with splitting plane:
    double t = (m_split_pos - r_origin[m_axis]) * inv_dir[m_axis];

    if(t < 0 || t > tmax)
    {
        // near only
        if(near->traverse(hit, inv_dir, tmin, tmax))
        {
            return true;
        }
    } else if(t < tmin)
    { // && t >= 0, which holds true
        // far only
        if(far->traverse(hit, inv_dir, tmin, tmax))
        {
            return true;
        }
    } else
    {
        // both
        if(near->traverse(hit, inv_dir, tmin, t))
        {
            return true;
        }

        if(far->traverse(hit, inv_dir, t, tmax))
        {
            return true;
        }
    }

    return false;
}

/*
 * @brief Traverse the kd-tree in order to find intersections with primitives stored in leaf nodes
 * @param ray a given ray to test for intersection with the kd-tree and (possibly) contained primitives
 * @param inv_dir the inverted direction of the input ray
 * @param tmin intersection entry point between the ray and this node
 * @param tmax intersection exit point between the ray and this node
 * @param tMax maximum allowed ray parameter
 * @return true if an intersection can be found in this node or anywhere further down the tree
 */
bool CPU::KDNode::traverse(KIRK::Ray *ray, const glm::vec3 &inv_dir, float tmin, float tmax, float tMax)
{
    if(tmin > tMax || tmax < 0.f)
    {
        return false;
    }

    // Leaf node
    // Since bounding boxes are not allowed to overlap in a kd-tree, there is no need to test for intersections in any other node
    if(m_axis == -1)
    {
        if(!hasCandidates())
        {
            return false;
        }
        return isIntersectionWithCandidates(ray, tMax);
    }

    // Branch node
    const glm::vec3 &r_origin = ray->m_origin;
    KDNode *near;
    KDNode *far;

    // Check which node is closer to the ray origin
    if(m_split_pos > r_origin[m_axis])
    {
        near = LEFT_CHILD;
        far = RIGHT_CHILD;
    } else
    {
        near = RIGHT_CHILD;
        far = LEFT_CHILD;
    }

    // Ray intersection with splitting plane:
    double t = (m_split_pos - r_origin[m_axis]) * inv_dir[m_axis];

    if(t < 0 || t > tmax)
    {
        // near only
        if(near->traverse(ray, inv_dir, tmin, tmax, tMax))
        {
            return true;
        }
    } else if(t < tmin)
    { // && t >= 0, which holds true
        // far only
        if(far->traverse(ray, inv_dir, tmin, tmax, tMax))
        {
            return true;
        }
    } else
    {
        // both
        if(near->traverse(ray, inv_dir, tmin, t, tMax))
        {
            return true;
        }

        if(far->traverse(ray, inv_dir, t, tmax, tMax))
        {
            return true;
        }
    }

    return false;
}

/*
 * @brief Create a leaf node, adding objects to its underlying container
 * @param objects contains all objects in this node
 * @param p KD-Tree properties such as depth or amount of nodes
 * @param depth current recursion depth
 */
void CPU::KDNode::createLeaf(std::vector<KIRK::Triangle *> &objects, TreeAccelProperties &p, int depth)
{
    ++(p.amnt_leaves);

    // Remember highest depth for debugging.
    if(depth > p.depth)
    {
        p.depth = depth;
    }

    for(auto &object : objects)
    {
        addCandidate(object);
    }
}

/*
 * @brief Create a leaf node, adding objects to its underlying container
 * @param objects contains all objects in the entire scene
 * @param planes Basically just objects represented by their bounds as splitting planes
 * @param p KD-Tree properties such as depth or amount of nodes
 * @param depth current recursion depth
 *
 * Alternative used for exact SAH-splitting.
 */
void CPU::KDNode::createLeaf(const std::vector<KIRK::Triangle *> &objects, const std::vector<Plane> ( &planes )[3],
                        TreeAccelProperties &p, int depth)
{
    ++(p.amnt_leaves);

    // Remember highest depth for debugging.
    if(depth > p.depth)
    {
        p.depth = depth;
    }

    // Get the actual pointer corresponding to each ID and store the object in the underlying container.
    for(const auto &plane : planes[0])
    {
        if(plane.type == Plane::Type::Min)
        {
            addCandidate(objects[plane.obj_index]);
        }
    }
}

/*
 * @brief Given a set of objects, create a partition splitting these objects in two optimal sets for the KD-Tree (using an binned version of the Surface Area Heuristic)
 * @param left_range after creating the partition, should contain objects to the left of the splitting plane
 * @param right_range after creating the partition, should contain objects to the right of the splitting plane
 * @param objects contains all objects contained within the current node
 * @param box bounding box of the current node; does NOT necessarily contain all objects!
 * @param left_box results from splitting the node's bounding box in two halves
                         (can then be used as parameter "box" for the recursive partitioning of the left part of this)
 * @param right_box results from splitting the node's bounding box in two halves
                          (can then be used as parameter "box" for the recursive partitioning of the right part of this)
 *
 * Any documentation in the function body assumes that you are familiar with the Surface Area Heuristic and what is needed to compute the SAH cost for a given split plane candidate.
 */
void CPU::KDNode::partition(std::vector<KIRK::Triangle *> &left_range, std::vector<KIRK::Triangle *> &right_range,
                       std::vector<KIRK::Triangle *> &objects, const BoundingBox &box, BoundingBox &left_box,
                       BoundingBox &right_box)
{
    /////////////////////////////////////////
    //   B I N N E D   S A H   S P L I T   //
    /////////////////////////////////////////
    float best = objects.size() * cKi;
    const int amnt_bins = 64;
    const int amnt_planes = amnt_bins - 1;
    unsigned int final_amnt_left = 0;
    unsigned int final_amnt_right = 0;

    for(int axis = 0; axis < 3; ++axis)
    {
        // Precompute some values that are constant for each axis
        // These are used in the formula used to assign bins to objects
        const float bmin = box.min()[axis];
        const float bmax = box.max()[axis];
        const float bdiff = bmax - bmin;
        const float k = (amnt_bins) / bdiff;
        std::vector<unsigned int> l(amnt_bins, 0); // Collects lower bounds per bin
        std::vector<unsigned int> h(amnt_bins, 0); // Collects higher bounds per bin

        // For each plane, one pass over all objects in the current node.
        // Project object's lower and higher bounds into bins on this axis
        for(auto &object : objects)
        {
            glm::vec3 *bb = object->getBounds();
            int bin_min = k * (bb[0][axis] - bmin);
            int bin_max = k * (bb[1][axis] - bmin);

            // Since our bounding boxes are not accurate, spanning only the distance between split positions,
            // the formula will yield results that are out of range for objects with bounding boxes lying partially or completely outside of the node bounds.
            // We can simply add these to the leftmost/rightmost bins.
            if(bin_min < 0)
            {
                ++(l[0]);
            } else
            {
                ++(l[bin_min]);
            }

            if(bin_max >= amnt_bins)
            {
                ++(h[amnt_bins - 1]);
            } else
            {
                ++(h[bin_max]);
            }
        }

        // Suffix sum
        for(int i = 0; i < amnt_bins; ++i)
        {
            for(int j = i + 1; j < amnt_bins; ++j)
            {
                l[i] += l[j];
            }
        }

        // Prefix sum
        for(int i = amnt_bins - 1; i >= 0; --i)
        {
            for(int j = i - 1; j >= 0; --j)
            {
                h[i] += h[j];
            }
        }

        // Iterate over all planes once.
        // Using our l and h lists, we can determine the amount of objects on either side of each plane without further object iteration.
        for(int plane = 0; plane < amnt_planes; ++plane)
        {
            float ll = (static_cast < float > ( plane + 1 )) / amnt_bins;
            float split_pos = bmin + ll * bdiff;

            // A rare case that we must avoid. Otherwise, no triangles would be shown if the camera had a position of 0 on the splitting axis.
            if(split_pos == 0.f)
            {
                split_pos += std::numeric_limits<float>::epsilon(); // Shift by a tiny tiny amount
            }

            unsigned int amnt_left =
                    objects.size() - l[plane + 1]; // l[ plane - 1 ] is the amount of objects ONLY on the right side
            unsigned int amnt_right =
                    objects.size() - h[plane]; // h[ plane ] is the amount of objects ONLY on the left side
            glm::vec3 extend = box.extend();
            float inverted_surface_area = 0.5f / (extend.x * extend.y + extend.x * extend.z + extend.y * extend.z);
            float other_extends[2] = {
                    extend[(axis + 1) % 3],
                    extend[(axis + 2) % 3]
            };

            float other_axis_square_extend = other_extends[0] * other_extends[1];
            float left_size_on_this_axis = split_pos - bmin; // not accurate, doesn't contain all objects
            float right_size_on_this_axis = bmax - split_pos; // "
            float left_bounding_box_area =
                    2 * (other_axis_square_extend + left_size_on_this_axis * (other_extends[0] + other_extends[1]));
            float right_bounding_box_area =
                    2 * (other_axis_square_extend + right_size_on_this_axis * (other_extends[0] + other_extends[1]));
            float p_left = left_bounding_box_area * inverted_surface_area;
            float p_right = right_bounding_box_area * inverted_surface_area;
            float cost = cKt + cKi * (p_left * amnt_left + p_right * amnt_right);

            if(cost < best)
            {
                best = cost;
                m_axis = axis;
                final_amnt_left = amnt_left;
                final_amnt_right = amnt_right;
                m_split_pos = split_pos;
            }
        }
    }

    if(m_axis != -1)
    {
        left_range.reserve(final_amnt_left);
        right_range.reserve(final_amnt_right);

        // So far, we know only the amount of objects on either side, but not the actual objects.
        // Create the ranges now; doing this once afterwards is faster than doing it for each plane and copying the best ranges in the end.
        for(auto &object : objects)
        {
            glm::vec3 *bb = object->getBounds();

            if(bb[0][m_axis] < m_split_pos)
            {
                left_range.push_back(object);
            }

            if(bb[1][m_axis] > m_split_pos)
            {
                right_range.push_back(object);
            }
        }

        // The resulting voxels are NOT accurate, i.e. not the entire object bounds are included.
        // Objects that may be partially on either side will have bounding boxes extending into both sides; this bounding box will not include them.
        glm::vec3 left_clamp = box.max();
        glm::vec3 right_clamp = box.min();
        left_clamp[m_axis] = m_split_pos;
        right_clamp[m_axis] = m_split_pos;
        left_box = BoundingBox(box.min(), left_clamp);
        right_box = BoundingBox(right_clamp, box.max());
    }
    // if no plane can be found that has a cost which is lower than objects.size () * cKi, we will make this a leaf node.
}

}
