#include "CPU_BVH.h"
#include "KIRK/Utils/Clock.h"

namespace KIRK {

CPU::BVH::BVH(const int leaf_threshold, const int max_depth) :
        TreeAccelBase(leaf_threshold, max_depth)
{
}

const char *CPU::BVH::printName() const
{
    return "BVH";
}

    void CPU::BVH::addBaseDataStructure(KIRK::CPU::Scene *scene)
{
	Clock<> clock;
    const std::vector<KIRK::Object *> objects = scene->getSceneObjects();

    // Instead of copying object pointers, refer to objects with their index in the original pointer array.
    // The ID of an object is its index in the object array
    // When we sort IDs, we can simply refer to objects contained in one node through begin/end indices
    std::vector<unsigned int> object_ids;
    object_ids.reserve(objects.size());

    // We are going to test centroids against planes instead of whole triangles.
    BoundingBox scene_centbox; // bounding box of all objects' centroids

    for(unsigned int id = 0; id < objects.size(); ++id)
    {
        scene_centbox.expandToInclude(objects[id]->getCentroid());
        object_ids.push_back(id); // Naturally IDs will be in ascending order now
    }

    m_root.split(IdRange(0, object_ids.size() - 1), objects, object_ids, scene_centbox, m_prop, 1);
	LOG_INFO("CPU BVH Construction took % milliseconds", clock.getElapsedTime());

    // This is needed for the visualizer; they would be set to 0 anyway, so we might as well use them
    m_minBound = m_root.m_bvol.min();
    m_maxBound = m_root.m_bvol.max();

	printDebugInfo();
}

/*
 * @brief Find an intersection between a ray and this BVH (node)
 * @param hit an object that includes a ray pointer, but will also store the intersection point if there is an intersection
 * @return true if an intersection can be found in this node or anywhere further down the tree
 */
bool CPU::BVH::closestIntersection(KIRK::Intersection *hit)
{
    const glm::vec3 &dir = (hit->m_ray.m_direction);
    glm::vec3 inv_dir = 1.f / dir; // precompute this as it's going to be used a lot
    int dir_sign[3] = {dir[0] < 0.0f, dir[1] < 0.0f,
                       dir[2] < 0.0f}; // Precomputing ray direction signs also helps speeding up intersection tests

    float tmin, tmax;

    // Check global bounding box first to quickly handle rays not colliding with the scene at all
    if(!m_root.m_bvol.intersects(hit->m_ray, inv_dir, dir_sign, tmin, tmax))
    {
        return false;
    }

    m_root.traverse(hit, inv_dir, dir_sign, tmin, tmax);

    return (hit->m_object != 0);
}

/*
 * @brief Only check whether there is an intersection between a ray and this BVH ( node)
 * @param ray a given ray to test for intersection with the BVH and ( possibly) contained primitives
 * @param tMax maximum allowed ray parameter
 * @return true if an intersection can be found in this node or anywhere further down the tree
 */
bool CPU::BVH::isIntersection(KIRK::Ray *ray, float tMax)
{
    const glm::vec3 &dir = ray->m_direction;
    glm::vec3 inv_dir = 1.f / dir; // precompute this as it's going to be used a lot
    int dir_sign[3] = {dir[0] < 0.0f, dir[1] < 0.0f,
                       dir[2] < 0.0f}; // Precomputing ray direction signs also helps speeding up intersection tests

    float tmin, tmax;

    // Check global bounding box first to quickly handle rays not colliding with the scene at all
    if(!m_root.m_bvol.intersects(*ray, inv_dir, dir_sign, tmin, tmax))
    {
        return false;
    }

    return m_root.traverse(ray, inv_dir, dir_sign, tmin, tmax, tMax);
}

void CPU::BVHNode::split(const CPU::IdRange &range, const std::vector<KIRK::Object *> &objects,
                    std::vector<unsigned int> &object_ids, const CPU::BoundingBox &centbox, TreeAccelProperties &p,
                    int depth)
{
    m_bvol.computeBoundaries(objects, object_ids,
                             range); // This will actually contain all bounding boxes for objects, not just their centroids
    ++(p.amnt_nodes);

    // Keep splitting nodes until minimum amount of objects per leaf is reached OR centroid bounding box is too small to be split
    if(depth < p.max_depth && range.second - range.first > p.leaf_threshold && centbox.worthSplitting())
    {
        IdRange left_range, right_range; // to be filled
        BoundingBox left_centbox, right_centbox; // to be filled

        // Compute splitting plane and decide which objects should be on either side of it
        partition(range, left_range, right_range, object_ids, objects, centbox, left_centbox, right_centbox);

        // Our BVH is a binary tree - each node has two children which might need to split recursively.
        // Upon splitting, each node will compute its bounding volume regardless of whether it's a leaf node or not.
        // Keep in mind the node's bounding box is NOT the same as the centroid box, as the centroid box takes into account just the objects' bounding box centroids instead of all vertices
        LEFT_CHILD = new BVHNode;
        LEFT_CHILD->split(left_range, objects, object_ids, left_centbox, p, depth + 1);
        RIGHT_CHILD = new BVHNode;
        RIGHT_CHILD->split(right_range, objects, object_ids, right_centbox, p, depth + 1);
    }
        // Leaf node; add objects here, no further splitting
        // A little bit of explanation here because the hierarchy is confusing: Our BVH is derived from Voxel, which in turn is derived from Container. This Container class implements the method "addCandidate" we are using here. Our tree itself is not concerned with storing and testing ray intersection with primitives, we delegate this responsibility to the underlying Container class
    else
    {
        ++(p.amnt_leaves);

        // Remember highest depth for debugging. Simply check whether this leaf node has the current highest depth
        if(depth > p.depth)
        {
            p.depth = depth;
        }

        // Get the actual pointer corresponding to each ID and store it.
        for(unsigned int id = range.first; id <= range.second; ++id)
        {
            addCandidate(objects[object_ids[id]]);
        }
    }
}

/*
 * @brief Traverse the BVH in order to find intersections with primitives stored in leaf nodes
 * @param hit an object that includes a ray pointer, but will also store the intersection point if there is an intersection
 * @param inv_dir the inverted direction of the input ray
 * @param dir_sign ray direction signs on each axis
 * @param tmin ray entry point in the current node's bounding box
 * @param tmax ray exit point in the current node's bounding box
 */
void CPU::BVHNode::traverse(KIRK::Intersection *hit, const glm::vec3 &inv_dir, int dir_sign[3], const float tmin, const float tmax)
{
    // If the entry point is already further away from the camera than the current best intersection, this can be skipped
    if(tmax < 0.f || tmin > hit->m_lambda)
    { // ... or if exit point is negative
        return;
    }

    // If this is a leaf node, check intersection with contained primitives
    if(hasCandidates())
    {
		KIRK::Intersection tmp_hit(
                hit->m_ray); // using temporary hit object so that we don't accidentally overwrite the current closest intersection if this is further away

        if(closestIntersectionWithCandidates(&tmp_hit, 0, tmax))
        {
            if(tmp_hit.m_lambda < hit->m_lambda)
            {
                *hit = tmp_hit; // Overwrite current closest intersection only if this is even closer
            }
        }
    } else
    {
        // Check for intersection with both children; if any child's bounding volume is not intersected, we can skip its entire subtree
        float left_tmin, left_tmax, right_tmin, right_tmax;
        bool left_hit = LEFT_CHILD->m_bvol.intersects(hit->m_ray, inv_dir, dir_sign, left_tmin, left_tmax);
        bool right_hit = RIGHT_CHILD->m_bvol.intersects(hit->m_ray, inv_dir, dir_sign, right_tmin, right_tmax);

        if(left_hit && right_hit)
        { // traverse both children
            if(left_tmin < right_tmin)
            { // left child is closer
                LEFT_CHILD->traverse(hit, inv_dir, dir_sign, left_tmin, left_tmax);
                RIGHT_CHILD->traverse(hit, inv_dir, dir_sign, right_tmin, right_tmax);
            } else
            { // right child is closer
                RIGHT_CHILD->traverse(hit, inv_dir, dir_sign, right_tmin, right_tmax);
                LEFT_CHILD->traverse(hit, inv_dir, dir_sign, left_tmin, left_tmax);
            }
        } else
        {
            if(left_hit)
            { // traverse left child only
                LEFT_CHILD->traverse(hit, inv_dir, dir_sign, left_tmin, left_tmax);
            } else if(right_hit)
            { // traverse right child only
                RIGHT_CHILD->traverse(hit, inv_dir, dir_sign, right_tmin, right_tmax);
            }
        }
    }
    // else, no intersection with this node's bounding volume, hence no intersection with any primitive contained in either this node or its children. Continue in another branch.
}

/*
 * @brief Traverse the BVH in order to find any intersection with primitives stored in leaf nodes
 * @param ray a given ray to test for intersection with the BVH and ( possibly) contained primitives
 * @param inv_dir the inverted direction of the input ray
 * @param dir_sign ray direction signs on each axis
 * @param tmin ray entry point in the current node's bounding box
 * @param tmax ray exit point in the current node's bounding box
 * @param tMax maximum allowed ray parameter
 * @return true if an intersection can be found in this node or anywhere further down the tree
 */
bool CPU::BVHNode::traverse(KIRK::Ray *ray, const glm::vec3 &inv_dir, int dir_sign[3], const float tmin, const float tmax,
                       const float tMax)
{
    // If the entry point is further away than the maximum allowed entry point, this node is not intersected.
    if(tmax < 0.f || tmin > tMax)
    { // ... or if exit point is negative
        return false;
    }

    // If this is a leaf node, check intersection with contained primitives
    if(hasCandidates())
    {
        if(isIntersectionWithCandidates(ray, tMax))
        {
            return true; // Any intersection is enough
        }
    } else
    {
        // Check for intersection with both children; if any child's bounding volume is not intersected, we can skip its entire subtree
        float left_tmin, left_tmax, right_tmin, right_tmax;
        bool left_hit = LEFT_CHILD->m_bvol.intersects(*ray, inv_dir, dir_sign, left_tmin, left_tmax);
        bool right_hit = RIGHT_CHILD->m_bvol.intersects(*ray, inv_dir, dir_sign, right_tmin, right_tmax);

        if(left_hit && right_hit)
        { // traverse both children
            if(left_tmin < right_tmin)
            { // left child is closer
                if(LEFT_CHILD->traverse(ray, inv_dir, dir_sign, left_tmin, left_tmax, tMax))
                    return true;
                if(RIGHT_CHILD->traverse(ray, inv_dir, dir_sign, right_tmin, right_tmax, tMax))
                    return true;
            } else
            { // right child is closer
                if(RIGHT_CHILD->traverse(ray, inv_dir, dir_sign, right_tmin, right_tmax, tMax))
                    return true;
                if(LEFT_CHILD->traverse(ray, inv_dir, dir_sign, left_tmin, left_tmax, tMax))
                    return true;
            }
        } else
        {
            if(left_hit)
            { // traverse left child only
                if(LEFT_CHILD->traverse(ray, inv_dir, dir_sign, left_tmin, left_tmax, tMax))
                    return true;
            } else if(right_hit)
            { // traverse right child only
                if(RIGHT_CHILD->traverse(ray, inv_dir, dir_sign, right_tmin, right_tmax, tMax))
                    return true;
            }
        }
    }
    // else, no intersection with this node's bounding volume, hence no intersection with any primitive contained in either this node or its children. Continue in another branch.

    return false;
}

void CPU::BVHNode::partition(const IdRange &id_range, IdRange &left_range, IdRange &right_range,
                        std::vector<unsigned int> &object_ids, const std::vector<KIRK::Object *> &objects,
                        const BoundingBox &centbox, BoundingBox &left_centbox, BoundingBox &right_centbox)
{
/*
    ////////////////////////////////////////////
    //   M  E  D  I  A  N     S  P  L  I  T   //
    ////////////////////////////////////////////
    const glm::vec3& min = centbox.min ();
    const glm::vec3& max = centbox.max ();
    glm::vec3 extent = max - min;

    int split_dim = 0;

    if ( extent.y > extent.x ) {
        split_dim = 1;
    }

    if ( extent.z > extent.y ) {
        split_dim = 2;
    }

    auto begin = std::next ( object_ids.begin (), id_range.first );
    auto end = std::next ( object_ids.begin (), id_range.second );

    std::sort ( begin, end, [ &objects, split_dim ] ( unsigned int a, unsigned int b ) { return objects[ a ]->getCentroid ()[ split_dim] < objects[ b ]->getCentroid ()[ split_dim ]; } );

    unsigned int median = ( id_range.second - id_range.first + 1 ) / 2;

    left_range.first = id_range.first;
    left_range.second = id_range.first + median;
    right_range.first = id_range.first + median + 1;
    right_range.second = id_range.second;

    for ( unsigned int id = left_range.first; id <= left_range.second; ++id ) {
        left_centbox.expandToInclude ( objects[ object_ids[ id ] ]->getCentroid () );
    }

    for ( unsigned int id = right_range.first; id <= right_range.second; ++id ) {
        right_centbox.expandToInclude ( objects[ object_ids[ id ] ]->getCentroid () );
    }

*/
/*
    //////////////////////////////////////////////////
    //   M  I  D  C  O  O  R  D     S  P  L  I  T   //
    //////////////////////////////////////////////////
    const glm::vec3& min = centbox.min ();
    const glm::vec3& max = centbox.max ();
    glm::vec3 extent = max - min;

    int split_dim = 0;

    if ( extent.y > extent.x ) {
        split_dim = 1;
    }

    if ( extent.z > extent.y ) {
        split_dim = 2;
    }

    float split_coord = .5f * ( min[ split_coord ] + max[ split_coord ] );

    unsigned int mid = id_range.first;

    for ( unsigned int id = id_range.first; id <= id_range.second; ++id ) {
        if ( centroids[ object_ids[ id ] ][ split_dim ] < split_coord ) {
            std::swap ( object_ids[ id ], object_ids[ mid ] );
            ++mid;
        }
    }

    if ( mid == id_range.first || mid == id_range.second ) {
        mid = id_range.first + ( ( id_range.second - id_range.first ) / 2 );
    }

    left_range.first = id_range.first;
    left_range.second = mid - 1;
    right_range.first = mid;
    right_range.second = id_range.second;

    for ( unsigned int id = left_range.first; id <= left_range.second; ++id ) {
        left_centbox.expandToInclude ( centroids[ object_ids[ id ] ] );
    }

    for ( unsigned int id = right_range.first; id <= right_range.second; ++id ) {
        right_centbox.expandToInclude ( centroids[ object_ids[ id ] ] );
    }

*/
    /////////////////////////////////////////
    //   B I N N E D   S A H   S P L I T   //
    /////////////////////////////////////////
    float best = std::numeric_limits<float>::max();
    int best_axis = 0;
    int best_plane = 0;
    const int amnt_bins = 16;
    const int amnt_planes = amnt_bins - 1;

    AxisConstants axis_constants[3]; // for each axis, there's precomputations that we can make to make things faster. We still need them later once we've found the best axis, hence we declare this before iterating the axis.

    // The algorithm used has a couple of "unusual properties" for maximum efficiency:
    // 1. Ignore perfect splits - only consider object AABB's centroids to determine which side an object is on
    // 2. Do not test all potential split planes, but K equidistantly placed ones (binned sah)
    // 3. The SAH cost function can be simplified, as some terms are common for all planes.
    // General idea: For each axis, project objects into K+1 bins formed by the K split planes. Each bin counts the number of objects that overlap it. After computing the SAH cost for each of these planes for each axis, the best one is selected, and a second linear pass over all objects generates the actual list of objects IDs for the left and right subtrees.
    // Conclusion: What we need is the amount of triangles to the right and to the left of each plane, as well as the surface area of the accumulated bounding box on each side (everything else can be simplified and/or precomputed)
    for(int axis = 0; axis < 3; ++axis)
    {
        // Precompute some values that are constant for each axis
        // The "epsilon" ensures that centroids on exactly the right bounds of "centbox" still project to bin K-1
        const float cbmin = centbox.min()[axis];
        const float cbmax = centbox.max()[axis];
        const float cbdiff = cbmax - cbmin;
        const float epsilon = 0.1;
        const float k = (amnt_bins * (1 - epsilon)) / cbdiff;
        axis_constants[axis] = AxisConstants(cbmin, cbmax, cbdiff,
                                             k); // stores these for later on (as these local values will run out of scope)

        BoundingBox bin_bounds[amnt_bins]; // bounding box of each bin
        unsigned int amnt_per_bin[amnt_bins]; // amount of objects in each bin

        // Initialize with 0 since we haven't started adding objects to bins
        for(unsigned int i = 0; i < amnt_bins; ++i)
        {
            amnt_per_bin[i] = 0;
        }

        // Project objects into their respective bins
        // Iterate over the range of objects contained in the current node that needs to be split
        for(unsigned int id = id_range.first; id <= id_range.second; ++id)
        {
            unsigned int object_id = object_ids[id];
            int bin_id = k * ((objects[object_id]->getCentroid())[axis] - cbmin); // truncate result
            bin_bounds[bin_id].expandToInclude(
                    objects[object_id]->getCentroid()); // calculate bounds of each bin containing centroids of objects belonging in this bin
            ++(amnt_per_bin[bin_id]);
        }

        // First pass: Incrementally accumulate the bounds and number of triangles on the left half
        unsigned int left_amnt_objects[amnt_planes];
        BoundingBox left_bounds[amnt_planes];

        // Check this out: the bin to the left of each plane has the same index as the plane itself:
        // 0 | 1 | 2   <-- bin indices
        //   0   1    <-- plane indices
        // Naturally, since there's one more bin than there's planes, the last plane will not consider the last bin (as it lies to the right)
        // Growing the left bounding box incrementally means we need to take into account the bin's bounding box as well as the bounding box of the previous plane
        // Same goes for the object count
        left_bounds[0].expandToInclude(bin_bounds[0]);
        left_amnt_objects[0] = amnt_per_bin[0];

        for(int plane = 1; plane < amnt_planes; ++plane)
        {
            left_bounds[plane].expandToInclude(left_bounds[plane - 1]);
            left_bounds[plane].expandToInclude(bin_bounds[plane]);
            left_amnt_objects[plane] = left_amnt_objects[plane - 1] + amnt_per_bin[plane];
        }

        // Second pass: Incrementally accumulate the bounds and number of triangles on the right half, starting on the right side.
        // Now we have information on both sides, so we can compute the SAH costs in the third pass as well.
        // As for incremental computation, refer to the image above once more.
        unsigned int right_amnt_objects[amnt_planes];
        BoundingBox right_bounds[amnt_planes];

        for(int plane = amnt_planes - 1; plane >= 0; --plane)
        {
            // plane+1 is actually a valid index if plane=amnt_planes-1 as there's one more bin than there are planes
            right_bounds[plane].expandToInclude(bin_bounds[plane + 1]);
            right_amnt_objects[plane] = amnt_per_bin[plane + 1];

            // There is no plane to the right of the rightmost plane...
            if(plane != amnt_planes - 1)
            {
                right_bounds[plane].expandToInclude(right_bounds[plane + 1]);
                right_amnt_objects[plane] += right_amnt_objects[plane + 1];
            }

            // At this point, we have computed the total amount of objects on either side of this plane (bounds, too)
            // We can now compute the sah cost.
            float surface_area_left = left_bounds[plane].surfaceArea();
            float surface_area_right = right_bounds[plane].surfaceArea();
            float cost = surface_area_left * left_amnt_objects[plane] + surface_area_right * right_amnt_objects[plane];

            if(cost < best)
            {
                best = cost;
                best_axis = axis;
                best_plane = plane;
                left_centbox = left_bounds[plane];
                right_centbox = right_bounds[plane];
            }
        }

    }

    // In order to rearrange the primitive IDs to match the new node layout, we use one iterator for each side (start/end) approaching one another, swapping elements if necessary.
    const float cbmin = axis_constants[best_axis].cbmin;
    const float k = axis_constants[best_axis].k;
    int left = id_range.first;
    int right = id_range.second;
    bool left_stopped = false;
    bool right_stopped = false;

    while(left < right)
    {
        if(!left_stopped)
        {
            int bin_id_left = k * ((objects[object_ids[left]]->getCentroid())[best_axis] - cbmin);

            // If the bin is to the right of the splitting plane, we need to get this object on the other side; stop and swap once the right iterator has found an unfitting object as well.
            if(bin_id_left > best_plane)
            {
                left_stopped = true;
            } else
            {
                ++left; // This object is on the right side, continue searching
            }
        }

        if(!right_stopped)
        {
            int bin_id_right = k * ((objects[object_ids[right]]->getCentroid())[best_axis] - cbmin);

            if(bin_id_right <= best_plane)
            {
                right_stopped = true;
            } else
            {
                --right;
            }
        }

        // If both iterators are stopped, swap elements so that they are both on the appropriate side
        if(left_stopped && right_stopped)
        {
            // unsigned int tmp = object_ids[ left ];
            // object_ids[ left ] = object_ids[ right ];
            // object_ids[ right ] = tmp;

            std::swap(object_ids[left], object_ids[right]);
            left_stopped = false;
            right_stopped = false;
            ++left;
            --right;
        }
    }

    left_range.first = id_range.first;
    right_range.second = id_range.second;

    // If the iterators passed each other by each incrementing at the same time, both must've been fine with the bin id at their previous location - each is now on the begin/end index of their opponent
    if(left > right)
    {
        left_range.second = right;
        right_range.first = left;
    } else
    { // left == right
        // If left is stopped, right is not stopped. The object at "left" belongs to the right side - hence it marks the beginning of the right range.
        if(left_stopped)
        {
            left_range.second = left - 1;
            right_range.first = left;
        }
            // The other way round
        else if(right_stopped)
        {
            left_range.second = right;
            right_range.first = right + 1;
        } else
        { // !left_stopped && !right_stopped
            // No iterator is stopped; we don't know what's up with the last index. Let's find out and set bounds correctly.
            int bin_id_left = k * ((objects[object_ids[left]]->getCentroid())[best_axis] - cbmin);

            if(bin_id_left > best_plane)
            {
                left_range.second = left - 1;
                right_range.first = left;
            } else
            {
                left_range.second = left;
                right_range.first = left + 1;
            }
        }
    }
}

}
