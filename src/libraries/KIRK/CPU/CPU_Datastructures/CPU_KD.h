#ifndef __CVK_DS_KD_H__
#define __CVK_DS_KD_H__

#include <bitset>
#include "TreeAccel.h"
#include "KIRK/Utils/DS_Visualizer.h"

namespace KIRK {
namespace CPU {

/*
 * @brief A splitting plane in a KD-Tree, which is either a minimum or a maximum of an object's bounding box (for exact SAH-splitting)
 *
 * This really is irrelevant if SplitType::Binned is used, and usually SplitType::Binned is what you want.
 */
struct Plane
{
    enum class Type
    {
        Min,
        Max
    };

    float coord;
    unsigned int obj_index; // remember which object's bounds are used
    Type type; // remember bound coordinate used
    static const std::bitset<2> cLeft;
    static const std::bitset<2> cRight;

    /*
     * @brief Plane constructor
     * @param pos the coordinate to split voxels with when testing how good this plane is for splitting a node
     * @param index the index of the object providing the splitting coordinate in the scene object array
     * @param type whether the coordinate used comes from minimum or maximum bounds of the object
     */
    Plane(const float pos, const unsigned int index, const Type type) : coord(pos), obj_index(index), type(type) {}

    /*
     * @return true if an object is on a certain side of a splitting plane
     * @param object_mask is either 1 (if the object is on the left side) or 2 (if the object is on the right side) or 3 (if the object is on both sides)
     * @param side either 1 (left side) or 2 (right side)
     */
    static bool objectOnSide(const std::bitset<2> object_mask, const std::bitset<2> side)
    {
        return (object_mask & side) == side;
    }

    /*
     * @brief Comparison operator overload for sorting splitting planes
     * @param rhs another splitting plane to compare to
     * @return true if this splitting plane is to the left of the other splitting plane
     */
    bool operator<(const Plane &rhs) const
    {
        return (coord == rhs.coord) ? static_cast < int > ( type ) < static_cast < int > ( rhs.type ) : coord <
                                                                                                        rhs.coord;
    }
};

/*
 * @brief A KD-Tree node
 */
class KDNode : public TreeAccelNode<KDNode, 2>
{
    template<typename T> friend
    class TreeVisualizer;

public:
    KDNode();
    void split(std::vector<KIRK::Triangle *> &objects, const BoundingBox &box, TreeAccelProperties &p, int depth);
    void split(const std::vector<KIRK::Triangle *> &objects, std::vector<Plane> ( &planes )[3], const BoundingBox &box,
               TreeAccelProperties &p, int depth);
    bool traverse(KIRK::Intersection *hit, const glm::vec3 &inv_dir, float tmin, float tmax);
    bool traverse(KIRK::Ray *ray, const glm::vec3 &inv_dir, float tmin, float tmax, float tMax);

    float getSplitPos() { return m_split_pos; }

    int getAxis() { return m_axis; }
    int m_axis;
    float m_split_pos;
    static const float cKt; // cost for traversal step
    static const float cKi; // cost for triangle intersection

private:
    void createLeaf(std::vector<KIRK::Triangle *> &objects, TreeAccelProperties &p, int depth);
    void createLeaf(const std::vector<KIRK::Triangle *> &objects, const std::vector<Plane> ( &planes )[3],
                    TreeAccelProperties &p, int depth);
    void partition(std::vector<KIRK::Triangle *> &left_range, std::vector<KIRK::Triangle *> &right_range,
                   std::vector<KIRK::Triangle *> &objects, const BoundingBox &box, BoundingBox &left_box,
                   BoundingBox &right_box);

};

/*
 * @brief A KD-Tree for raytracing acceleration
 */
class KDTree : public TreeAccelBase<KDNode>
{
    template<typename T> friend
    class TreeVisualizer;

public:
    // Algorithm used to split tree in child nodes
    enum class SplitMethod
    {
        Exact, // Fast traversal, slow build
        Binned, // Traversal usually as fast as exact, very fast build
    };

    KDTree(const int leaf_threshold = 1, const int max_depth = std::numeric_limits<int>::max(),
           const SplitMethod split_method = SplitMethod::Binned);
    virtual const char *printName() const override final;
    virtual void printAdditionalDebugInfo() const override final;
    virtual void addBaseDataStructure(KIRK::CPU::Scene *scene) override final;
    virtual bool closestIntersection(KIRK::Intersection *hit) override final;
    virtual bool isIntersection(KIRK::Ray *ray, float tMax) override final;
private:
    SplitMethod m_split_method;
    BoundingVolume m_bvol;
};
}}


#endif
