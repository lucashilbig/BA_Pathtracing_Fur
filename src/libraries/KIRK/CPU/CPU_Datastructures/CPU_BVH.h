#ifndef __CVK_DS_BVH_H__
#define __CVK_DS_BVH_H__

#include "TreeAccel.h"
#include "KIRK/Utils/DS_Visualizer.h"

namespace KIRK {
namespace CPU {
/*
 * @brief A Bounding Volume Hierarchy node
 */
class BVHNode : public TreeAccelNode<BVHNode, 2>
{
    friend class BVH;

    template<typename T> friend
    class TreeVisualizer;

public:
    void
    split(const IdRange &range, const std::vector<KIRK::Triangle *> &objects, std::vector<unsigned int> &object_ids,
          const BoundingBox &centbox, TreeAccelProperties &p, int depth);
    void traverse(KIRK::Intersection *hit, const glm::vec3 &inv_dir, int dir_sign[3], const float tmin, const float tmax);
    bool
    traverse(KIRK::Ray *ray, const glm::vec3 &inv_dir, int dir_sign[3], const float tmin, const float tmax, const float tMax);

    BoundingVolume &getBvol() { return m_bvol; }
    
    BoundingVolume m_bvol;
private:
    void
    partition(const IdRange &id_range, IdRange &left_range, IdRange &right_range, std::vector<unsigned int> &object_ids,
              const std::vector<KIRK::Triangle *> &objects, const BoundingBox &centbox, BoundingBox &left_centbox,
              BoundingBox &right_centbox);


    /*
     * @brief For keeping in mind precomputed values still needed when they go out of scope.
     * @see BVHNode::partition
     */
    struct AxisConstants
    {
        float cbmin;
        float cbmax;
        float cbdiff;
        float k;

        AxisConstants(float cbmin_, float cbmax_, float cbdiff_, float k_) : cbmin(cbmin_), cbmax(cbmax_),
                                                                             cbdiff(cbdiff_), k(k_) {}

        AxisConstants() {}
    };
};

/*
 * @brief A Bounding Volume Hierarchy
 */
class BVH : public TreeAccelBase<BVHNode>
{
    template<typename T> friend
    class TreeVisualizer;

public:
    BVH(const int leaf_threshold = 1, const int max_depth = std::numeric_limits<int>::max());
    virtual const char *printName() const override final;
    virtual void addBaseDataStructure(KIRK::CPU::Scene *scene) override final;
    virtual bool closestIntersection(KIRK::Intersection *hit) override final;
    virtual bool isIntersection(KIRK::Ray *ray, float tMax) override final;
};
}}

#endif
