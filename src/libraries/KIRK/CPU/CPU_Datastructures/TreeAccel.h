#ifndef __CVK_DS_TREEACCEL_H__
#define __CVK_DS_TREEACCEL_H__

#include "BoundingBox.h"
#include "CPU_DataStructure.h"

#define LEFT_CHILD m_children[ 0 ]
#define RIGHT_CHILD m_children[ 1 ]

namespace KIRK {
namespace CPU {

/*
 * @brief A base for nodes in raytracing acceleration tree datastructures
 */
template<typename T, int N>
class TreeAccelNode : public Container
{
public:

	/*
	* @brief Destructor recursively destroys children
	*/
    ~TreeAccelNode();
    int getSizeInBytes() const;

    T **getChildren() { return m_children; }
    T *m_children[N];

protected:

	/*
	* @brief A base for nodes in raytracing acceleration tree datastructures
	*/
    TreeAccelNode();

};

/*
* @brief Properties applying to the entire datastructure, filled during construction
* @param lt amount of objects required to stop splitting nodes
* @param md maximum allowed total depth of tree
* @param d actual depth of the tree
* @param an total amount of nodes
* @param al total amount of leaves
*/
struct TreeAccelProperties
{
    int leaf_threshold;
    int max_depth;
    int depth;
    int amnt_nodes;
    int amnt_leaves;
    TreeAccelProperties(int lt, int md, int d, int an, int al);
};

/*
* @brief A base for raytracing acceleration tree datastructures
* @param leaf_threshold amount of objects required to stop splitting nodes
* @param max_depth maximum allowed total depth of tree
*/
template<typename T>
class TreeAccelBase : public CPU_DataStructure
{
    template<typename U> friend
    class TreeVisualizer;

public:

	/*
	* @brief Virtual destructor
	*/
    virtual ~TreeAccelBase();
    virtual const char *printName() const = 0;
    virtual void addBaseDataStructure(KIRK::CPU::Scene *scene) override = 0;
    virtual bool closestIntersection(KIRK::Intersection *hit) override = 0;
    virtual bool isIntersection(KIRK::Ray *ray, float tMax) override = 0;
    virtual void printAdditionalDebugInfo() const;
    virtual int getSizeInBytes() override final;

	/*
	* @brief Print some debugging information in console. This will usually be called right after construction
	*/
    virtual void printDebugInfo() override final;

    T &getRootNode() { return m_root; }

    int getDepth() { return m_prop.depth; }
    
    T m_root;
    TreeAccelProperties m_prop;
protected:
    TreeAccelBase(const int leaf_threshold, const int max_depth);

};
}}

#endif
