#include "KIRK/Utils/Log.h"
#include "TreeAccel.h"
#include "CPU_BVH.h"
#include "CPU_KD.h"

namespace KIRK {
namespace CPU {

TreeAccelProperties::TreeAccelProperties(int lt, int md, int d, int an, int al) :
        leaf_threshold(lt),
        max_depth(md),
        depth(d),
        amnt_nodes(an),
        amnt_leaves(al)
{
}

template<typename T>
TreeAccelBase<T>::TreeAccelBase(const int leaf_threshold, const int max_depth) :
        m_prop(TreeAccelProperties(leaf_threshold, max_depth, 0, 0, 0))
{

    if(leaf_threshold < 1)
    {
		LOG_ERROR("Leaf threshold is too small! (Choose 1 or higher; default: 1)");
        throw std::runtime_error("Leaf threshold is too small! (Choose 1 or higher; default: 1)");
    }

    if(max_depth < 1)
    {
		LOG_ERROR("Maximum depth cannot be less than 1!");
        throw std::runtime_error("Maximum depth cannot be less than 1!");
    }
}

template<typename T>
TreeAccelBase<T>::~TreeAccelBase()
{
}

template<typename T>
void TreeAccelBase<T>::printAdditionalDebugInfo() const
{
}

/*
 * @brief Template specialization for KDNode
 * @return total size of this tree (recurse in children nodes)
 */
template<>
int TreeAccelBase<KDNode>::getSizeInBytes()
{
    return sizeof(KDTree) + m_root.getSizeInBytes();
}

/*
 * @brief Template specialization for BVHNode
 * @return total size of this tree (recurse in children nodes)
 */
template<>
int TreeAccelBase<BVHNode>::getSizeInBytes()
{
    return sizeof(BVH) + m_root.getSizeInBytes();
}

template<typename T>
void TreeAccelBase<T>::printDebugInfo()
{
	LOG_DEBUG("% depth: %", printName(), m_prop.depth);
	LOG_DEBUG("Memory size: % bytes.", getSizeInBytes());
	LOG_DEBUG("# of nodes: %, # of leaves: %", m_prop.amnt_nodes, m_prop.amnt_leaves);
    printAdditionalDebugInfo();
}

template<typename T, int N>
TreeAccelNode<T, N>::TreeAccelNode()
{
    for(int n = 0; n < N; ++n)
    {
        m_children[n] = 0;
    }
}

template<typename T, int N>
TreeAccelNode<T, N>::~TreeAccelNode()
{
    for(int n = 0; n < N; ++n)
    {
        if(m_children[n])
        {
            delete m_children[n];
        }
    }
}

/*
 * @brief Recursively determines size of tree made up of nodes of this type
 * @return size of tree in bytes
 */
template<typename T, int N>
int TreeAccelNode<T, N>::getSizeInBytes() const
{
    int result = sizeof(T); // T is the type used for the nodes
    result += m_candidateList.size() * sizeof(KIRK::Triangle *);

    for(int n = 0; n < N; ++n)
    {
        if(m_children[n])
        {
            result += m_children[n]->getSizeInBytes();
        }
    }

    return result;
}

    template TreeAccelBase<BVHNode>::TreeAccelBase(const int, const int);
template TreeAccelBase<BVHNode>::~TreeAccelBase();
template void TreeAccelBase<BVHNode>::printDebugInfo();
template TreeAccelBase<KDNode>::TreeAccelBase(const int, const int);
template TreeAccelBase<KDNode>::~TreeAccelBase();
template void TreeAccelBase<KDNode>::printDebugInfo();
template TreeAccelNode<BVHNode, 2>::TreeAccelNode();
template TreeAccelNode<BVHNode, 2>::~TreeAccelNode();
template int TreeAccelNode<BVHNode, 2>::getSizeInBytes() const;
template TreeAccelNode<KDNode, 2>::TreeAccelNode();
template TreeAccelNode<KDNode, 2>::~TreeAccelNode();
template int TreeAccelNode<KDNode, 2>::getSizeInBytes() const;

}
}
