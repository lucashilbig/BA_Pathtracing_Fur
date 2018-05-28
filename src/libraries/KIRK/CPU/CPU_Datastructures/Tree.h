#ifndef __CVK_DS_TREE_H__
#define __CVK_DS_TREE_H__

#include <vector>

namespace KIRK {
namespace CPU {

template<class T>
class Tree
{
public:
    Tree() {}

    virtual ~Tree()
    {
        for(T *t : m_children)
            delete t;
    }

    bool hasChildren() { return !m_children.empty(); }

    bool isLeaf() { return m_children.empty(); }

    int numTotalChildren()
    {
        int result = 1;
        for(T *c : m_children)
            result += c->numTotalChildren();
        return result;
    }

    int numTotalLeafs()
    {
        if(isLeaf())
            return 1;
        int result = 0;
        for(T *c : m_children)
            result += c->numTotalLeafs();
        return result;
    }

    int numTotalBranches()
    {
        if(isLeaf())
            return 0;
        int result = 1;
        for(T *c : m_children)
            result += c->numTotalBranches();
        return result;
    }

    std::vector<T *> *getChildren() { return &m_children; }

    int getRecursionDepth() { return m_recursionDepth; }

protected:
    int m_recursionDepth;
    virtual void subdivide() = 0;

    void addChild(T *c) { m_children.push_back(c); }

    std::vector<T *> m_children;
};

}}

#endif //__CVK_DS_TREE_H__
