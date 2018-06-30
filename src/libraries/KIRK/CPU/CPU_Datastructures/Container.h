#ifndef __CVK_DS_CONTAINER_H__
#define __CVK_DS_CONTAINER_H__


#include <vector>

#include "KIRK/Utils/Log.h"
#include "KIRK/Common/Ray.h"
#include "KIRK/Common/Object.h"
#include "KIRK/Common/Intersection.h"

namespace KIRK {
namespace CPU {

class Container
{
public:
    Container();
	virtual ~Container();

    bool closestIntersectionWithCandidates(KIRK::Intersection *hit, float tMin = 0.0f, float tMax = FLT_MAX);
    bool isIntersectionWithCandidates(KIRK::Ray *ray, float tMax = FLT_MAX);

    inline int getNumberOfCandidates() const { return m_candidateList.size(); }

    inline bool hasCandidates() const { return m_hasCandidates; }

    inline void addCandidate(KIRK::Triangle *o)
    {
        if(m_candidateList.empty())
            m_candidateList.push_back(o);
        else if(m_candidateList.back() != o)
            m_candidateList.push_back(o);
        m_hasCandidates = true;
    }

    inline std::vector<KIRK::Triangle *> *getCandidates() { return &m_candidateList; }

protected:
    std::vector<KIRK::Triangle *> m_candidateList;
    bool m_hasCandidates;

};

}}

#endif //__CVK_DS_CONTAINER_H__
