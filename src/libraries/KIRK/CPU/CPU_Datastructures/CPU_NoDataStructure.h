#ifndef __CVK_DS_NODATASTRUCTURE_H
#define __CVK_DS_NODATASTRUCTURE_H

#include "CPU_DataStructure.h"

namespace KIRK {
namespace CPU {

class NoDataStructure : public CPU_DataStructure
{
public:
    NoDataStructure();
    ~NoDataStructure();

    virtual void addBaseDataStructure(KIRK::CPU::Scene *scene);

    bool closestIntersection(KIRK::Intersection *hit);
    bool isIntersection(KIRK::Ray *ray, float tMax);

};
}}

#endif
