#ifndef __CVK_PT_LIGHTBOUNCE_H
#define __CVK_PT_LIGHTBOUNCE_H

//#include "KIRK/CPU/CPU_Datastructures/CPU_DataStructure.h"
#include <glm/glm.hpp>

namespace KIRK {

class Intersection;

namespace CPU {

/**
* \brief Holds the intersection of the LightBounce and the Importance
*/
struct LightBounce
{
	Intersection* hit;
	glm::vec3 importance;
	float distance; //!< Distance from light source
	float att_lin; //!< linear attenuation
	float att_quad; //!< quadratic attenuation
	int mat_flags;
};
}
}

#endif
