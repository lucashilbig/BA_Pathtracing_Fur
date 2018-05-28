#ifndef __CVK_TEAPOT_H
#define __CVK_TEAPOT_H

#include "CVK_Defs.h"
#include "CVK_Geometry.h"

namespace CVK {

/**
* The Teapot is one of the example geometry classes in the CVK.
* @brief Teapot class as Geometry
*/
class Teapot : public CVK::Geometry
{
public:
    /**
     * Standard Constructor for the Teapot
     */
    Teapot();
    /**
     * Standard Destructor for the Teapot
     */
    ~Teapot();

private:
    /**
    * Create the Teapot and the buffers with the given attributes
    * @brief Create the Teapot and the buffers
    */
    void create();
};

};

#endif /* __CVK_TEAPOT_H */

