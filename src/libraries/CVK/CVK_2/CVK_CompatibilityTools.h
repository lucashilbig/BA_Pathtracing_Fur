#ifndef __CVK_COMPATIBILITYTOOLS_H
#define __CVK_COMPATIBILITYTOOLS_H

#include "CVK_Defs.h"


namespace CVK {
/**
 * Use OpenGL Core Profile if using Apple
 */
void useOpenGL33CoreProfile();
void useOpenGLCoreProfile(int major, int minor);
/**
 * check the compatibility of the hardware by giving useful output
 */
void checkCompatibility();
}

#endif /* __CVK_COMPATIBILITYTOOLS_H */
