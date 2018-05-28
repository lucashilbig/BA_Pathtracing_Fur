#ifndef __CVK_RT_Templates
#define __CVK_RT_Templates

/*
* A class which provides template names.
* These names will be used as paramater T which has to be specified by the user.
* The given struct templates will return a string which represents the type of the parameter T
**/

namespace KIRK {

template<typename T>
struct CLTypes
{
    static const char *getName() { return "Unsupported"; }
};

template<>
struct CLTypes<int>
{
    static const char *getName() { return "int"; }
};

template<>
struct CLTypes<float>
{
    static const char *getName() { return "float"; }
};

template<>
struct CLTypes<unsigned>
{
    static const char *getName() { return "unsigned"; }
};

}
#endif /* __CVK_RT_Templates */
