#ifndef KIRK_SHADER_H
#define KIRK_SHADER_H

#include "KIRK/CPU/CPU_Raytracer/CPU_LightBounce.h"
#include <vector>

namespace KIRK {

	class Intersection;
	class Ray;
	namespace CPU {
		class PathTracer;
		class Bounce;
	}

	using LightBounceIterator = std::vector <KIRK::CPU::LightBounce>::iterator;

	class Shader
	{
	public:
        Shader(const std::string name) : m_name(name){}
		virtual ~Shader() {};
		virtual void shade(KIRK::CPU::PathTracer& pathtracer, const KIRK::Intersection& hit, KIRK::CPU::Bounce& resultBounce, KIRK::Ray& resultRay) {throw std::logic_error("This shader is not implemented for a pathtracer");}
        const std::string m_name;
    };
}

#endif
