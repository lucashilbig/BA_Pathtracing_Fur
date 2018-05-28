#ifndef __SAMPLER_H
#define __SAMPLER_H

#include <random>
#include <glm/glm.hpp>

namespace KIRK
{
	template<typename Distribution = std::uniform_real_distribution<double>>
	class Sampler
	{
	public:
		Sampler();
		float sample1D();
		glm::vec2 sample2D();
	private:
		std::mt19937 m_gen;
		Distribution m_dis;
	};

	using UniformSampler = Sampler<>;

	template <typename Distribution>
	Sampler<Distribution>::Sampler()
	{
		std::random_device rd;
		m_gen = std::mt19937(rd());
	}

	template <typename Distribution>
	float Sampler<Distribution>::sample1D()
	{
		return m_dis(m_gen);
	}

	template <typename Distribution>
	glm::vec2 Sampler<Distribution>::sample2D()
	{
		return {sample1D(), sample1D()};
	}
}

#endif