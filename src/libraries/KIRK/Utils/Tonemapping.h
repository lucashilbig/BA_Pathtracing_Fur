#pragma once
#include <glm/glm.hpp>
#include <vector>
#include <cmath>
#include "KIRK/Utils/Gui/Gui.h"

namespace KIRK {

class Tonemapper : public KIRK::GuiElement
{
public:
	static constexpr float epsilon = 1e-06;
	static constexpr float log05 = -0.693147; // log(0.5). log has e as it's base.
	static constexpr float log2 = 0.693147; // log(2)
	static constexpr float e = 2.718281828;

	Tonemapper() : GuiElement() {};
	~Tonemapper() {};
	void map(std::vector<glm::vec3>& image, int width, int height);

	void onGui() override;

	float m_exposure = 0.0f; // Scene exposure
	float m_biasParam = 0.85f; // Bias Param
	float m_gammaval = 1.0f; // Gamma value should use 2.2
	float m_contParam = 0.0f; // Contrast improvement
	float m_white = 1.0f; // Max lum value of display - White 0.98
	float m_black = 0.0f; // Min lum value of display - Black
	float m_kernel_multiplier = 0.125; // Kernel size multiplier for center-weighted algorithm
	int m_center_x = -1;
	int m_center_y = -1;
	bool m_use_rec_gamma = false; // Whether to use Rec. 709 Gamma calculation
	bool m_center_weight = false; // Whether to calculate world luminance from the center of the image using gauss filtering

	glm::vec3 RGB2Yxy[3] = {{0.5141364, 0.3238786, 0.16036376},
		{0.265068, 0.67023428, 0.06409157},
		{0.0241188, 0.1228178, 0.84442666}
	};
	glm::vec3 Yxy2RGB[3] = {{2.5651, -1.1665, -0.3986},
		{-1.0217, 1.9777, 0.0439},
		{0.0753, -0.2543, 1.1892}
	};
private:
	float gui_world_luminance = 0;
	float gui_max_luminance = 0;

	void RGB_to_Yxy(std::vector<glm::vec3>& image, float& maximum, float& world_lum) const;
	void Yxy_to_RGB(std::vector<glm::vec3>& image) const;
	void tonemapping(std::vector<glm::vec3>& image, float exposure, float max_lum, float world_lum) const;

	static float bias(double b, double x)
	{
		return pow(x, b); //pow(x, log(b) / log(0.5));
	}

	void gamma_calc(std::vector<glm::vec3>& image) const;

	// Rec. 709 Gamma calculation
	void rec_gamma_calc(std::vector<glm::vec3>& image) const;

	void luminance_from_center(std::vector<glm::vec3>& image, int width, int height, int center_x, int center_y, float& world_lum);
};
}
