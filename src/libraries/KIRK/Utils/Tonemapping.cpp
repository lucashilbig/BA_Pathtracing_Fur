#include "Tonemapping.h"
#include <algorithm>
#include <numeric>

/* PAPER: http://resources.mpi-inf.mpg.de/tmo/logmap/logmap.pdf */

namespace KIRK {

void Tonemapper::map(std::vector<glm::vec3>& image, int width, int height)
{
	if(m_center_x < 0)
	{
		m_center_x = width / 2;
	}
	if(m_center_y < 0)
	{
		m_center_y = height / 2;
	}
	float max_luminance;
	float world_luminance;
	float exposure = pow(2, m_exposure); //default exposure is 1, 2^0
	RGB_to_Yxy(image, max_luminance, world_luminance);
	if(m_center_weight)
	{
		luminance_from_center(image, width, height, m_center_x, m_center_y, world_luminance);
	}
	gui_world_luminance = world_luminance;
	gui_max_luminance = max_luminance;
	tonemapping(image, exposure, max_luminance, world_luminance);
	Yxy_to_RGB(image);
	if(m_gammaval != 1.)
	{
		if(m_use_rec_gamma == true)
			rec_gamma_calc(image);
		else
			gamma_calc(image);
	}
	if(m_white != 1.0f || m_black != 0.0f)
	{
		for(auto& vec : image)
		{
			vec = glm::clamp(vec, m_black, m_white);
		}
	}
}

void Tonemapper::onGui()
{
	if (ImGui::CollapsingHeader("Tonemapper")) {
		ImGui::Text("Log. World Luminance: %.6f", gui_world_luminance);
		ImGui::Text("Max Luminance: %.6f", gui_max_luminance);
		ImGui::DragFloat("Exposure", &m_exposure, 1e-2f, -10.0f, 10.0f, "%.2f");
		ImGui::DragFloat("Bias", &m_biasParam, 1e-2f, 0.0f, 1.0f, "%.2f");
		ImGui::DragFloat("Gamma", &m_gammaval, 1e-2f, 0.0f, 4.0f, "%.2f");
		ImGui::DragFloat("Contrast", &m_contParam, 1e-2f, 0.0f, 100.0f, "%.2f");
		ImGui::DragFloat("White", &m_white, 1e-2f, 0.0f, 1.0f, "%.2f");
		ImGui::DragFloat("Black", &m_black, 1e-2f, 0.0f, 1.0f, "%.2f");
		ImGui::Checkbox("Rec. 709 Gamma", &m_use_rec_gamma);
		ImGui::Checkbox("Centered world luminance", &m_center_weight);
		ImGui::DragFloat("Kernel Multiplier", &m_kernel_multiplier, 1e-3f, 0.0f, 1.0f);
		ImGui::DragInt("Center X", &m_center_x);
		ImGui::DragInt("Center Y", &m_center_y);
	}
}

void Tonemapper::RGB_to_Yxy(std::vector<glm::vec3>& image, float& maximum, float& world_lum) const
{
	float W;
	float max;
	float sum = 0.0;
	glm::vec3 result;

	max = epsilon;
	for(auto& vec : image)
	{
		result = {glm::dot(RGB2Yxy[0], vec), glm::dot(RGB2Yxy[1], vec), glm::dot(RGB2Yxy[2], vec)};

		if((W = dot(result, {1, 1, 1})) > 0.0f)
		{
			//         Y         x             y
			vec = {result.y, result.x / W, result.y / W};
		}
		else
			vec = {0, 0, 0};

		max = glm::max(max, vec.x); // Max Luminance in Scene
		sum += log(2.3e-5 + vec.x); //Contrast constant (Tumblin paper)
	}
	maximum = max;
	world_lum = sum / image.size(); // Average logarithmic luminance
}

void Tonemapper::Yxy_to_RGB(std::vector<glm::vec3>& image) const
{
	glm::vec3 result;
	float X, Y, Z;

	for(auto& vec : image)
	{
		Y = vec.x; // Y 
		result.y = vec.y; // x 
		result.z = vec.z; // y 
		if(Y > epsilon && result.y > epsilon && result.z > epsilon)
		{
			X = result.y * Y / result.z;
			Z = X / result.y - X - Y;
		}
		else
			X = Z = epsilon;
		vec = {X, Y, Z};
		result = {glm::dot(Yxy2RGB[0], vec), glm::dot(Yxy2RGB[1], vec), glm::dot(Yxy2RGB[2], vec)};
		vec = result;
	}
}

void Tonemapper::tonemapping(std::vector<glm::vec3>& image, float exposure, float max_lum, float world_lum) const
{
	float Lmax, divider, av_lum, interpol, biasP, contP, exp_adapt;

	exp_adapt = 1; // pow(m_biasParam, 5);
	av_lum = exp(world_lum) / exp_adapt;

	biasP = log(m_biasParam) / log05;
	contP = 1 / m_contParam;
	Lmax = max_lum / av_lum;
	divider = log10(Lmax + 1);

	// Normal tone mapping of every pixel
	for(auto& vec : image)
	{
		// inverse gamma function to enhance contrast
		// Not in paper
		if(m_contParam)
		{
			vec.x = pow(vec.x, contP);
		}
		vec.x /= av_lum;
		vec.x *= exposure;

		interpol = log(2 + bias(biasP, vec.x / Lmax) * 8);

		vec.x = log(vec.x + 1) / interpol / divider;
	}
}

void Tonemapper::gamma_calc(std::vector<glm::vec3>& image) const
{
	float inv_gamma = 1 / m_gammaval;

	for(auto& vec : image)
	{
		vec = glm::pow(vec, glm::vec3(inv_gamma));
	}
}

void Tonemapper::rec_gamma_calc(std::vector<glm::vec3>& image) const
{
	float inv_gamma = 0.45 / m_gammaval * 2;
	float slope = 4.5;
	float start = 0.018;

	if(m_gammaval >= 2.1)
	{
		start = 0.018 / ((m_gammaval - 2) * 7.5);
		slope = 4.5 * ((m_gammaval - 2) * 7.5);
	}
	else if(m_gammaval <= 1.9)
	{
		start = 0.018 * ((2 - m_gammaval) * 7.5);
		slope = 4.5 / ((2 - m_gammaval) * 7.5);
	}

	for(auto& vec : image)
	{
		vec.x = vec.x <= start ?
			        vec.x * slope : 1.099 * pow(vec.x, inv_gamma) - 0.099;
		vec.y = vec.y <= start ?
			        vec.y * slope : 1.099 * pow(vec.y, inv_gamma) - 0.099;
		vec.z = vec.z <= start ?
			        vec.z * slope : 1.099 * pow(vec.z, inv_gamma) - 0.099;
	}
}

void Tonemapper::luminance_from_center(std::vector<glm::vec3>& image, int width, int height, int center_x, int center_y, float& world_lum)
{
	{
		// kernel size is determined by the smaller dimension of our image
		int kernel_size = width < height ? width * m_kernel_multiplier : height * m_kernel_multiplier;

		// clamp kernel_size in case the user tries out ridiculous kernel multipliers
		if(kernel_size > width || kernel_size > height)
		{
			kernel_size = glm::min(width, height);
		}
		else if(kernel_size < 1)
		{
			kernel_size = 1;
		}

		// we don't want to errorcheck here, as kernel_size is dynamic, so we ensure ourselves that the kernel_size is an odd number
		if(kernel_size % 2 == 0)
			kernel_size -= 1;

		int half_kernel_size = floor(kernel_size * 0.5);

		// we do not extend the image with a black border, because that would mess up our luminance
		// thus, we move it inside the image if the kernel would extend the image
		int x_start, y_start;
		if(center_x + half_kernel_size > width)
			x_start = width - kernel_size;
		else if(center_x - half_kernel_size < 0)
			x_start = 0;
		else
			x_start = center_x - half_kernel_size;

		if(center_y + half_kernel_size > height)
			y_start = height - kernel_size;
		else if(center_y - half_kernel_size < 0)
			y_start = 0;
		else
			y_start = center_y - half_kernel_size;

		// initialise our kernel
		std::vector<double> mask(kernel_size * kernel_size);

		/*
		* Build the Gaussian kernel
		*/
		// Playing around with the STL... this could just be a nested for-loop, but who cares?
		std::generate(mask.begin(), mask.end(), [idx = 0, kernel_size, half_kernel_size]() mutable
		              {
			              int x = idx % kernel_size - half_kernel_size;
			              int y = idx / kernel_size - half_kernel_size;
			              float r = sqrt(x * x + y * y);
			              ++idx;
			              return exp(-log(2) * pow(r / half_kernel_size, 2));
		              });
		double mean = kernel_size * kernel_size / std::accumulate(mask.begin(), mask.end(), 0.0);
		double sum = 0.0f;
		for(int x = x_start, i = 0; x < x_start + kernel_size; x++ , i++)
			for(int y = y_start, j = 0; y < y_start + kernel_size; y++ , j++)
			{
				int i1 = x * (y_start + kernel_size) + y;
				int i2 = j * kernel_size + i;

				sum += log(2.3e-5 + image[i1].x * mask[i2] * mean);
			}
		world_lum = float(sum / (kernel_size * kernel_size)); // Average logarithmic luminance
	}
}
}
