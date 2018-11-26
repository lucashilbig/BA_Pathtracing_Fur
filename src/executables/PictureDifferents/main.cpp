#include "KIRK/Common/Texture.h"



////////////////////////////////////////////////////////////////////////////////////////////
//
//  MAIN:
//  ARGUMENTS:
//		-s <scene_path>		Load a scene from file system with relative path  (Prioritized)
//		-r <scene_res>		Load a scene from RESOURCES_PATH
//		-w <img_width_px>	Set a render image width
//		-h <img_height_px>	Set a render image height
//
////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
	//textures we use for rendering the images
	KIRK::Texture input_1("C:/Users/LucasHilbig/Documents/Uni/Bachelorarbeit/Ressourcen/Bilder_Ausarbeitung/mar_250s_cam1.png");
	KIRK::Texture input_2("C:/Users/LucasHilbig/Documents/Uni/Bachelorarbeit/Ressourcen/Bilder_Ausarbeitung/deon_250s_cam1.png");

	glm::vec2 in1_size = input_1.getSize();
	glm::vec2 in2_size = input_2.getSize();

	if (in1_size.x != in2_size.x || in1_size.y != in2_size.y)
		return 0;

	int width = in1_size.x;
	int height = in1_size.y;

	KIRK::Texture output(input_2.getSize().x, input_2.getSize().y, 1);

	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			if (i == 500)
				float test;
			KIRK::Color::RGBA color = glm::vec4(100) * glm::abs(input_1.getColor(i, j) - input_2.getColor(i, j));
			output.setPixel(i, j, color);
			
		}
	}


	output.glDownload();
	output.saveTo(SCREENSHOTS_PATH "/diffPic.png");

	return 0;
}
