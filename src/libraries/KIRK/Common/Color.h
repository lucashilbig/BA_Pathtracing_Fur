#ifndef KIRK_COLOR_H
#define KIRK_COLOR_H

#include <glm/glm.hpp>

namespace KIRK
{
	namespace Color 
	{
		typedef glm::vec4 RGBA; //!< Giving the glm::vec4 color a more descriptive name and format. */

		//A bunch of predefined colors
		const RGBA RED = RGBA{ 1,0,0,1 };
		const RGBA GREEN = RGBA{ 0,1,0,1 };
		const RGBA BLUE = RGBA{ 0,0,1,1 };
		const RGBA YELLOW = RGBA{ 1,1,0,1 };
		const RGBA MAGENTA = RGBA{ 1,0,1,1 };
		const RGBA CYAN = RGBA{ 0,1,1,1 };

		const RGBA WHITE = RGBA{ 1,1,1,1 };
		const RGBA BLACK = RGBA{ 0,0,0,1 };
		const RGBA CLEAR = RGBA{ 0,0,0,0 };

		const RGBA SKY_BLUE = RGBA{ 0.7f, 0.9f, 1.f, 1.f };
		const RGBA BLOOD_RED = RGBA{ 138.0f/255.0f,7.0f/255.0f,7.0f/255.0f, 1.f };
        
		/**
		Interpolates with a linear interpolation between two colors.
		@param first Interpolation start color
		@param secont Interpolatino end color
		@param alpha Interpolation alpha (as progress: float value starting at 0 and going to 1)
		@return a color_rgba value with the given interpolation values.
		*/
        inline RGBA mix(const RGBA &first, const RGBA &second, const float alpha)
		{
            return (1-alpha)*first + alpha*second;
        }
	}
}

#endif //KIRK_COLOR_H