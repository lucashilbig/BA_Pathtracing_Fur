#ifndef KIRK_MATERIAL_H
#define KIRK_MATERIAL_H

#include <glm/glm.hpp>
#include <string>
#include <memory>

#include "Color.h"
#include "Texture.h"

#include "Ray.h"
#include "KIRK/Utils/Gui/Gui.h"
#include "Shading/Shader.h"
#include "Shading/Bsdf.h"

namespace KIRK {

/**
A collection of possible texture/color map types.
*/
	enum class MatParamType
	{
		DIFFUSE,
		SPECULAR,
		NORMAL,
		VOLUME,
		EMISSION,
		TRANSPARENCY,
		REFLECTIVITY,
		ROUGHNESS,
		BUMP,
		ALPHA_LOBESHIFT,//longditudinal shift of a hair lobe. Between -10 and -5 degrees (marschner hairmodel property)
		BETA_LOBEWIDTH//longditudinal width (standart derivation) of a hair lobe. Between 5 and 10 degrees (marschner hairmodel property)
	};

	/**
	Small unified container for storing a value and a corresponding texture.
	*/
	template<typename BaseType>
	struct MatParam
	{
		BaseType value;

		//For checking if there is a texture, just do if(texture){...}
		std::shared_ptr<Texture> texture;
	};

	// Available MatParam types are...
	using MatParamFloat = MatParam<float>;
	using MatParamColor = MatParam<Color::RGBA>;

	/**
	A material is just a collection of values defining how an object is being shaded. No (enhanced) calculations should be found here.
	*/
	class Material : public KIRK::GuiElement
	{
	public:

		/**
		@param name The unique material name.
		*/
		Material(std::string name);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////
		//////	MATERIAL PARAMETERS.
		//////	For each of the following parameters, there should be a specialization of the 
		//////	fetchParameterColor or fetchParameterFloat method.
		//////
		MatParamColor m_diffuse = { Color::WHITE };				//!< Diffuse color of the material.
		MatParamColor m_specular = { Color::WHITE };				//!< Specular color of the reflection
		MatParamColor m_normal = { Color::GREEN };				//!< Normal map. Color::GREEN is equivalent to the up-vector (0,1,0)
		MatParamColor m_volume = { Color::WHITE };				//!< Color which will be used for refractions.
		MatParamColor m_emission = { Color::BLACK };				//!< Emissive color. Only used in Pathtracers.

		MatParamFloat m_transparency = { 0.f };					//!< Material transparency. (0 is fully opaque, 1 fully transparent)
		MatParamFloat m_reflectivity = { 0.f };					//!< Material reflectivity. Nothing to say here.
		MatParamFloat m_roughness = { 1.f };					//!< Material roughness. Used for reflections and refractions.
		MatParamFloat m_bump = { 0.f };					//!< Bump map and it's strength.

		//Marschner hair properties
		MatParamFloat m_alpha_shift = { 0.f };	//!< longditudinal shift of a hair lobe. Between -10 and -5 degrees. Stored in radians.
		MatParamFloat m_beta_width = { 0.f };	//!< longditudinal width (standart derivation) of a hair lobe. Between 5 and 10 degrees. Stored in radians.
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		std::string name = "";		//!< Material name. This one is unique!
		float m_ior = 1.52f;		//!< Index of refraction for this Material.

        std::shared_ptr<KIRK::Shader> m_shader;            //!< Shader
        std::shared_ptr<KIRK::BSDF> m_bsdf;			//!< BSDF for pathtracing.

		//ImGui parameters
		int m_current_bsdf = 0;

		/**
		* Fetches a color from a MatParamColor object. If it has a texture, the texcoord will be used to get
		* the color from the texture, otherwise it returns just the base value.
		*
		* CAUTION: Has to be specialized for every MatParamColor there is in the Material!
		*/
		template<MatParamType ParamType>
		Color::RGBA fetchParameterColor(const glm::vec2 &texcoord) const
		{
			throw std::invalid_argument("This parameter type is not valid or has no template implementation yet.");
		}

		/**
		* Fetches a float value from a MatParamFloat object. If it has a texture, the texcoord will be used to get
		* the color from the texture and the pixel value color vector length will be used as the value, otherwise
		* it returns just the base value.
		*
		* CAUTION: Has to be specialized for every MatParamFloat there is in the Material!
		*/
		template<MatParamType ParamType>
		float fetchParameterFloat(const glm::vec2 &texcoord) const
		{
			throw std::invalid_argument("This parameter type is not valid or has no template implementation yet.");
		}

		/**
		@return The material shinyness as 1/roughness (with a little safety epsilon. We don't want to divide by 0)
		*/
		float shininess()
		{
			return 1 / glm::max(m_roughness.value, 1e-3f);
		}

		/**
		CVK::Controllable override
		*/
		void onGui() override;

		/**
		Output stream operator for simplifying the material debugging.
		*/
		friend std::ostream& operator<<(std::ostream& os, const Material* material);

	private:

		/**
		Determines, if the MatParamColor has a texture. If so, the pixel color at the texture coordinate will be returned, else, only the base color.
		@param colorParam A MatParamColor object.
		@param texcoord A vec2 determining the texture coordinate.
		@return The pixel color at texcoord if there is a texture, otherwise just the colorParam base color.
		*/
		Color::RGBA getFromParam(const MatParamColor &colorParam, const glm::vec2 &texcoord) const;

		/**
		Determines, if the MatParamFloat has a texture. If so, the pixel color at the texture coordinate will be returned, else, only the base color.
		@param floatParam A MatParamFloat object.
		@param texcoord A vec2 determining the texture coordinate.
		@return The pixel color vector length at texcoord if there is a texture, otherwise just the floatParam base value.
		*/
		float getFromParam(const MatParamFloat &floatParam, const glm::vec2 &texcoord) const;

	};
	// DIFFUSE COLOR
	template<>
	inline Color::RGBA Material::fetchParameterColor<MatParamType::DIFFUSE>(const glm::vec2 &texcoord) const
	{
		return getFromParam(m_diffuse, texcoord);
	}

	// SPECULAR COLOR
	template<>
	inline Color::RGBA Material::fetchParameterColor<MatParamType::SPECULAR>(const glm::vec2 &texcoord) const
	{
		return getFromParam(m_specular, texcoord);
	}

	// NORMAL MAP
	template<>
	inline Color::RGBA Material::fetchParameterColor<MatParamType::NORMAL>(const glm::vec2 &texcoord) const
	{
		return getFromParam(m_normal, texcoord);
	}

	// VOLUME COLOR
	template<>
	inline Color::RGBA Material::fetchParameterColor<MatParamType::VOLUME>(const glm::vec2 &texcoord) const
	{
		return getFromParam(m_volume, texcoord);
	}

	// EMISSION COLOR
	template<>
	inline Color::RGBA Material::fetchParameterColor<MatParamType::EMISSION>(const glm::vec2 &texcoord) const
	{
		return getFromParam(m_emission, texcoord);
	}

	// TRANSPARENCY VALUE
	template<>
	inline float Material::fetchParameterFloat<MatParamType::TRANSPARENCY>(const glm::vec2 &texcoord) const
	{
		return getFromParam(m_transparency, texcoord);
	}

	//REFLECTIVITY VALUE
	template<>
	inline float Material::fetchParameterFloat<MatParamType::REFLECTIVITY>(const glm::vec2 &texcoord) const
	{
		return getFromParam(m_reflectivity, texcoord);
	}

	//ROUGHNESS VALUE
	template<>
	inline float Material::fetchParameterFloat<MatParamType::ROUGHNESS>(const glm::vec2 &texcoord) const
	{
		return getFromParam(m_roughness, texcoord);
	}

	//BUMP MAP TEXTURE AND BUMP STRENGTH
	template<>
	inline float Material::fetchParameterFloat<MatParamType::BUMP>(const glm::vec2 &texcoord) const
	{
		return getFromParam(m_bump, texcoord);
	}	

	//ALPHA_LOBESHIFT value
	template<>
	inline float Material::fetchParameterFloat<MatParamType::ALPHA_LOBESHIFT>(const glm::vec2 &texcoord) const
	{
		return getFromParam(m_alpha_shift, texcoord);
	}

	//BETA_LOBEWIDTH value
	template<>
	inline float Material::fetchParameterFloat<MatParamType::BETA_LOBEWIDTH>(const glm::vec2 &texcoord) const
	{
		return getFromParam(m_beta_width, texcoord);
	}
}

#endif //KIRK_MATERIAL_H
