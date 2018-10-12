#include "Material.h"
#include "Shading/Bsdf.h"
#include "Shading/SimpleShader.h"
#include "Shading/MarschnerHairShader.h"
#include "Shading/ChiangHairShader.h"

namespace KIRK
{

    KIRK::Material::Material(std::string name) : KIRK::GuiElement(), name(name), m_bsdf(std::make_shared<BSDF>("LambertianReflectionBSDF", LambertianReflectionBSDF::localSample, LambertianReflectionBSDF::evaluateLight)), m_shader(new SimpleShader("default")), m_current_bsdf(0)
	{}

	Material::Material(std::string name, int bsdf) : KIRK::GuiElement(), name(name), m_bsdf(initBSDF(bsdf))
	{}

	KIRK::Color::RGBA KIRK::Material::getFromParam(const MatParamColor &colorParam, const glm::vec2 &texcoord) const
	{
		return colorParam.texture ? colorParam.texture->getColor(texcoord.x, texcoord.y) : colorParam.value;
	}

	float KIRK::Material::getFromParam(const MatParamFloat &floatParam, const glm::vec2 &texcoord) const
	{
		return floatParam.texture ? glm::length(floatParam.texture->getColor(texcoord.x, texcoord.y)) : floatParam.value;
	}

	std::shared_ptr<KIRK::BSDF> Material::initBSDF(int bsdf)
	{
		std::shared_ptr<KIRK::BSDF> bsdf_ptr;
		int safeBsdf = (bsdf < 0 || bsdf > 8) ? 0 : bsdf;
		switch (safeBsdf)//set our bsdf
		{
		case 0:
			bsdf_ptr = std::make_shared<BSDF>("LambertianReflectionBSDF",
				LambertianReflectionBSDF::localSample,
				LambertianReflectionBSDF::evaluateLight);
			break;
		case 1:
			bsdf_ptr = std::make_shared<BSDF>("SpecularReflectionBSDF", SpecularReflectionBSDF::localSample,
				SpecularReflectionBSDF::evaluateLight);
			break;
		case 2:
			bsdf_ptr = std::make_shared<BSDF>("GlassBSDF", GlassBSDF::localSample, GlassBSDF::evaluateLight);
			break;
		case 3:
			bsdf_ptr = std::make_shared<BSDF>("LambertianTransmissionBSDF",
				LambertianTransmissionBSDF::localSample,
				LambertianTransmissionBSDF::evaluateLight);
			break;
		case 4:
			bsdf_ptr = std::make_shared<BSDF>("EmissionBSDF", EmissionBSDF::localSample,
				EmissionBSDF::evaluateLight);
			break;
		case 5:
			bsdf_ptr = std::make_shared<BSDF>("TransparentBSDF", TransparentBSDF::localSample,
				TransparentBSDF::evaluateLight);
			break;
		case 6:
			bsdf_ptr = std::make_shared<BSDF>("MarschnerHairBSDF", MarschnerHairBSDF::localSample,
				MarschnerHairBSDF::evaluateLight);
			break;
		case 7:
			bsdf_ptr = std::make_shared<BSDF>("DEonHairBSDF", DEonHairBSDF::localSample,
				DEonHairBSDF::evaluateLight);
			break;
		case 8:
			bsdf_ptr = std::make_shared<BSDF>("ChiangHairBSDF", ChiangHairBSDF::localSample,
				ChiangHairBSDF::evaluateLight);
			break;
		default:
			bsdf_ptr = std::make_shared<BSDF>("LambertianReflectionBSDF",
				LambertianReflectionBSDF::localSample,
				LambertianReflectionBSDF::evaluateLight);
			break;
		}

		//set the shader which needs to be used for the selected bsdf
		if (safeBsdf >= 0 && safeBsdf < 6)
			m_shader = KIRK::ShaderFactory::getInstance().getShader("SimpleShader");
		else if (safeBsdf == 6 || safeBsdf == 7)
			m_shader = KIRK::ShaderFactory::getInstance().getShader("MarschnerHairShader");
		else if (safeBsdf == 8)
			m_shader = KIRK::ShaderFactory::getInstance().getShader("ChiangHairShader");

		m_current_bsdf = safeBsdf;
		//finaly return our bsdf
		return bsdf_ptr;
	}
		
	void imguicolor(const char* name, KIRK::Color::RGBA *color)
	{
		float dc[4] = { color->r, color->g, color->b, color->a };
		ImGui::ColorEdit4(name, dc);
		*color = { dc[0], dc[1], dc[2], dc[3] };
	}

	void KIRK::Material::onGui()
	{
		if (ImGui::TreeNode(name.c_str()))
        {
            ImGui::BeginGroup();

            ImGui::PushID(("mat_id_" + name).c_str());

            int oldcb = m_current_bsdf;
            const char *names[] = {"Diffuse", "Glossy", "Glass", "Translucent", "Emission", "Transparent", "MarschnerHair", "DEonHair", "ChiangHair" };
            ImGui::Combo("BSDF", &m_current_bsdf, names, 9);

            if(m_current_bsdf != oldcb)
            {
                switch(m_current_bsdf)
                {
                    case 0:
                        m_bsdf = initBSDF(0);
                        break;
                    case 1:
                        m_bsdf = initBSDF(1);
                        break;
                    case 2:
                        m_bsdf = initBSDF(2);
                        break;
                    case 3:
                        m_bsdf = initBSDF(3);
                        break;
                    case 4:
                        m_bsdf = initBSDF(4);
                        break;
                    case 5:
                        m_bsdf = initBSDF(5);
                        break;
					case 6:
						m_bsdf = initBSDF(6);
						break;
					case 7:
						m_bsdf = initBSDF(7);
						break;
					case 8:
						m_bsdf = initBSDF(8);
						break;
                    default:
						m_bsdf = initBSDF(0);
                        break;
                }
            }

            ImGui::SliderFloat("Reflectivity", &m_reflectivity.value, 0.f, 1.f);
            ImGui::SliderFloat("Roughness", &m_roughness.value, 0.f, 1.f);
            ImGui::SliderFloat("Transparency", &m_transparency.value, 0.f, 1.f);

            imguicolor("Diffuse Color", &m_diffuse.value);
            imguicolor("Reflection Color", &m_specular.value);
            imguicolor("Volume Color", &m_volume.value);
            imguicolor("Emission Color", &m_emission.value);

            ImGui::PopID();

            ImGui::EndGroup();
            ImGui::TreePop();
        }
	}

	std::ostream& operator<<(std::ostream& os, const KIRK::Material* material)
	{
		if (material)
		{
			os << "Printing Material: " << std::endl
				<< "Name: " << material->name << std::endl
				<< "-----" << std::endl

				<< "diffuse: ("
				<< material->m_diffuse.value.x << ", "
				<< material->m_diffuse.value.y << ", "
				<< material->m_diffuse.value.z << ")" << std::endl
				<< "-----" << std::endl

				<< "specular: ("
				<< material->m_specular.value.x << ", "
				<< material->m_specular.value.y << ", "
				<< material->m_specular.value.z << ")" << std::endl
				<< "-----" << std::endl

				<< "absortion: ("
				<< material->m_volume.value.x << ", "
				<< material->m_volume.value.y << ", "
				<< material->m_volume.value.z << ")" << std::endl
				<< "-----" << std::endl

				<< "emission: ("
				<< material->m_emission.value.x << ", "
				<< material->m_emission.value.y << ", "
				<< material->m_emission.value.z << ")" << std::endl
				<< "-----" << std::endl

				<< "ior: ("
				<< material->m_ior << ")" << std::endl
				<< "-----" << std::endl

				<< "transparency: ("
				<< material->m_transparency.value << ")" << std::endl
				<< "-----" << std::endl

				<< "reflectivity: ("
				<< material->m_reflectivity.value << ")" << std::endl
				<< "-----" << std::endl

				<< "roughness: ("
				<< material->m_roughness.value << ")" << std::endl
				<< "-----" << std::endl;

			return os << std::endl;
		}
		else
		{
			return os << "Well ... seems like there is no material :/" << std::endl;
		}
	}
	
}
