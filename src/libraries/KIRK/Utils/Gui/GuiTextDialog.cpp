#include "GuiTextDialog.h"

namespace KIRK
{
	GuiTextDialog::GuiTextDialog(const std::string &title, const std::string &message)
		: GuiNamedElement(title), m_content(message)
	{}

	void GuiTextDialog::onGui()
	{
		if (!m_opened) {
			ImGui::OpenPopup(m_title.c_str());
			m_opened = true;
		}

		ImGui::SetNextWindowContentWidth(400);
		if (ImGui::BeginPopupModal(m_title.c_str(), nullptr, ImGuiWindowFlags_AlwaysAutoResize))// ImGui::Begin(m_title.c_str(), nullptr, ImVec2(400, 400), -1, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse))
		{
			ImGui::TextWrapped(m_content.c_str());

			ImVec2 button_large = ImVec2(ImGui::GetContentRegionAvailWidth() / 2 - 4, 32);

			if (ImGui::Button(m_negative.c_str(), button_large))
			{
				if (m_onclick_negative) m_onclick_negative();
				ImGui::CloseCurrentPopup();
				hintDeletion();
			}
			ImGui::SameLine();

			if (ImGui::Button(m_positive.c_str(), button_large))
			{
				if (m_onclick_positive) m_onclick_positive();
				ImGui::CloseCurrentPopup();
				hintDeletion();
			}

			ImGui::EndPopup();
		}
	}

	GuiTextDialog &GuiTextDialog::addPositive(std::string positive, std::function<void()> onclick)
	{
		m_has_positive = true;
		m_positive = positive;
		m_onclick_positive = onclick;
		return *this;
	}

	GuiTextDialog &GuiTextDialog::addNegative(std::string negative, std::function<void()> onclick)
	{
		m_has_negative = true;
		m_negative = negative;
		m_onclick_negative = onclick;
		return *this;
	}
}