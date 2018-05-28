#include "GuiFileDialog.h"
#include "../Log.h"
#include "GuiTextDialog.h"

namespace KIRK
{
	GuiFileDialog::GuiFileDialog(const std::string &title, const std::string &positive_action, const std::string &negative_action, const std::string &start_path, bool check_if_exists, dialog_submit_callback on_submit)
		: GuiNamedElement(title), m_title(title), m_submit_callback(on_submit), m_positive(positive_action), m_negative(negative_action), m_check_if_exists(check_if_exists)
	{
		m_current_path = start_path;
	}

	GuiFileDialog::~GuiFileDialog()
	{
		
	}

	void GuiFileDialog::onGui()
	{
		if (!m_opened) {
			ImGui::OpenPopup(m_title.c_str());
			m_opened = true;
		}

		if (ImGui::BeginPopupModal(m_title.c_str(), nullptr, ImGuiWindowFlags_AlwaysAutoResize))// ImGui::Begin(m_title.c_str(), nullptr, ImVec2(400, 400), -1, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse))
		{
			ImVec2 button_small = ImVec2(ImGui::GetContentRegionAvailWidth(), 24);
			ImVec2 button_large = ImVec2(ImGui::GetContentRegionAvailWidth()/2-4, 32);
			char* buff = &m_filename[0];
			ImGui::PushItemWidth(ImGui::GetContentRegionAvailWidth());

			ImGui::TextColored(getGui()->Defaults.color_primary, m_current_path.string().c_str());

			ImGui::BeginChild("FileSystem", ImVec2(ImGui::GetContentRegionAvailWidth(), 400 - 114));

			ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(1, 1, 1, 0));
			if (m_current_path.has_parent_path() && ImGui::Button("..", button_small))
			{
				m_filename.clear();
				m_current_path = m_current_path.parent_path();
			}
			ImGui::PopStyleColor();

			m_file_exists = false;

			for(auto &c_path : fs::directory_iterator(m_current_path))
			{

				fs::path c(c_path);

				if (m_filename == c.filename().string())
				{
					ImGui::PushStyleColor(ImGuiCol_Button, getGui()->Defaults.color_accent);
				}
				else {
					ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(1, 1, 1, 0));
				}

				ImGui::PushStyleColor(ImGuiCol_ButtonHovered, getGui()->Defaults.color_accent_light);
				ImGui::PushStyleColor(ImGuiCol_ButtonActive, getGui()->Defaults.color_accent);
				if (fs::is_directory(c))
				{
					ImGui::PushStyleColor(ImGuiCol_ButtonHovered, getGui()->Defaults.color_primary_light);
					ImGui::PushStyleColor(ImGuiCol_ButtonActive, getGui()->Defaults.color_primary);
					if (ImGui::Button(c.filename().string().c_str(), button_small)) {
						m_filename.clear();
						m_current_path = c_path;
					}
					ImGui::PopStyleColor();
					ImGui::PopStyleColor();
				} 
				else if(c.extension() == ".json")
				{
					if (c.filename().string() == m_filename)
						m_file_exists = true;
					if (ImGui::Button(c.filename().string().c_str(), button_small))
					{
						if (m_filename != c.filename().string())
						{
							m_filename = c.filename().string();
						}
						else
						{
							if (m_submit_callback) m_submit_callback(*this, c);
							ImGui::CloseCurrentPopup();
							hintDeletion();
						}
					}
				}
				else if (c.extension() == ".obj")
				{
					if (c.filename().string() == m_filename)
						m_file_exists = true;
					if (ImGui::Button(c.filename().string().c_str(), button_small))
					{
						if (m_filename != c.filename().string())
						{
							m_filename = c.filename().string();
						}
						else
						{
							if (m_submit_callback) m_submit_callback(*this, c);
							ImGui::CloseCurrentPopup();
							hintDeletion();
						}
					}
				}
				ImGui::PopStyleColor();
				ImGui::PopStyleColor();
				ImGui::PopStyleColor();
			}

			ImGui::EndChild();

			ImGui::InputText("", buff, 128);
			if (ImGui::Button(m_negative.c_str(), button_large))
			{
				ImGui::CloseCurrentPopup();
				hintDeletion();
			}
			ImGui::SameLine();

			if (!(!m_check_if_exists || m_file_exists) || buff[0] == '\0') {
				ImGui::PushStyleColor(ImGuiCol_Button, getGui()->Defaults.color_frames);
				ImGui::PushStyleColor(ImGuiCol_ButtonHovered, getGui()->Defaults.color_frames);
				ImGui::PushStyleColor(ImGuiCol_ButtonActive, getGui()->Defaults.color_frames);
				ImGui::Button(m_positive.c_str(), button_large);
				ImGui::PopStyleColor();
				ImGui::PopStyleColor();
				ImGui::PopStyleColor();
			} else if (ImGui::Button(m_positive.c_str(), button_large))
			{
				if (m_submit_callback) m_submit_callback(*this, m_current_path.append(buff));
				ImGui::CloseCurrentPopup();
				hintDeletion();
			}
			ImGui::PopItemWidth();

			ImGui::EndPopup();
		}
	}
}