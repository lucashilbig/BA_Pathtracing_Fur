#include "GuiMenu.h"

namespace KIRK
{
	MenuButton::MenuButton(const std::string &title, MenuItemCallback onclick, bool enabled)
		: m_title(title), m_onclick_callback(onclick), m_enabled(enabled)
	{}

	void MenuButton::draw()
	{
		if(ImGui::MenuItem(m_title.c_str(), nullptr, false, m_enabled))
		{
			m_onclick_callback();
		}
	}

	void MenuButton::setEnabled(bool enabled)
	{
		m_enabled = enabled;
	}

	void MenuSeparator::draw()
	{
		ImGui::Separator();
	}

	Menu::Menu(const std::string &title)
		: m_title(title)
	{}

	void Menu::draw()
	{
		if(ImGui::BeginMenu(m_title.c_str()))
		{
			for (auto &&item : m_menu_items)
				item->draw();
			ImGui::EndMenu();
		}
	}


	GuiMainMenuBar::GuiMainMenuBar()
		: GuiNamedElement("_main_menu")
	{
		
	}

	GuiMainMenuBar::~GuiMainMenuBar()
	{
		
	}

	void GuiMainMenuBar::onGui()
	{
		ImGui::BeginMainMenuBar();

		for (auto &&item : m_menus)
			item->draw();


		if (m_enable_view)
		{
			if (ImGui::BeginMenu(m_view_title.c_str())) {
				getGui()->forAll<GuiWindow>([](std::shared_ptr<KIRK::GuiWindow> window)
				{
					if (!window->alwaysOpen() && ImGui::MenuItem(window->getTitle().c_str()))
					{
						window->setEnabled(true);
					}
				});
				ImGui::EndMenu();
			}
		}

		ImGui::EndMainMenuBar();
	}
}