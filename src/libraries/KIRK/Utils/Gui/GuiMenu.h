#ifndef __KIRK_GUI_MENU_H
#define __KIRK_GUI_MENU_H

#include <vector>
#include <functional>

#include "Gui.h"

namespace KIRK {

	using MenuItemCallback = std::function<void()>;

	class MenuItem
	{
	public:
		virtual void draw() = 0;
	};

	class MenuButton : public MenuItem
	{
	public:
		MenuButton(const std::string &title, MenuItemCallback onclick, bool enabled = true);
		void draw() override;

		void setEnabled(bool enabled);

	private:
		std::string m_title;
		MenuItemCallback m_onclick_callback;
		bool m_enabled;
	};

	class MenuSeparator : public MenuItem
	{
	public:
		void draw() override;
	};

	class Menu
	{
	public:
		Menu(const std::string &title);
		void draw();

		template<typename T = MenuItem, typename... Args>
		void makeItem(Args... args)
		{
			m_menu_items.push_back(std::make_unique<T>(args...));
		}

		const std::string &getTitle() const
		{
			return m_title;
		}

	private:
		std::string m_title;
		std::vector<std::unique_ptr<MenuItem>> m_menu_items;
	};

	class GuiMainMenuBar : public GuiNamedElement
	{
	public:
		GuiMainMenuBar();
		~GuiMainMenuBar();

		void onGui() override;

		void addItem(const std::string &menu_title, const std::string &item, MenuItemCallback onclick)
		{
			for(auto &&menu : m_menus)
			{
				if(menu->getTitle() == menu_title)
				{
					menu->makeItem<MenuButton>(item, onclick);
					return;
				}
			}

			m_menus.push_back(std::make_unique<Menu>(menu_title));
			m_menus.back()->makeItem<MenuButton>(item, onclick);
		}

		void addSeparator(const std::string &menu_title)
		{
			for (auto &&menu : m_menus)
			{
				if (menu->getTitle() == menu_title)
				{
					menu->makeItem<MenuSeparator>();
					return;
				}
			}

			m_menus.push_back(std::make_unique<Menu>(menu_title));
			m_menus.back()->makeItem<MenuSeparator>();
		}

		void enableViewMenu(const std::string &title)
		{
			m_view_title = title;
			m_enable_view = true;
		}

	private:
		bool m_enable_view = false;
		std::string m_view_title;
		std::vector<std::unique_ptr<Menu>> m_menus;
	};
}

#endif // !__KIRK_GUI_MENU_H
