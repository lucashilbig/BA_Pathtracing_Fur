#ifndef __KIRK_GUI_H
#define __KIRK_GUI_H

#include <map>
#include <vector>
#include <functional>
#include <memory>

#include "KIRK/Utils/Window.h"
#include <externals/ImGui/imgui.h>
#include <externals/ImGui/imgui_impl_glfw_gl3.h>
#include "../../../externals/ImGui/imgui_impl_glfw_gl3.h"

namespace KIRK {

	class Callbacks
	{
	public:
		Callbacks() = delete;
		Callbacks(Callbacks &c) = delete;
		Callbacks(Callbacks &&c) = delete;
		~Callbacks() {};

		static void init(const KIRK::Window &cvk_window);

		static size_t addKeyCallback(std::function<void(GLFWwindow*, int, int, int, int)> callback);
		static size_t addCharCallback(std::function<void(GLFWwindow*, unsigned)> callback);
		static size_t addMouseButtonCallback(std::function<void(GLFWwindow*, int, int, int)> callback);
		static size_t addScrollCallback(std::function<void(GLFWwindow*, double, double)> callback);
		static size_t addFileDropCallback(std::function<void(GLFWwindow*, int, const char**)> callback);

	private:
		static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
		static void charCallback(GLFWwindow* window, unsigned c);
		static void mouseButtonCallback(GLFWwindow* window, int key, int action, int mods);
		static void scrollCallbacks(GLFWwindow* window, double x, double y);
		static void fileDropCallback(GLFWwindow* window, int count, const char** paths);

		static GLFWwindow* m_window;

		static std::vector<std::function<void(GLFWwindow*, unsigned)>> m_charCallbacks;
		static std::vector<std::function<void(GLFWwindow*, int, int, int, int)>> m_keyCallbacks;
		static std::vector<std::function<void(GLFWwindow*, int, int, int)>> m_mouseButtonCallbacks;
		static std::vector<std::function<void(GLFWwindow*, double, double)>> m_scrollCallbacks;
		static std::vector<std::function<void(GLFWwindow*, int, const char**)>> m_dropCallbacks;
		
	};

	class Gui;

class GuiElement
{
public:
	GuiElement() : m_enabled(true) {};
	virtual ~GuiElement() {};

	void draw()
	{
		if (m_enabled)
			onGui();
	}

	void setEnabled(bool enabled)
	{
		m_enabled = enabled;
	}

	bool isEnabled() const
	{
		return m_enabled;
	}

	void hintDeletion()
	{
		m_delete = true;
	}

	bool hintedDeletion() const
	{
		return m_delete;
	}

	virtual void onGui() = 0;

protected:
	bool m_enabled;
	bool m_delete = false;
};

class GuiNamedElement : public GuiElement
{
public:
	GuiNamedElement(const std::string &title);
	virtual ~GuiNamedElement() override;

	const std::string &getTitle() const;

	void setGui(Gui *gui)
	{
		m_gui = gui;
	}

	Gui* getGui()
	{
		return m_gui;
	}

protected:
	//Raw ptr because of unpredictability of Gui.
	Gui* m_gui;
	std::string m_title;
};

class GuiNamedLambdaElement : public GuiNamedElement
{
public:
	GuiNamedLambdaElement(const std::string &name, std::function<void()> on_gui);
	~GuiNamedLambdaElement();

	void onGui() override;

private:
	std::function<void()> m_on_gui;
};

class GuiWindow : public GuiNamedElement
{
public:
	GuiWindow(const std::string &title, bool start_open, bool always_open = false);
	~GuiWindow();

	unsigned add(GuiElement *element);

	/**
	TODO: CAUTION! Memory leak because never deleted.
	*/
	template<typename Type, typename... Args>
	unsigned make(Args... args)
	{
		auto ptr = new Type(args...);
		return add(ptr);
	}
	
	template<typename Type = GuiElement>
	Type *get(unsigned index)
	{
		if(index >= m_elements.size())
		{
            //LOG_WARN("Index out of range (id % of % elements).", index, m_elements.size());
			return nullptr;
		}

		return m_elements[index];
	}

	void detachChildren()
	{
		/*for (auto it = m_elements.begin(); it != m_elements.end(); ++it)
		{
			if(*it)
				delete *it;
		}*/
		m_elements.clear();
	}

	size_t getChildCount() const
	{
		return m_elements.size();
	}

	bool alwaysOpen() const
	{
		return m_always_open;
	}

	void onGui() override;

private:
	std::vector<GuiElement*> m_elements;
	bool m_always_open;
};

class Gui : public std::enable_shared_from_this<Gui>
{
public:
    Gui(const KIRK::Window &cvk_window);
	~Gui();

	template<typename Type>
	void add(std::shared_ptr<Type> named_element)
	{
		named_element->setGui(this);
		m_named_elements.emplace(named_element->getTitle(), std::dynamic_pointer_cast<GuiNamedElement>(named_element));
	}

	template<typename Type, typename... Args>
	std::shared_ptr<Type> make(Args... args)
	{
		auto ptr = std::make_shared<Type>(args...);
		ptr->setGui(this);
		add(ptr);
		return ptr;
	}

	void remove(const std::string &title);

	template<typename Type = GuiNamedElement>
	std::shared_ptr<Type> get(const std::string &title)
	{
		if (m_named_elements.count(title) != 0)
			return std::dynamic_pointer_cast<Type>(m_named_elements.at(title));
		else
			return std::shared_ptr<Type>(nullptr);
	}

	template<typename Type = GuiNamedElement>
	void forAll(std::function<void(std::shared_ptr<Type>)> do_for)
	{
		for (auto &&element : m_named_elements)
		{
			if (std::shared_ptr<Type> converted = std::dynamic_pointer_cast<Type>(element.second))
			{
				do_for(converted);
			}
		}
	}

    void render();
    bool wantsKeyboardInput();
    bool wantsMouseInput();
	
	struct
	{
		ImVec4 color_background;
		ImVec4 color_menu;
		ImVec4 color_bar;
		ImVec4 color_frames;

		ImVec4 color_primary;
		ImVec4 color_primary_light;

		ImVec4 color_accent;
		ImVec4 color_accent_light;

		ImVec4 color_modal_darkening;
	} Defaults;

	GLFWwindow* getWindow() const
	{
		return m_window;
	}

private:

	GLFWwindow* m_window;

	ImVec4 fromHex(unsigned hex)
	{
		return ImVec4(
			((hex >> 16) & 0xff) / 255.0,	//R
			((hex >> 8) & 0xff) / 255.0,	//G
			((hex >> 0) & 0xff) / 255.0,	//B
			((hex >> 24) & 0xff) / 255.0	//A
		);
	}
	
	std::map<std::string, std::shared_ptr<GuiNamedElement>> m_named_elements;
};

}
#endif // __KIRK_GUI_H
