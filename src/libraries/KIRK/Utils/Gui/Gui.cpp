#include "Gui.h"

namespace KIRK
{
	GLFWwindow* Callbacks::m_window;

	std::vector<std::function<void(GLFWwindow*, unsigned)>> Callbacks::m_charCallbacks;
	std::vector<std::function<void(GLFWwindow*, int, int, int, int)>> Callbacks::m_keyCallbacks;
	std::vector<std::function<void(GLFWwindow*, int, int, int)>> Callbacks::m_mouseButtonCallbacks;
	std::vector<std::function<void(GLFWwindow*, double, double)>> Callbacks::m_scrollCallbacks;
	std::vector<std::function<void(GLFWwindow*, int, const char**)>> Callbacks::m_dropCallbacks;

	void Callbacks::init(const KIRK::Window &cvk_window)
	{
		m_window = cvk_window.getGLFWWindow();
		glfwSetKeyCallback(m_window, Callbacks::keyCallback);
		glfwSetCharCallback(m_window, Callbacks::charCallback);
		glfwSetMouseButtonCallback(m_window, Callbacks::mouseButtonCallback);
		glfwSetScrollCallback(m_window, Callbacks::scrollCallbacks);
		glfwSetDropCallback(m_window, Callbacks::fileDropCallback);
	}

	size_t Callbacks::addKeyCallback(std::function<void(GLFWwindow*, int, int, int, int)> callback)
	{
		m_keyCallbacks.push_back(callback);
		return m_keyCallbacks.size() - 1;
	}

	size_t Callbacks::addCharCallback(std::function<void(GLFWwindow*, unsigned)> callback)
	{
		m_charCallbacks.push_back(callback);
		return m_charCallbacks.size() - 1;
	}

	size_t Callbacks::addMouseButtonCallback(std::function<void(GLFWwindow*, int, int, int)> callback)
	{
		m_mouseButtonCallbacks.push_back(callback);
		return m_mouseButtonCallbacks.size() - 1;
	}

	size_t Callbacks::addScrollCallback(std::function<void(GLFWwindow*, double, double)> callback)
	{
		m_scrollCallbacks.push_back(callback);
		return m_scrollCallbacks.size() - 1;
	}

	size_t Callbacks::addFileDropCallback(std::function<void(GLFWwindow*, int, const char**)> callback)
	{
		m_dropCallbacks.push_back(callback);
		return m_dropCallbacks.size() - 1;
	}

	void Callbacks::charCallback(GLFWwindow* window, unsigned c)
	{
		if(ImGui::GetIO().WantCaptureKeyboard)
		{
			ImGui_ImplGlfwGL3_CharCallback(window, c);
			return;
		}

		for (auto &&callback : m_charCallbacks)
			callback(window, c);
	}

	void Callbacks::keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
	{
		if (ImGui::GetIO().WantCaptureKeyboard)
		{
			ImGui_ImplGlfwGL3_KeyCallback(window, key, scancode, action, mods);
			return;
		}

		for (auto &&callback : m_keyCallbacks)
			callback(window, key, scancode, action, mods);
	}

	void Callbacks::mouseButtonCallback(GLFWwindow* window, int key, int action, int mods)
	{
		if(ImGui::GetIO().WantCaptureMouse)
		{
			ImGui_ImplGlfwGL3_MouseButtonCallback(window, key, action, mods);
			return;
		}

		for (auto &&callback : m_mouseButtonCallbacks)
			callback(window, key, action, mods);
	}

	void Callbacks::scrollCallbacks(GLFWwindow* window, double x, double y)
	{
		if (ImGui::GetIO().WantCaptureMouse)
		{
			ImGui_ImplGlfwGL3_ScrollCallback(window, x, y);
			return;
		}

		for (auto &&callback : m_scrollCallbacks)
			callback(window, x, y);
	}

	void Callbacks::fileDropCallback(GLFWwindow* window, int count, const char** files)
	{
		for (auto &&callback : m_dropCallbacks)
			callback(window, count, files);
	}


GuiNamedElement::GuiNamedElement(const std::string &title)
	: m_title(title)
{}

GuiNamedElement::~GuiNamedElement()
{}

const std::string &GuiNamedElement::getTitle() const
{
	return m_title;
}

GuiNamedLambdaElement::GuiNamedLambdaElement(const std::string &name, std::function<void()> on_gui)
	: GuiNamedElement(name), m_on_gui(on_gui)
{}

GuiNamedLambdaElement::~GuiNamedLambdaElement()
{
	
}

void GuiNamedLambdaElement::onGui()
{
	m_on_gui();
}

GuiWindow::GuiWindow(const std::string &title, bool start_open, bool always_open)
	: GuiNamedElement(title), m_always_open(always_open)
{
	m_enabled = start_open;
}

GuiWindow::~GuiWindow()
{}

void GuiWindow::onGui()
{
	bool always_open_lock = true;
	bool *opened = m_always_open ? &always_open_lock : &m_enabled;
	ImGui::Begin(getTitle().c_str(), opened);

	for (auto iterator = m_elements.begin(); iterator != m_elements.end();)
	{
		if (!(*iterator))
			iterator = m_elements.erase(iterator);
		else
		{
			(*iterator)->draw();
			++iterator;
		}
	}

	ImGui::End();
	
}

unsigned GuiWindow::add(GuiElement *element)
{
	m_elements.push_back(element);
	return unsigned(m_elements.size() - 1);
}

void Gui::remove(const std::string &title)
{
	if(m_named_elements.count(title) != 0)
		m_named_elements.erase(title);
}

Gui::Gui(const KIRK::Window &window)
{
	m_window = window.getGLFWWindow();
    // Init ImGui
    ImGui_ImplGlfwGL3_Init(window.getGLFWWindow(), true);

    // Load fonts
    ImGuiIO &io = ImGui::GetIO();
    //io.Fonts->AddFontDefault();
	io.Fonts->AddFontFromFileTTF(GUI_RESOURCES_PATH "/Montserrat-SemiBold.ttf", 14);

	Defaults.color_background = fromHex(0xff303030);
	Defaults.color_menu = fromHex(0xff000000);
	Defaults.color_bar = fromHex(0xff212121);
	Defaults.color_frames = fromHex(0xff424242);

	Defaults.color_primary = fromHex(0xff2196f3);
	Defaults.color_primary_light = fromHex(0xff90caf9);

	Defaults.color_accent = fromHex(0xffff5722);
	Defaults.color_accent_light = fromHex(0xffff8a65);

	Defaults.color_modal_darkening = fromHex(0x80212121);

	ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 2);
	ImGui::PushStyleVar(ImGuiStyleVar_FrameRounding, 2);
	ImGui::PushStyleVar(ImGuiStyleVar_ChildWindowRounding, 2);

	ImGui::PushStyleVar(ImGuiStyleVar_IndentSpacing, 8);
	ImGui::PushStyleVar(ImGuiStyleVar_ItemInnerSpacing, ImVec2(4, 4));
	ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(8, 8));
	ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(8, 8));

	ImGui::PushStyleColor(ImGuiCol_WindowBg, Defaults.color_background);

	ImGui::PushStyleColor(ImGuiCol_Header, Defaults.color_bar);
	ImGui::PushStyleColor(ImGuiCol_HeaderActive, Defaults.color_bar);
	ImGui::PushStyleColor(ImGuiCol_HeaderHovered, Defaults.color_bar);

	ImGui::PushStyleColor(ImGuiCol_TitleBg, Defaults.color_bar);
	ImGui::PushStyleColor(ImGuiCol_TitleBgActive, Defaults.color_bar);
	ImGui::PushStyleColor(ImGuiCol_TitleBgCollapsed, Defaults.color_bar);

	ImGui::PushStyleColor(ImGuiCol_CloseButton, Defaults.color_accent);
	ImGui::PushStyleColor(ImGuiCol_CloseButtonActive, Defaults.color_accent);
	ImGui::PushStyleColor(ImGuiCol_CloseButtonHovered, Defaults.color_accent_light);

	ImGui::PushStyleColor(ImGuiCol_ScrollbarGrab, Defaults.color_primary);
	ImGui::PushStyleColor(ImGuiCol_ScrollbarGrabHovered, Defaults.color_primary_light);
	ImGui::PushStyleColor(ImGuiCol_ScrollbarGrabActive, Defaults.color_primary);
	ImGui::PushStyleColor(ImGuiCol_ScrollbarBg, Defaults.color_frames);

	ImGui::PushStyleColor(ImGuiCol_CheckMark, Defaults.color_primary);

	ImGui::PushStyleColor(ImGuiCol_FrameBg, Defaults.color_frames);

	ImGui::PushStyleColor(ImGuiCol_SliderGrab, Defaults.color_primary);
	ImGui::PushStyleColor(ImGuiCol_SliderGrabActive, Defaults.color_primary_light);

	ImGui::PushStyleColor(ImGuiCol_Button, Defaults.color_primary);
	ImGui::PushStyleColor(ImGuiCol_ButtonActive, Defaults.color_primary);
	ImGui::PushStyleColor(ImGuiCol_ButtonHovered, Defaults.color_primary_light);

	ImGui::PushStyleColor(ImGuiCol_PopupBg, Defaults.color_background);

	ImGui::PushStyleColor(ImGuiCol_Border, Defaults.color_bar);

	ImGui::PushStyleColor(ImGuiCol_MenuBarBg, Defaults.color_menu);
	ImGui::PushStyleColor(ImGuiCol_ModalWindowDarkening, Defaults.color_modal_darkening);

	Callbacks::init(window);
}

Gui::~Gui()
{
	for (int i = 0; i < 23; i++)
		ImGui::PopStyleColor();
	for (int j = 0; j < 7; j++)
		ImGui::PopStyleVar();

	ImGui_ImplGlfwGL3_Shutdown();
}

void Gui::render()
{
    // ImGui new frame
    ImGui_ImplGlfwGL3_NewFrame();

    for(auto iterator = m_named_elements.begin(); iterator != m_named_elements.end();)
    {
		if (!(iterator->second))
			iterator = m_named_elements.erase(iterator);
		else
		{
			iterator->second->draw();
			if((iterator->second->hintedDeletion()))
				iterator = m_named_elements.erase(iterator);
			++iterator;
		}
    }

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    ImGui::Render();
}

bool Gui::wantsKeyboardInput()
{
    ImGuiIO &io = ImGui::GetIO();
	return (io.WantCaptureKeyboard);
}

bool Gui::wantsMouseInput()
{
    ImGuiIO &io = ImGui::GetIO();
	return (io.WantCaptureMouse);
}
}
