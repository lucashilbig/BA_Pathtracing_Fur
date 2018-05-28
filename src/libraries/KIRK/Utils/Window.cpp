#include "Window.h"

#include <CVK/CVK_2/CVK_Framework.h>
#include <KIRK/Utils/Log.h>

namespace KIRK
{

	Window::Window(std::string title, int width, int height, unsigned gl_version_major, unsigned gl_version_minor)
	{
		glfwInit();

        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, gl_version_major);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, gl_version_minor);
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
		glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

        glewExperimental = GL_TRUE;

		m_width = width;
		m_height = height;
		m_window = glfwCreateWindow(width, height, title.c_str(), nullptr, nullptr);
		glfwMakeContextCurrent(m_window);

		//Now check and print a warning if
		int v_major, v_minor;
		glGetIntegerv(GL_MAJOR_VERSION, &v_major);
		glGetIntegerv(GL_MINOR_VERSION, &v_minor);
		m_ogl_version_code = 100 * v_major + v_minor;
		uint32_t hinted_gl_version = 100 * gl_version_major + gl_version_minor;

		if (m_ogl_version_code < hinted_gl_version) {
			LOG_WARN("Hinted a higher OpenGL version than supported by your system. (hinted %.%, available %.%)", gl_version_major, gl_version_minor, v_major, v_minor);
		}

		glfwSetWindowPos(m_window, 300, 50);

		//init opengl 3 extensions
		glewInit();
		setBackgroundColor(Black);
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		CVK::State::getInstance()->updateSceneSettings(darkgrey, NO_FOG, white, 0, 0, 1);
		CVK::State::getInstance()->setWindow(m_window);
	}

	Window::~Window()
	{
		glfwDestroyWindow(m_window);
		glfwTerminate();
	}

	void Window::endFrame() const
	{
		// show what's been drawn
		glfwSwapBuffers(m_window);
		//if (omp_get_thread_num() == 0)
		glfwPollEvents();
	}

	void Window::setBackgroundColor(glm::vec3 color)
	{
		CVK::State::getInstance()->setBackgroundColor(color);
		glClearColor(color.r, color.g, color.b, 0.0);
	}

	bool Window::extensionsSupported(std::vector<std::string> extensions)
	{
		bool all_supported = true;
		for(auto &&extension : extensions)
		{
			if(!glfwExtensionSupported(extension.c_str()))
			{
				LOG_ERRORs() << "Extension not supported: " << extension;
				all_supported = false;
			} 
			else
			{
				LOG_INFOs() << "Extension supported: " << extension;
			}
		}
		return all_supported;
	}
}
