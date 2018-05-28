#ifndef __KIRK_WINDOW_H
#define __KIRK_WINDOW_H

#include <string>
#include <GL/glew.h>
#define GLFW_INCLUDE_GLEXT
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <vector>

namespace KIRK 
{
class Window
{
public:
	Window(std::string title, int width, int height, unsigned gl_version_major = 3, unsigned gl_version_minor = 3);
	virtual ~Window();

	void endFrame() const;
	static void setBackgroundColor(glm::vec3 color);

	bool shouldClose() const { return glfwWindowShouldClose(m_window); }

	void close() const { glfwSetWindowShouldClose(m_window, GL_TRUE); }

	GLFWwindow *getGLFWWindow() const { return m_window; }

	int getWidth() const { return m_width; }

	int getHeight() const { return m_height; }

	int getVersionCode() const { return m_ogl_version_code; }

	static bool extensionsSupported(std::vector<std::string> extensions);

private:
	GLFWwindow *m_window;
	int m_width, m_height;

	int m_ogl_version_code;
};
}

#endif //!__KIRK_WINDOW_H