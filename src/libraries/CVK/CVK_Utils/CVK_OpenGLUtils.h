
#ifndef __OPENGL_UTILS_H__
#define __OPENGL_UTILS_H__

#include <CVK/CVK_2/CVK_Framework.h>

namespace CVK {

namespace OpenGLUtils {
void checkError(const char *info = "");

class Window
{
public:
    Window(std::string title, int width, int height, uint32_t gl_version_major = 3, uint32_t gl_version_minor = 3);
    virtual ~Window();

    void endFrame();
    void setBackgroundColor(glm::vec3 color);

    inline bool shouldClose() { return glfwWindowShouldClose(m_window); }

	inline void close() { glfwSetWindowShouldClose(m_window, GL_TRUE); }

    inline GLFWwindow *getGLFWWindow() const { return m_window; }

    inline int getWidth() { return m_width; }

    inline int getHeight() { return m_height; }
    
    inline int getVersionCode() { return m_ogl_version_code; }

private:
    GLFWwindow *m_window;
    int m_width, m_height;
    
    int m_ogl_version_code;
};

};


};

#endif //__OPENGL_UTILS_H__
