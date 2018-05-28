#include "KIRK/Utils/Log.h"
#include "CVK_OpenGLUtils.h"

void CVK::OpenGLUtils::checkError(const char *info)
{
    for(GLenum err; (err = glGetError()) != GL_NO_ERROR;)
    {
        std::string error;
        switch(err)
        {
            case GL_INVALID_ENUM:
                error = "GL_INVALID_ENUM";
                break;
            case GL_INVALID_VALUE:
                error = "GL_INVALID_VALUE";
                break;
            case GL_INVALID_OPERATION:
                error = "GL_INVALID_OPERATION";
                break;
            case GL_STACK_OVERFLOW:
                error = "GL_STACK_OVERFLOW";
                break;
            case GL_STACK_UNDERFLOW:
                error = "GL_STACK_UNDERFLOW";
                break;
            case GL_OUT_OF_MEMORY:
                error = "GL_OUT_OF_MEMORY";
                break;
            case GL_INVALID_FRAMEBUFFER_OPERATION:
                error = "GL_INVALID_FRAMEBUFFER_OPERATION";
                break;
            default:
                error = "Unknown OpenGL Error";
                break;
        }

        std::string additional = "";
        if(std::string(info) != "")
            additional = " in " + std::string(info);

        std::cerr << "ERROR: " << error << " (" << err << ")" << additional << std::endl;
    }
}

CVK::OpenGLUtils::Window::Window(std::string title, int width, int height, uint32_t gl_version_major, uint32_t gl_version_minor)
{
    glfwInit();

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, gl_version_major);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, gl_version_minor);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	glewExperimental = GL_TRUE;

    m_width = width;
    m_height = height;
    m_window = glfwCreateWindow(width, height, title.c_str(), 0, 0);
    glfwMakeContextCurrent(m_window);

    //Now check and print a warning if
    int v_major, v_minor;
    glGetIntegerv(GL_MAJOR_VERSION, &v_major);
    glGetIntegerv(GL_MINOR_VERSION, &v_minor);
    m_ogl_version_code = 100 * v_major + v_minor;
	uint32_t hinted_gl_version = 100 * gl_version_major + gl_version_minor;

    if(m_ogl_version_code < hinted_gl_version){
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

CVK::OpenGLUtils::Window::~Window()
{
    // cleanup
    LOG_INFO("deletion ------------------------------------------");
    glfwDestroyWindow(m_window);
    glfwTerminate();
}

void CVK::OpenGLUtils::Window::endFrame()
{
    // show what's been drawn
    glfwSwapBuffers(m_window);
    //if (omp_get_thread_num() == 0)
    glfwPollEvents();
}

void CVK::OpenGLUtils::Window::setBackgroundColor(glm::vec3 color)
{
    CVK::State::getInstance()->setBackgroundColor(color);
    glClearColor(color.r, color.g, color.b, 0.0);
}
