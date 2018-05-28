#include "CVK_FreeFlight.h"
#include "KIRK/Utils/KeyMapping.h"

namespace CVK {
	FreeFlight::FreeFlight(int width, int height, CVK::Projection *projection)
	{
		m_cameraPos = glm::vec3(0.0f, 0.0f, 0.0);
		m_positiondefault = m_cameraPos;
		m_direction = glm::vec3(0.0f, 0.0f, -1.0f);
		m_dirdefault = m_direction;
		m_up = glm::vec3(0.0f, 1.0f, 0.0f);

		m_width = width;
		m_height = height;

		m_sensitivity = 0.01f;
		m_theta = glm::pi<float>() / 2.0f;
		m_phi = 0.f;

		m_altitudespeed = 0.002;
		m_flyspeed = 0.002;

		m_viewmatrix = glm::lookAt(m_cameraPos, m_direction, m_up);

		m_oldX = width / 2.f;
		m_oldY = height / 2.f;

		if (projection == nullptr) projection = new CVK::Perspective(width / static_cast <float>(height));
		setProjection(projection);
	}

	FreeFlight::~FreeFlight()
	{
	}

	void FreeFlight::update(GLFWwindow *window)
	{
		auto currentTime = std::chrono::steady_clock::now();
		auto deltaTime = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - m_clock).count();
		m_clock = currentTime;

		/*
		* HOME: Reset
		*/
		if (glfwGetKey(window, KIRK_KEY_CAMERA_RESET) ||
			glfwGetKey(window, KIRK_KEY_CAMERA_RESET_ALT) == GLFW_PRESS) {
			m_direction = m_dirdefault;
			m_cameraPos = m_positiondefault;
			m_theta = glm::pi<float>() / 2.0f;
			m_phi = 0.f;
		}

		/*
		* L-CTRL & L-SHIFT speed modifier toggles
		*/
		if (glfwGetKey(window, KIRK_KEY_CAMERA_SPEEDMODIFIER_SLOWER) ||
			glfwGetKey(window, KIRK_KEY_CAMERA_SPEEDMODIFIER_SLOWER_ALT) == GLFW_PRESS) {
			deltaTime *= 0.25f;
		}
		if (glfwGetKey(window, KIRK_KEY_CAMERA_SPEEDMODIFIER_FASTER) ||
			glfwGetKey(window, KIRK_KEY_CAMERA_SPEEDMODIFIER_FASTER_ALT) == GLFW_PRESS) {
			deltaTime *= 4.0f;
		}

		double x, y;
		glfwGetCursorPos(window, &x, &y);

		/*
		* Phi and Theta rotation angles
		*/
		if (glfwGetMouseButton(window, KIRK_MOUSE_CAMERA_ROTATION_CENTER) ||
			glfwGetMouseButton(window, KIRK_MOUSE_CAMERA_ROTATION_CAMERA) == GLFW_PRESS)
		{
			float changeX = (static_cast <float>(x) - m_oldX) * m_sensitivity;
			float changeY = (static_cast <float>(y) - m_oldY) * m_sensitivity;
			glm::vec3 xAxis = glm::normalize(glm::cross(m_direction, m_up));
			m_direction = glm::vec3(glm::rotate(glm::mat4(1.0f), -changeX, m_up) * glm::vec4(m_direction, 0.0f));
			glm::vec3 tempcenter = glm::vec3(glm::rotate(glm::mat4(1.0f), -changeY, xAxis) * glm::vec4(m_direction, 0.0f));
			if(!(glm::abs(glm::dot(tempcenter, m_up)) > epsilon))
				m_direction = tempcenter;
		}

		float altitudeAddend = m_altitudespeed * deltaTime;
		if (glfwGetKey(window, KIRK_KEY_CAMERA_STRAFE_UP) || glfwGetKey(window, KIRK_KEY_CAMERA_STRAFE_UP_ALT) == GLFW_PRESS)
			m_cameraPos += glm::normalize(glm::cross(glm::cross(m_direction, m_up), m_direction)) * altitudeAddend;
		if (glfwGetKey(window, KIRK_KEY_CAMERA_STRAFE_DOWN) || glfwGetKey(window, KIRK_KEY_CAMERA_STRAFE_DOWN_ALT) == GLFW_PRESS)
			m_cameraPos -= glm::normalize(glm::cross(glm::cross(m_direction, m_up), m_direction)) * altitudeAddend;

		float flyAddend = m_flyspeed * deltaTime;
		if (glfwGetKey(window, KIRK_KEY_CAMERA_STRAFE_LEFT) || glfwGetKey(window, KIRK_KEY_CAMERA_STRAFE_LEFT_ALT) == GLFW_PRESS)
			m_cameraPos -= glm::normalize(glm::cross(m_direction, m_up)) * flyAddend;
		if (glfwGetKey(window, KIRK_KEY_CAMERA_STRAFE_RIGHT) || glfwGetKey(window, KIRK_KEY_CAMERA_STRAFE_RIGHT_ALT) == GLFW_PRESS)
			m_cameraPos += glm::normalize(glm::cross(m_direction, m_up)) * flyAddend;
		if (glfwGetKey(window, KIRK_KEY_CAMERA_MOVE_FORWARD) || glfwGetKey(window, KIRK_KEY_CAMERA_MOVE_FORWARD_ALT) == GLFW_PRESS)
			m_cameraPos += m_direction * flyAddend;
		if (glfwGetKey(window, KIRK_KEY_CAMERA_MOVE_BACKWARD) || glfwGetKey(window, KIRK_KEY_CAMERA_MOVE_BACKWARD_ALT) == GLFW_PRESS) {
			/*
			* Don't move backwards if Control is pressed to prevent SaveTexture camera Movement.
			*/
			if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_RELEASE && glfwGetKey(window, GLFW_KEY_RIGHT_CONTROL) == GLFW_RELEASE)
				m_cameraPos -= m_direction * flyAddend;
		}

		m_oldX = static_cast <float>(x);
		m_oldY = static_cast <float>(y);

		m_viewmatrix = glm::lookAt(m_cameraPos, m_cameraPos + m_direction, m_up);
	}

	void FreeFlight::setUpvector(glm::vec3 up)
	{
		m_up = glm::normalize(up);
	}

	void FreeFlight::setPosition(glm::vec3 pos)
	{
		m_cameraPos = pos;
		m_positiondefault = pos;
	}

}
