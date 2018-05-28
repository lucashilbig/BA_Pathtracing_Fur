#include "Camera.h"

KIRK::Camera::~Camera()
{}

void KIRK::Camera::applyParameters()
{
	//Go up through node tree to calculate actual world transform
	glm::mat4 transform = calculateTransform();

	m_position = glm::vec3(transform * glm::vec4(m_local_position, 1));
	m_look_at = glm::vec3(transform * glm::vec4(m_local_look_at, 0));
	m_up = glm::vec3(transform * glm::vec4(m_local_up, 0));

    m_aperture = m_focal_length / m_f_stop;

    //Does nothing yet, but could be useful in the future when creating nicer Bokeh effects or something
    m_circle_of_confusion_limit = (m_focal_length * m_aperture) / (m_focus_distance - m_focal_length);

    // Creating a coordinate system for camera space
    m_axis_z = glm::normalize(-m_look_at);
    m_axis_x = glm::normalize(glm::cross(m_up, m_axis_z));
    m_axis_y = glm::normalize(glm::cross(m_axis_z, m_axis_x));

    //The field of view is no dedicated value but dependant from sensor size and focal length
    float sensor_diameter = glm::sqrt((m_sensor_size.x * m_sensor_size.x) + (m_sensor_size.y * m_sensor_size.y));
    m_field_of_view = 2 * glm::atan(sensor_diameter/ (2 * m_focal_length));

    // height sy to m_height is proproportional to sx to m_width
    float sy = glm::tan(0.5f * m_field_of_view);
    float sx = sy * m_aspect;

    m_pixel_size = (2.f * sx / (float)m_resolution.x);

    // calculating the position of the bottom left corner of the image
    m_bottom_left = m_position - m_axis_z - sy * m_axis_y - sx * m_axis_x;
}

KIRK::Ray KIRK::Camera::transformToDof(const glm::vec3 &direction) const
{
    glm::vec3 focus_point = m_position + m_focus_distance * direction;
    // Get random vector on disk
    glm::vec2 disk_rnd = glm::diskRand(m_aperture*3);
    // creating new camera position(or ray start using jittering)
    glm::vec3 start = m_position + disk_rnd.x*m_axis_x + disk_rnd.y*m_axis_y;
    // getting the new direction of ray
    glm::vec3 ndir = focus_point - start;
    glm::normalize(ndir);

    return KIRK::Ray (start, ndir);
}

////////////////////////////////////////////////////////////
///////////////
///////////////		GETTERS
///////////////
////////////////////////////////////////////////////////////

KIRK::Ray KIRK::Camera::getRayFromPixel(const float x, const float y, const float subpixel_delta_x, const float subpixel_delta_y) const
{
	glm::vec3 direction = m_bottom_left
		+ (x + subpixel_delta_x) * m_pixel_size * m_axis_x
		+ (y + subpixel_delta_y) * m_pixel_size * m_axis_y
		- m_position;
	return KIRK::Ray(m_position, direction);
}

float KIRK::Camera::getPixelSize() const
{
	return m_pixel_size;
}

glm::vec3 KIRK::Camera::getStepX() {
	return m_pixel_size * m_axis_x;
}

glm::vec3 KIRK::Camera::getStepY() {
	return m_pixel_size * m_axis_y;
}


////////////////////////////////////////////////////////////
///////////////
///////////////		SETTERS
///////////////
////////////////////////////////////////////////////////////

void KIRK::Camera::setResolution(const glm::vec2 &resolution)
{
    m_resolution = resolution;
    m_aspect = m_resolution.x / m_resolution.y;
}

void KIRK::Camera::setTransform(const glm::vec3 &position, const glm::vec3 &look_at, const glm::vec3 &up)
{
    m_local_position = position;
    m_local_look_at = look_at;
    m_local_up = up;
}

void KIRK::Camera::setSensorSize(const glm::vec2 &size)
{
    m_sensor_size = size;
}

void KIRK::Camera::setFocalLength(const float focal_length)
{
    m_focal_length = focal_length;
}

void KIRK::Camera::setFStop(const float f_stop)
{
    m_f_stop = f_stop;
}

void KIRK::Camera::setFocus(const float &focus_distance) 
{
	m_focus_distance = focus_distance;
}


////////////////////////////////////////////////////////////
///////////////
///////////////		OPERATORS
///////////////
////////////////////////////////////////////////////////////

namespace KIRK
{
	std::ostream& operator<<(std::ostream& os, const KIRK::Camera* camera)
	{
		if (camera)
		{
			//Not sure if we need those
			os << "Printing Camera: " << std::endl
			   << "Focal Length: " << camera->m_focal_length << std::endl
			   << "F Stop: " << camera->m_f_stop << std::endl
			   << "-----" << std::endl

			   << "World Position: ("
			   << camera->m_position.x << ", "
			   << camera->m_position.y << ", "
			   << camera->m_position.z <<")" << std::endl
			   << "-----" << std::endl

			   << "Local Position: ("
			   << camera->m_local_position.x << ", "
			   << camera->m_local_position.y << ", "
			   << camera->m_local_position.z <<")" << std::endl
			   << "-----" << std::endl

			   << "World Look at: ("
			   << camera->m_look_at.x << ", "
			   << camera->m_look_at.y << ", "
			   << camera->m_look_at.z <<")" << std::endl
				  << "-----" << std::endl

			   << "Local Look at: ("
			   << camera->m_local_look_at.x << ", "
			   << camera->m_local_look_at.y << ", "
			   << camera->m_local_look_at.z <<")" << std::endl
			   << "-----" << std::endl

			   << "World up vector: ("
			   << camera->m_up.x << ", "
			   << camera->m_up.y << ", "
			   << camera->m_up.z <<")" << std::endl
			   << "-----" << std::endl

			   << "Local up vector: ("
			   << camera->m_local_up.x << ", "
			   << camera->m_local_up.y << ", "
			   << camera->m_local_up.z <<")" << std::endl
			   << "-----" << std::endl;

			return os <<std::endl;
		}
		else
			return os << "Well ... seems like there is no camera :/" << std::endl;
	}
}
