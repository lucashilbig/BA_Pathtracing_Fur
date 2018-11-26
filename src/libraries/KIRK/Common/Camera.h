#ifndef KIRK_CAMERA_H
#define KIRK_CAMERA_H

#include <iostream>
#include <memory>

#include "KIRK/Common/Ray.h"

#include "SceneNode.h"

namespace KIRK
{
	/**
	This class can be seen as a virtual "real" camera with several parameters which all influence the rendered image.
	We let all ray calculations be done by the camera before tracing, instead of bloating up several ray-/pathtracer classes all with the same lines of code.

	@brief This class represents a virtual camera with all ray logic.
	*/
    class Camera : public NodeDataObject
    {
	public:

		/** Default constructor. Parameters should be set one by one for the sake of not having a constructor with like a thousand parameters*/
		Camera() {}

		/** Destructor. Doesn't do much here. */
        ~Camera();

		/**
		Applies all parameters set by the setters. Please call this after any setter. Otherwise it might lead to weird results.
		*/
		virtual void applyParameters();

		/**
		Calculates a ray in world coordinates which will go through a given pixel of the camera sensor.
		@param x x-pixel-coordinate on the sensor. Includes half pixel offsets.
		@param y y-pixel-cordinate on the sensor. Includes half pixel offsets.
		@param subpixel_delta_x x-pixel-offset from the bottom left of the pixel.
		@param subpixel_delta_y x-pixel-offset from the bottom left of the pixel.
		*/
		KIRK::Ray getRayFromPixel(const float x, const float y, const float subpixel_delta_x = 0.5f, const float subpixel_delta_y = 0.5f) const;

		/**
		Generates a ray from a general direction, calculating a depth of field jittering effect.<br/>
		For that, a random point on the virtual lense will be calculated and from there, the ray will look towards the general focus point of the old direction.
		@param direction General direction. Will be used for calculating the focus point here.
		@return A PixelRay that has been jittered to include depth of field effects.
		*/
		KIRK::Ray transformToDof(const glm::vec3 &direction) const;

		//Setters
		/** Call applyParameters() after setting any value. */
        void setResolution(const glm::vec2 &resolution);
		/** Call applyParameters() after setting any value. */
		void setTransform(const glm::vec3 &position, const glm::vec3 &look_at, const glm::vec3 &up);
		/** Call applyParameters() after setting any value. */
		void setFocus(const float &focus_distance);
		/** Call applyParameters() after setting any value. */
		void setSensorSize(const glm::vec2 &size);
		/** Call applyParameters() after setting any value. */
		void setFocalLength(const float focal_length);
		/** Call applyParameters() after setting any value. */
		void setFStop(const float f_stop);

		/**
		@return This camera's image resolution width.
		*/
        inline int getImageWidth() const { return m_resolution.x; }

		/**
		@return This camera's image resolution height.
		*/
        inline int getImageHeight() const { return m_resolution.y; }

		/**
		@return This camera's sensor width.
		*/
		inline float getSensorWidth() const { return m_sensor_size.x; }

		/**
		@return This camera's sensor height.
		*/
		inline float getSensorHeight() const { return m_sensor_size.y; }

		/**
		@return This camera's calculated horizontal field of view.
		*/
        inline float getFieldOfView() const { return m_field_of_view; }

		/**
		@return This camera's image aspect ratio.
		*/
        inline float getAspect() const { return m_aspect; }

		/** @return The camera pixel size.*/
		float getPixelSize() const;

        /**
        @return This camera's local Position.
        */
        inline glm::vec3 getLocalPosition() const { return m_local_position; }

        /**
        @return This camera's local look at vector.
        */
        inline glm::vec3 getLocalLookAt() const { return m_local_look_at; }

        /**
        @return This camera's local up vector.
        */
        inline glm::vec3 getLocalUp() const { return m_local_up; }

		/** Stream operator for printing out camera data. */
        friend std::ostream& operator<<(std::ostream& os, const Camera* camera);

        float m_focus_distance	= 11;						//!< Camera distance to sharpest focus plane (in meters) */

		glm::vec3 getBottomLeft() { return m_bottom_left; }

		glm::vec3 getPosition() { return m_position; }

		glm::vec3 getLookAt() { return m_look_at; }

		glm::vec3 getUp() { return m_up; }

		glm::vec3 getStepX();

		glm::vec3 getStepY();

	private:
		//The following preset values default to a vertical field of view of 60 degrees
		glm::vec2 m_sensor_size = glm::vec2(0.036, 0.024);	/*!< Sensor size (in meters) */
        float m_focal_length	= 0.0415f;					/*!< Distance from sensor to camera lense (in meters). Can be used for zooming. */
		float m_f_stop			= 1.8f;						/*!< F-Stop value. Indicator for the aperture and thus the depth of field size. (smaller = larger DOF) */
		float m_aspect;										/*!< Aspect ratio of the rendered image. Don't set it yourself! */
		float m_aperture;									/*!< Camera aperture. It's a relation between the aperture opening diameter and the focal length. Don't set it yourself! */
		float m_circle_of_confusion_limit;					/*!< Largest blur spot that will still be perceived as a point. Don't set it yourself! */
		//float m_hyperfocal_distance;						/*!< Minimal focus distance where far objects will be perceived sharp enough. Don't set it yourself! */
        float m_field_of_view;                              /*!< Field of view. Will be calculated on applyParameters(). Don't set it yourself!*/

		glm::vec3 m_local_position;							/*!< Camera position in local space. */
		glm::vec3 m_local_look_at;							/*!< Camera forward vector (aka looking direction) in local space. */
		glm::vec3 m_local_up;								/*!< Camera up vector in local space. */

		glm::vec3 m_position;								/*!< Camera position in world space. */
		glm::vec3 m_look_at;								/*!< Camera forward vector (aka looking direction) in world space. */
		glm::vec3 m_up;										/*!< Camera up vector in world space. */

		glm::vec3 m_axis_x;									/*!< Camera coordinate system x axis. Don't set it yourself! */
		glm::vec3 m_axis_y;									/*!< Camera coordinate system y axis. Don't set it yourself! */
		glm::vec3 m_axis_z;									/*!< Camera coordinate system z axis. Don't set it yourself! */

        glm::vec2 m_resolution;								/*!< The resolution in which the image will be rendered. */
		glm::vec3 m_bottom_left;							/*!< The bottom left most world point seen by the camera field of view. Don't set it yourself! */
		float m_pixel_size;									/*!< The size of a single pixel (in meters). Don't set it yourself! */
    };

    std::ostream& operator<<(std::ostream& os, const KIRK::Camera* camera);
}

#endif //KIRK_CAMERA_H
