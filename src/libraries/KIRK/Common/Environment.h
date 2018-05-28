#ifndef KIRK_ENVIRONMENT_H
#define KIRK_ENVIRONMENT_H

#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include "Color.h"
#include "Texture.h"

#include "Shading/ShaderFactory.h"
#include "Shading/Shader.h"

namespace KIRK
{
	/**
	An environment is used for having a CubeMap, SphereMap or just a Color as background for rendering.
	*/
	class Environment
	{
	public:
		std::shared_ptr<KIRK::Shader> m_eshader;
		/**
		Defines a type of which an Environment can be.
		*/
		enum Type
		{
            COLOR,
            CUBE_MAP,
            SPHERE_MAP
        };
        /**
         * @brief Default constructor
         */
        Environment();

        /** Destructor which runs a cleanUp function found in Environment.cpp */
        ~Environment();


		/**
         Copy Constructor to copy the textures
         */
        Environment(Environment& obj);

		/**
		Loads a CubeMap by stitching together the 6 side textures and replaces the current environment.
		@param posx Positive X side texture path
		@param posy Positive Y side texture path
		@param posz Positive Z side texture path
		@param negx Negative X side texture path
		@param negy Negative Y side texture path
		@param negz Negative Z side texture path
		*/
		void loadCubeMap(const char* posx, const char* posy, const char* posz, const char* negx, const char* negy, const char* negz);

		/**
		Loads a SphereMap and replaces the current environment.
		@param map SphereMap texture path
		*/
		void loadSphereMap(const char* map);

		/**
		Remove all maps (if there are any) and just set a color as background.
		@color The background colorc
		*/
        void setColor(const KIRK::Color::RGBA color);

        /**
         * @brief sets the ambient light color
         * @param ambient the color of the ambient light
         */
        void setAmbientLight(Color::RGBA ambient){m_ambient = ambient;};

        /**
        * @brief returns the ambient light color
        */
        Color::RGBA getAmbientLight() {return m_ambient;};

		/**
		Retrieves the color for a given ray direction. All calculations for CubeMaps and SphereMaps are contained within this method.
		@param ray_direction Direction of the ray.
		*/
        Color::RGBA getColor(const glm::vec3 &ray_direction) const;

		/** Deletes all cube map textures and the sphere map texture. At least if they exist. */
		void cleanUp();

        Type m_type = COLOR;                    //!< Current type of environment. COLOR, CUBE_MAP or SPHERE_MAP

        /**
        * @brief returns the environment color
        */
        Color::RGBA getBackgroundColor() {return m_color;};

        Texture* m_cubemap_textures[6];			//!< All CubeMap textures, sorted like the parameters in loadCubeMap.
        Texture *m_spheremap_texture = nullptr;	//!< The current SphereMap texture or nullptr if there is none

        struct GPU_Environment
        {
            glm::vec4 backgroundColor;
            glm::vec4 ambient;
            GLuint64 cubeMapIDs[6];
            GLuint64 sphereMap;

            GLuint type;
            GLuint hasCubeMap;
            GLuint hasSphereMap;
            GLuint padding;
        };

        GPU_Environment getGPUEnvironment()
        {
            GLuint64 cubeMapIDs[6];
            for (int cmt = 0; cmt < 6; cmt++)
            {
                if (m_has_cubemap)
                {
                    m_cubemap_textures[cmt]->glUpload(true);
                    cubeMapIDs[cmt] = m_cubemap_textures[cmt]->getBindlessTextureID();
                }
                else
                {
                    cubeMapIDs[cmt] = -1;
                }
            }

            GLuint64 sphereMapID;
            if (m_has_spheremap)
            {
                m_spheremap_texture->glUpload(true);
                sphereMapID = m_spheremap_texture->getBindlessTextureID();
            }
            else
            {
                sphereMapID = -1;
            }

            return GPU_Environment{
                    m_color,
                    m_ambient,
                {cubeMapIDs[0],
                    cubeMapIDs[1],
                    cubeMapIDs[2],
                    cubeMapIDs[3],
                    cubeMapIDs[4],
                    cubeMapIDs[5]},
                    sphereMapID,

                    (GLuint)m_type,
                    (GLuint)m_has_cubemap,
                    (GLuint)m_has_spheremap,
            };
        }

    private:
        Color::RGBA m_color;				//!< Current background color

        Color::RGBA m_ambient{ 0.1f, 0.1f, 0.1f, 1.0f }; //!< ambient light color

        bool m_has_cubemap = false;				//!< true, if the CubeMap textures are set.
        bool m_has_spheremap = false;			//!< true, it the SphereMap texture is set.

    };
}

#endif //KIRK_ENVIRONMENT_H
