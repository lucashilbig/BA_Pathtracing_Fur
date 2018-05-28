#ifndef KIRK_TEXTURE_H
#define KIRK_TEXTURE_H

#include "externals/stb_image/stb_image.h"
#include "Color.h"
#include "KIRK/Utils/Log.h"
#include <GL/glew.h>
#include <algorithm>
#include <cstring>

#include <externals/lodepng/lodepng.h>

#define INVALID_OGL_VALUE 0xFFFFFFFF

namespace KIRK {

	/**
	This class is a wrapper for combining OpenGL and Raytracer Textures into one object.
	It can work and be fetched without OpenGL. The only methods which require previous initialization of OpenGL are:
	\code{.cpp}
	glUpload();
	glBind();
	\endcode
	*/
	class Texture
	{
	public:
		//Flags
		const unsigned char FLAG_TEXTURE_HAS_DATA		= 0x01;	/*!< Has this Texture any data assigned? */
		const unsigned char FLAG_GL_TEXTURE_VALID		= 0x02;	/*!< Has this Texture already been loaded to an OpenGL Buffer? */
		const unsigned char FLAG_GL_TEXTURE_UP_TO_DATE	= 0x04;	/*!< Is the OpenGL buffer up to date with the local copy? */
		const unsigned char FLAG_GL_TEXTURE_IS_BINDLESS	= 0x08;	/*!< Is the OpenGL buffer up to date with the local copy? */

        const unsigned char TEXTURE_WRAP_MODE_CLAMP = 0;
        const unsigned char TEXTURE_WRAP_MODE_TILE = 1;

		constexpr const static unsigned TEX_RGBA_UB = 4;
		constexpr const static unsigned TEX_RGB_UB = 3;
		constexpr const static unsigned TEX_RG_UB = 2;
		constexpr const static unsigned TEX_R_UB = 1;
		
		unsigned glSizedTypeForChannels(unsigned channels)
		{
			switch (channels)
			{
			case TEX_RGBA_UB: return GL_RGBA8;
			case TEX_RGB_UB: return GL_RGB8;
			case TEX_RG_UB: return GL_RG8;
			case TEX_R_UB: return GL_R8;
			default: return 4;
			}
		}

		unsigned glTypeForChannels(unsigned channels)
		{
			switch (channels)
			{
			case TEX_RGBA_UB: return GL_RGBA;
			case TEX_RGB_UB: return GL_RGB;
			case TEX_RG_UB: return GL_RG;
			case TEX_R_UB: return GL_RED;
			default: return 4;
			}
		}

		unsigned glChannelsForType(int type)
		{
			switch (type)
			{
			case GL_RGBA: return TEX_RGBA_UB;
			case GL_RGB: return TEX_RGB_UB;
			case GL_RG: return TEX_RG_UB;
			case GL_RED: return TEX_R_UB;
			default: return 4;
			}
		}

		/**
		Constructor for creating a Texture from loading a file.
		@param path Filename to load.
		*/
		Texture(const char* path, unsigned type = 0);

		/**
		Constructor for creating a Texture by giving general size, format and content parameters.
		@param width Texture width in pixels
		@param height Texture height in pixels
		@param alpha Whether the Texture should have an alpha channel or not
		@param data The data to initialize the Texture from. Can be left out, but keep in mind to assign it later before further accessing.
		*/
		Texture(const unsigned int width, const unsigned int height, unsigned channels, unsigned char *data = NULL);

		/**
		Constructor for creating a new Texture object from an existing OpenGL Texture. Caution: The textureID will be used for the newly created object and can interfere with other extisting Textures.
		@param textureID The texture ID of the OpenGL Texture
		*/
		Texture(const GLuint textureID);

		/**
		Copy constructor
		*/
		Texture(const Texture &texture);

		/**
		Default Destructor.
		*/
		~Texture();

		//METHODS

		/**
		If the texture is in RGB format, this method will convert it to RGBA, adding an alpha channel.
		*/
		void convertToRGBA();

		/**
		Use this method to load the texture into OpenGL buffers.
		*/
		void glUpload(bool bindless = false);

		/**
		Calls glBindTexture on this Texture ID, but only if the Texture is already loaded in an OpenGL buffer via glUpload().
		@param as_compute_storage True, if you want to use the OpenGL texture for writing to or reading from a compute shader.
		@param unit Texture unit to bind to.
		@param access Access mode for compute storage.
		*/
		void glBind(bool as_compute_storage = false, GLuint unit = 0, GLenum access = GL_WRITE_ONLY);

		/**
		Downloads the Texture stored on the GPU via OpenGL.
		*/
		void glDownload();

		/**
		Fetches the color of this texture on a given Texture coordinate.
		@param x Texture coordinate x value.
		@param y Texture coordinate y value.
		@return The color at the mapped Texture coordinate (x, y)
		*/
		KIRK::Color::RGBA getColor(const float x, const float y) const;

		/**
		Saves this Texture to a file.
		@param file The file name of the saved image.
		@return true, if the saving process has been successful.
		*/
		bool saveTo(const char* file) const ;

		/**
		Resizes the Texture. By calling this method, the Texture will also be invalidated, so it has to be updated and (if necesssary) uploaded to OpenGL for further use.
		*/
		void resize(const int width, const int height);

		/**
		Sets the given pixel to the given color.
		@param x Pixel x coordinate in texture image from the bottom right
		@param y Pixel y coordinate in texture image from the bottom right
		@param color The new rgba color for the pixel.
		*/
		void setPixel(const int x, const int y, const KIRK::Color::RGBA &color);

        /**
         Sets the Texture wrap mode to either TEXTURE_WRAP_MODE_CLAMP or TEXTURE_WRAP_MODE_TILE.
         @param mode The target wrap mode. Can be either TEXTURE_WRAP_MODE_CLAMP or TEXTURE_WRAP_MODE_TILE.
         */
        void setWrapMode(int mode);

		// OPERATORS
		/**
		Let's say that the Texture and it's data are contextually equivalent.
		*/
		operator unsigned char*() const {
			return m_texture_data;
		}

		// INLINES
		/**
		Updates the Texture data (the image). Does not (!) do any OpenGL interactions. For uploading the new data to the buffer, call <code>glUpload()</code> afterwards.
		@param newData The new image Data. Should have the same parameters (size, offsets etc.) as the old one.
		*/
		inline void updateData(unsigned char* newData) {
			//delete old data
			if (m_flags & FLAG_TEXTURE_HAS_DATA)
				delete[] m_texture_data;

			m_texture_data = newData;
			//Invalidate Texture
			m_flags |= FLAG_TEXTURE_HAS_DATA;
			m_flags &= ~FLAG_GL_TEXTURE_UP_TO_DATE;
		}

		/**
		@return The OpenGL Texture ID. Can be invalid, if there is none assigned yet.
		*/
		inline unsigned int id() const
        { return m_texture_id; }

		/**
		@return the raw texture data as unsigned char*
		*/
        unsigned char* getData() const
        { return m_texture_data; }

		/**
		@return true if this texture has any data assigned to it.
		*/
        bool textureHasData() const
		{ return m_flags & FLAG_TEXTURE_HAS_DATA; }

		/**
		@return The current texture size.
		*/
        glm::vec2 getSize() const
		{ return glm::vec2(m_width, m_height); }

		/**
		Creates and returns an empty texture data array sized with the current texture size.
		*/
        unsigned char* emptyData() const;

		GLuint64 getBindlessTextureID()
		{
			if (m_flags & FLAG_GL_TEXTURE_IS_BINDLESS)
			{
				return m_bindless_texture_id;
			}

			LOG_WARN("Accessing the bindless texture ID of a non-bindless texture.");
			return -1;
		}

        unsigned char getFlags() {return m_flags;}

	private:
		/**
		Initializes the Texture and sets the according data.
		(See the Texture Constructor for more details)
		*/
		void init(const int width, const int height, unsigned char *data);

		/**
		Converts the given data to an 4-component rgba format.
		@param The rgb (!) formatted data.
		*/
		void convertRGBA(unsigned char* &data);

		/**
		Converts a float to a byte clamped to 0..255
		@param f float value
		@return f*255 clamped to 0..255
		*/
		inline unsigned char toByte(const float f) const {
			return std::max(std::min(f * 255.f, 255.f), 0.0f);
		}

		unsigned char m_flags = 0; /*!< Giving the glm::vec4 color a more descriptive name and format. */

        unsigned char m_texture_wrap_mode = TEXTURE_WRAP_MODE_TILE; /*!< Determines, how the texture should be wrapped on exceeding it's uv-bounds. */
		std::vector<unsigned char> data2;

		//General data
		int m_width;					/*!< Texture width. */
		int m_height;					/*!< Texture height */
		int m_channels;
		unsigned char* m_texture_data;	/*!< Texture content */

		unsigned int m_texture_id;		/*!< OpenGL Texture ID */
		GLuint64 m_bindless_texture_id;

	};
};

#endif //KIRK_TEXTURE_H
