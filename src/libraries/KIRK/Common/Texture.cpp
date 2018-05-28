#include "Texture.h"


KIRK::Texture::Texture(const char* path, unsigned type)
{
	m_texture_id = INVALID_OGL_VALUE;

	int request_channels = type; 

	//load image from path
	unsigned char *data = stbi_load(path, &m_width, &m_height, &m_channels, request_channels);

	m_channels = request_channels == 0 ? m_channels : request_channels;

	data2 = std::vector<unsigned char>(data, data + m_width * m_height * m_channels);

	//flip image vertically
	unsigned char* s = data;
	for (int y = 0; y<m_height / 2; y++)
	{
		unsigned char* e = data + (m_height - y - 1)*m_width*m_channels;
		for (int x = 0; x<m_width*m_channels; x++)
		{
			unsigned char temp = *s;
			*s = *e;
			*e = temp;
			s++;
			e++;
		}
	}

	init(m_width, m_height, data);
}

KIRK::Texture::Texture(const unsigned int width, const unsigned int height, unsigned channels, unsigned char *data)
{
	m_texture_id = INVALID_OGL_VALUE;
	m_channels = channels;
	init(width, height, data);
}

KIRK::Texture::Texture(const GLuint textureID)
{
	m_texture_id = textureID;

	glDownload();
}

KIRK::Texture::Texture(const Texture &texture)
{
    m_texture_id = INVALID_OGL_VALUE;
    m_height = texture.m_height;
    m_width = texture.m_width;
	m_channels = texture.m_channels;
    m_texture_data = new unsigned char [m_width*m_height*m_channels];
    for (int i = 0; i < m_width*m_height*m_channels; i++)
        m_texture_data[i] = texture.m_texture_data[i];
    m_flags |= FLAG_TEXTURE_HAS_DATA;
}

KIRK::Texture::~Texture()
{
	if(m_flags & FLAG_TEXTURE_HAS_DATA)
		delete[] m_texture_data;
	if (m_flags & FLAG_GL_TEXTURE_VALID)
		glDeleteTextures(1, &m_texture_id);
}

void KIRK::Texture::convertToRGBA() {
	if (m_channels != 3)
		return;

	m_channels = 4;
	convertRGBA(m_texture_data);

	init(m_width, m_height, m_texture_data);

	m_flags &= ~FLAG_GL_TEXTURE_UP_TO_DATE;
	m_flags |= FLAG_TEXTURE_HAS_DATA;
}

void KIRK::Texture::convertRGBA(unsigned char* &data) {
	int data_size = m_width*m_height * 4;
	int inserted_alphas = 0;

	unsigned char* new_data = new unsigned char[data_size];

	for (int i = 0; i < data_size; i++) {
		if ((i + 1) % 4 == 0) {
			new_data[i] = 255;
			inserted_alphas++;
		}
		else {
			new_data[i] = data[i - inserted_alphas];
		}
	}
	delete[] data;
	data = new_data;
}

void KIRK::Texture::init(const int width, const int height, unsigned char* data)
{
	m_width = width;
	m_height = height;
    m_texture_data = data;

	//Do we have any data now? If not, say so.
    if (!data)
	{
        m_texture_data = emptyData();
    }

    m_flags |= FLAG_TEXTURE_HAS_DATA;
}

void KIRK::Texture::resize(const int width, const int height)
{
	//Invalidate Texture
	m_flags &= ~FLAG_GL_TEXTURE_UP_TO_DATE;

	//Set new measurement
	m_width = width;
	m_height = height;
    updateData(emptyData());
}

void KIRK::Texture::glUpload(bool bindless)
{
	//We cannot upload anything non-existent to OpenGL
	if (!(m_flags & FLAG_TEXTURE_HAS_DATA))
	{
		LOG_ERROR("Cannot load a Texture without data.");
		return;
	}

	//Do we have a valid textureID?
	if (!(m_flags & FLAG_GL_TEXTURE_VALID))
	{
		glGenTextures(1, &m_texture_id);
	}

	//Is this texture up to date? If not, we can finally upload the new data.
	if (!(m_flags & FLAG_GL_TEXTURE_UP_TO_DATE))
	{
		glBindTexture(GL_TEXTURE_2D, m_texture_id);

		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glPixelStorei(GL_PACK_ALIGNMENT, 1);

		auto type = glTypeForChannels(m_channels);

		if(m_channels == 2)
		{
			GLint swizzleMask[] = { GL_RED, GL_RED, GL_RED, GL_ALPHA };
			glTexParameteriv(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_RGBA, swizzleMask);
		}
		if (m_channels == 1)
		{
			GLint swizzleMask[] = { GL_RED, GL_RED, GL_RED, GL_ONE };
			glTexParameteriv(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_RGBA, swizzleMask);
		}

		glTexImage2D(GL_TEXTURE_2D, 0, type, m_width, m_height, 0, type, GL_UNSIGNED_BYTE, m_texture_data);
		glGenerateMipmap(GL_TEXTURE_2D);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		glBindTexture(GL_TEXTURE_2D, 0);
	}

	if (bindless)
	{
		m_flags |= FLAG_GL_TEXTURE_IS_BINDLESS;
		m_bindless_texture_id = glGetTextureHandleARB(m_texture_id);
		if(!glIsTextureHandleResidentARB(m_bindless_texture_id))
			glMakeTextureHandleResidentARB(m_bindless_texture_id);
	}
	else
	{
		m_flags &= ~FLAG_GL_TEXTURE_IS_BINDLESS;
	}

	m_flags |= FLAG_GL_TEXTURE_UP_TO_DATE;
	m_flags |= FLAG_GL_TEXTURE_VALID;
}

void KIRK::Texture::glBind(bool as_compute_storage, GLuint unit, GLenum access)
{
	if (m_flags & FLAG_GL_TEXTURE_VALID) {
		glActiveTexture(GL_TEXTURE0 + unit);
		if (as_compute_storage) {
			auto type = glSizedTypeForChannels(m_channels);
			glBindImageTexture(unit, m_texture_id, 0, GL_FALSE, 0, access, type);
		}
		else {
			glBindTexture(GL_TEXTURE_2D, m_texture_id);
		}
	}
}

void KIRK::Texture::glDownload() {
	glBindTexture(GL_TEXTURE_2D, m_texture_id);

	int format;
	glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &m_width);
	glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &m_height);
	glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_INTERNAL_FORMAT, &format);

	//Load Texture from OpenGL.
	m_channels = glChannelsForType(format);
    delete [] m_texture_data;
	m_texture_data = new unsigned char[m_width*m_height*m_channels];
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glGetTexImage(GL_TEXTURE_2D, 0, format, GL_UNSIGNED_BYTE, m_texture_data);
	
	//We should definitively have an OpenGL Texture ID and at least something as data.
	m_flags |= FLAG_TEXTURE_HAS_DATA;
	m_flags |= FLAG_GL_TEXTURE_VALID;
}

void  KIRK::Texture::setPixel(const int x, const int y, const Color::RGBA &color)
{
    //If texture has not been loaded properly, don't try to set a pixel in the first place.
    if(!(m_flags & FLAG_TEXTURE_HAS_DATA))
	{
        LOG_ERROR("No data available.");
        return;
    }

	int base_index = m_channels * (y * m_width + x);
	m_texture_data[base_index] = toByte(color.r);
	if (m_channels > 1)
		m_texture_data[base_index+1] = toByte(color.g);
	if (m_channels > 2)
		m_texture_data[base_index+2] = toByte(color.b);
	if(m_channels > 3)
		m_texture_data[base_index+3] = toByte(color.a);

	m_flags &= ~FLAG_GL_TEXTURE_UP_TO_DATE;
}

KIRK::Color::RGBA KIRK::Texture::getColor(const float x, const float y) const
{
    //If texture has not been loaded properly, show an error-like red color instead.
    if(!(m_flags & FLAG_TEXTURE_HAS_DATA))
        return Color::RED;

    if(std::isnan(x) || std::isnan(y))
	{
		//Error color.
        return Color::RED;
    }

    float uv_x = x;
    float uv_y = y;

    /*
    Clamp:
    x>=1.0 -> x = 1.0f;

    Wrap:
    x>=1.0 -> x = x-floor(x);
     */

    if(uv_x > 1.0 || uv_x < 0.0)
    {
        uv_x = glm::pow(x-glm::floor(x), m_texture_wrap_mode);
    }

    if(uv_y > 1.0 || uv_y < 0.0)
    {
        uv_y = glm::pow(y-glm::floor(y), m_texture_wrap_mode);
    }

	int sx = (int)(uv_x * (m_width-1));
	int sy = (int)(uv_y * (m_height-1));
	unsigned char *texel = m_texture_data + m_channels * (sy * m_width + sx);

	// c=2 -> rrra
	// c=1 -> rrr1
	// c=3 -> rgb1
	// c=4 -> rgba

	switch (m_channels)
	{
	case TEX_RGBA_UB: return { *texel / 255.f , *(texel + 1) / 255.f, *(texel + 2) / 255.f, *(texel + 3) / 255.f };
	case TEX_RGB_UB: return { *texel / 255.f , *(texel + 1) / 255.f, *(texel + 2) / 255.f, 1.f };
	case TEX_RG_UB: return { *texel / 255.f , *texel / 255.f, *texel / 255.f, *(texel + 1) / 255.f };
	case TEX_R_UB: return { *texel / 255.f , *texel / 255.f, *texel / 255.f, *texel / 255.f };
	default: return Color::MAGENTA;
	}
}

void KIRK::Texture::setWrapMode(int mode)
{
    m_texture_wrap_mode = mode;
}

unsigned char* KIRK::Texture::emptyData() const
{
    return new unsigned char[m_width*m_height*m_channels];
}

bool KIRK::Texture::saveTo(const char* file) const
{
    if (!(m_flags & FLAG_TEXTURE_HAS_DATA))
	{
        LOG_ERROR("Encoder error: No Texture data available in Texture.");
		return false;
	}

	LodePNGColorType format;
	
	switch (m_channels)
	{
	case TEX_RGBA_UB: format = LCT_RGBA; break;
	case TEX_RGB_UB: format = LCT_RGB; break;
	case TEX_RG_UB: format = LCT_GREY_ALPHA; break;
	case TEX_R_UB: format = LCT_GREY; break;
	}

	int format_multiplier = m_channels;
	int wh = m_width * m_height;

    //flip image vertically
    unsigned char* cpy_data = new unsigned char[wh*format_multiplier];
    std::memcpy(cpy_data, m_texture_data, wh*format_multiplier);
    unsigned char* s = cpy_data;
    for (int y = 0; y<m_height / 2; y++)
    {
        unsigned char* e = cpy_data + (m_height - y - 1)*m_width*format_multiplier;
        for (int x = 0; x<m_width*format_multiplier; x++)
        {
            unsigned char temp = *s;
            *s = *e;
            *e = temp;
            s++;
            e++;
        }
    }

	//write png
	std::vector<unsigned char> png;
	unsigned error = lodepng::encode(png, cpy_data, m_width, m_height, format);

    delete[] cpy_data;

    if (!error)
	{
        LOG_INFO("Saving image to %", file);
		return lodepng::save_file(png, file)==0;
	}

    LOG_ERROR("Encoder error %: %", error, lodepng_error_text(error));
	return false;
}
