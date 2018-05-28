#include "Environment.h"
#include "Shading/EnvironmentShader.h"

KIRK::Environment::~Environment() {
	cleanUp();
}

KIRK::Environment::Environment()
{
	m_eshader = KIRK::ShaderFactory::getInstance().getShader("EnvironmentShader");
}

KIRK::Environment::Environment(Environment& obj)
{	
	m_eshader = KIRK::ShaderFactory::getInstance().getShader("EnvironmentShader");
    if (obj.m_has_cubemap)
    {
        for (int i = 0; i < 6; i++)
            m_cubemap_textures[i] = new KIRK::Texture(*obj.m_cubemap_textures[i]);
        m_has_cubemap = true;
        m_type = KIRK::Environment::CUBE_MAP;
    }
    else if (obj.m_has_spheremap)
    {
        m_spheremap_texture = new KIRK::Texture(*obj.m_spheremap_texture);
        m_has_spheremap = true;
        m_type = KIRK::Environment::SPHERE_MAP;
    }
    else
    {
        m_color = obj.getBackgroundColor();
        m_type = KIRK::Environment::COLOR;
    }
    setAmbientLight(obj.getAmbientLight());
}

void KIRK::Environment::cleanUp() {
	if (m_has_cubemap) {
		for (int i = 0; i < 6; i++) {
			delete m_cubemap_textures[i];
		}
		m_has_cubemap = false;
	}
	if (m_has_spheremap) {
		delete m_spheremap_texture;
		m_has_spheremap = false;
	}
}

void KIRK::Environment::loadCubeMap(const char* posx, const char* posy, const char* posz, const char* negx, const char* negy, const char* negz)
{
	cleanUp();
    
    LOG_DEBUG("Setting Cubemap as Environment...");
    
    //Let's sort the cubemap textures strategically to shorten the getColor method.
    m_cubemap_textures[0] = new Texture(posx); //x
    m_cubemap_textures[1] = new Texture(posy); //y
    m_cubemap_textures[2] = new Texture(posz); //z
    
    m_cubemap_textures[3] = new Texture(negx); //-x
    m_cubemap_textures[4] = new Texture(negy); //-y
    m_cubemap_textures[5] = new Texture(negz); //-z
    
    m_has_cubemap = true;
    m_type = CUBE_MAP;
}

void KIRK::Environment::setColor(KIRK::Color::RGBA color)
{
	cleanUp();
    
    LOG_DEBUG("Setting Color as Environment...");

    m_color = color;
    m_type = COLOR;
}

void KIRK::Environment::loadSphereMap(const char* map)
{
	cleanUp();
    
    LOG_DEBUG("Setting Spheremap as Environment...");
    
    m_spheremap_texture = new Texture(map);
    
    m_type = SPHERE_MAP;
    m_has_spheremap = true;
}

KIRK::Color::RGBA KIRK::Environment::getColor(const glm::vec3 &ray_direction) const
{
    glm::vec3 direction = glm::normalize(ray_direction);
    
    switch (m_type) {
        case CUBE_MAP:
        {
            glm::vec3 signs(glm::sign(direction.x), glm::sign(direction.y), glm::sign(direction.z));
            glm::vec3 absolutes(glm::abs(direction.x), glm::abs(direction.y), glm::abs(direction.z));
            float max = glm::max(absolutes.x, absolutes.y, absolutes.z);
            
            int side;
            glm::vec2 uv;
            
            if(max == absolutes.x){
                //Take either right or left map
                side = 0 + 1.5f - 1.5f * signs.x;
                uv = {(direction.z/direction.x+1)/2, (direction.y/absolutes.x+1)/2};
            } else if(max == absolutes.y) {
                //Take either top or bottom map
                side = 1 + 1.5f - 1.5f * signs.y;
                uv = {(direction.x/absolutes.y+1)/2, (direction.z/direction.y+1)/2};
            } else {
                //Take either front or back map
                side = 2 + 1.5f + 1.5f * signs.z;
                uv = {-(direction.x/direction.z+1)/2, (direction.y/absolutes.z+1)/2};
            }
            
            return m_cubemap_textures[side]->getColor(uv.x, uv.y);
        }
        case SPHERE_MAP:
        {
            float m = 2.0 * sqrt( direction.x*direction.x + direction.y*direction.y + (direction.z+1.0)*(direction.z+1.0) );
            glm::vec2 uv;
            uv.x = direction.x/m + 0.5;
            uv.y = direction.y/m + 0.5;
            return m_spheremap_texture->getColor(uv.x, uv.y);
        }
        case COLOR:
        default:
            return m_color;
    }
}
