#include "CVK_Material.h"

int CVK::Material::count = 0;

CVK::Material::Material(glm::vec3 diffuse, glm::vec3 specular, float shininess, float ior)
{
    init(1.f, diffuse, 1.f, specular, shininess, ior);
}

CVK::Material::Material(float kd, glm::vec3 diffuse, float ior)
{
    init(kd, diffuse, 0.f, black, 0.f, ior);
}

CVK::Material::Material(float kd, glm::vec3 diffuse, float ks, glm::vec3 specular, float shininess, float ior,
                        float alpha, int illum)
{
    init(kd, diffuse, ks, specular, shininess, ior, alpha, illum);
}

CVK::Material::Material(glm::vec3 kd, glm::vec3 diffuse, glm::vec3 ks, glm::vec3 specular, float shininess, float ior,
                        float alpha, int illum)
{
    init(kd, diffuse, ks, specular, shininess, ior, alpha, illum);
}

CVK::Material::Material(const char *colorTexturePath, float kd, float ior)
{
    init(kd, white, 0.f, black, 0.f, ior);
    setTexture(COLOR_TEXTURE, colorTexturePath);
}

CVK::Material::Material(GLuint colorTextureID, float kd, float ior)
{
    init(kd, white, 0.f, black, 0.f, ior);
    setTexture(COLOR_TEXTURE, colorTextureID);
}

CVK::Material::Material(const char *colorTexturePath, float kd, float ks, glm::vec3 specular, float shininess,
                        float ior, float alpha, int illum)
{
    init(kd, white, ks, specular, shininess, ior, alpha, illum);
    setTexture(COLOR_TEXTURE, colorTexturePath);
}


CVK::Material::Material(GLuint colorTextureID, float kd, float ks, glm::vec3 specular, float shininess, float ior)
{
    init(kd, white, ks, specular, shininess, ior);
    setTexture(COLOR_TEXTURE, colorTextureID);
}

void
CVK::Material::init(float kd, glm::vec3 diffuse, float ks, glm::vec3 specular, float shininess, float ior, float alpha,
                    int illum)
{
    // we have a vec and a float for the kd, ks values because it can differ for each canal
    m_kd = kd;
    m_ks = ks;
    m_kdVec = glm::vec3(kd);
    m_ksVec = glm::vec3(ks);

    // value for transparency - default it is set to 1.0
    m_alpha = alpha;

    // illumination model - default it is set to 0
    m_illum = illum;

    // the two values for reflection and refraction - default is ior set to 1.0
    m_shininess = shininess;
    m_ior = ior;

    // dont know what the value is for so far
    m_kt = 0.f;

    // those are implemented before from somebody else
    m_diffColor = diffuse;
    m_specColor = specular;

    m_colorTexturePath = "";
    m_normalTexturePath = "";

    m_colorTexture = NULL;
    m_normalTexture = NULL;

    m_id = count;
    count++;
}

void CVK::Material::init(glm::vec3 kd, glm::vec3 diffuse, glm::vec3 ks, glm::vec3 specular, float shininess, float ior,
                         float alpha, int illum)
{
    // we have a vec and a float for the kd, ks values because it can differ for each canal
    // if we get a vec with different values in each canal, we take the max value
    m_kd = glm::max(kd.r, kd.g);
    m_kd = glm::max(m_kd, kd.b);
    m_ks = glm::max(ks.r, ks.g);
    m_ks = glm::max(m_ks, ks.b);
    m_kdVec = kd;
    m_ksVec = ks;

    // the two values for reflection and refraction
    m_shininess = shininess;
    m_ior = ior;

    // value for transparency
    m_alpha = alpha;

    // illumination model
    m_illum = illum;

    // dont know what the value is for so far
    m_kt = 0.f;

    // those are implemented before from somebody else
    m_diffColor = diffuse;
    m_specColor = specular;

    m_colorTexturePath = "";
    m_normalTexturePath = "";

    m_colorTexture = NULL;
    m_normalTexture = NULL;

    m_id = count;
    count++;
}

CVK::Material::~Material()
{
}

void CVK::Material::setKd(float kd)
{
    m_kd = kd;
}

float CVK::Material::getKd()
{
    return m_kd;
}

glm::vec3 CVK::Material::getKdVec()
{
    return m_kdVec;
}

glm::vec3 *CVK::Material::getKdVecP()
{
    return &m_kdVec;
}

void CVK::Material::setKs(float ks)
{
    m_ks = ks;
}

float CVK::Material::getKs()
{
    return m_ks;
}

glm::vec3 CVK::Material::getKsVec()
{
    return m_ksVec;
}

glm::vec3 *CVK::Material::getKsVecP()
{
    return &m_ksVec;
}

void CVK::Material::setKt(float kt)
{
    m_kt = kt;
}

float CVK::Material::getKt()
{
    return m_kt;
}

void CVK::Material::setdiffColor(glm::vec3 col)
{
    m_diffColor = col;
}

glm::vec3 *CVK::Material::getdiffColor()
{
    return &m_diffColor;
}

void CVK::Material::setspecColor(glm::vec3 col)
{
    m_specColor = col;
}

glm::vec3 *CVK::Material::getspecColor()
{
    return &m_specColor;
}

void CVK::Material::setEmittedColor(const glm::vec3 &col)
{
    m_emitColor = col;
}

const glm::vec3 &CVK::Material::getEmittedColor() const
{
    return m_emitColor;
}

void CVK::Material::setShininess(float shininess)
{
    m_shininess = shininess;
}

float CVK::Material::getShininess()
{
    return m_shininess;
}

void CVK::Material::setIor(float ior)
{
    m_ior = ior;
}

float CVK::Material::getIor()
{
    return m_ior;
}

void CVK::Material::setIllum(int illum)
{
    m_illum = illum;
}

int CVK::Material::getIllum()
{
    return m_illum;
}

void CVK::Material::setAlpha(float alpha)
{
    m_alpha = alpha;
}

float CVK::Material::getAlpha()
{
    return m_alpha;
}

void CVK::Material::setTexture(TextureType type, const char *fileName)
{
    switch(type)
    {
        case COLOR_TEXTURE:
            m_colorTexturePath = fileName;
            if(m_colorTexture)
                m_colorTexture->load(fileName);
            else
                m_colorTexture = new Texture(fileName);
            break;
        case NORMAL_TEXTURE:
            m_normalTexturePath = fileName;
            if(m_normalTexture)
                m_normalTexture->load(fileName);
            else
                m_normalTexture = new Texture(fileName);
            break;
    }
}

void CVK::Material::setTexture(TextureType type, GLuint textureID)
{
    switch(type)
    {
        case COLOR_TEXTURE:
            if(m_colorTexture)
                m_colorTexture->setTexture(textureID);
            else
                m_colorTexture = new Texture(textureID);
            break;
        case NORMAL_TEXTURE:
            if(m_normalTexture)
                m_normalTexture->setTexture(textureID);
            else
                m_normalTexture = new Texture(textureID);
            break;
    }
}

bool CVK::Material::hasTexture(TextureType type)
{
    switch(type)
    {
        case COLOR_TEXTURE:
            return m_colorTexture != NULL;
        case NORMAL_TEXTURE:
            return m_normalTexture != NULL;
        default:
            return false;
    }
}

CVK::Texture *CVK::Material::getTexture(TextureType type)
{
    switch(type)
    {
        case COLOR_TEXTURE:
            return m_colorTexture;
        case NORMAL_TEXTURE:
            return m_normalTexture;
        default:
            return 0;
    }
}


