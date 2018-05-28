#include "KIRK/Utils/Log.h"
#include "CVK_Node.h"
#include "CVK_ShaderMinimal.h"
//#include "CVK_Utils/MTLLoader.h"

CVK::Node::Node()
{
    init();
}

CVK::Node::Node(std::string name)
{
    init();
    m_name = name;
}

CVK::Node::Node(std::string name, std::string path)
{
    init();
    m_name = name;
    load(path);
}

CVK::Node::~Node()
{
}

void CVK::Node::init()
{
    m_name = "";
    m_geometry = 0;
    m_material = 0;
    m_modelMatrix = glm::mat4(1.0f);
    m_parent = 0;
}

void CVK::Node::load(std::string path)
{
    // load a scene with ASSIMP
    Assimp::Importer importer;
    const aiScene *scene = importer.ReadFile(path, aiProcess_GenSmoothNormals | aiProcess_Triangulate);
    if(scene)
        LOG_INFO("SUCCESS: Loaded Model from Path: \"%\"", path.c_str());
    else
    {
        LOG_ERROR("Loading failed from Path: \"%\" With error message: %\n", path.c_str(),
               importer.GetErrorString());
        return;
    }

    //get string of the path from the mtl file
    std::string pathMTL = path.replace(path.length() - 3, 3, "mtl");
    //cout << "[FilePath]: " << pathMTL << endl;
    //creating the loader
    //MTLLoader mtlLoader = MTLLoader(pathMTL);

    //reading the values of illums
    //mtlLoader.readIllum();


    std::vector<std::shared_ptr<CVK::Material>> materials;

    // load all materials in this file
    for(unsigned int i = 0; i < scene->mNumMaterials; i++)
    {

        /* information about .mtl values
         * For all three, if there is just a value for R, G and B are the same like R.
         * Kd / AI_MATKEY_COLOR_DIFFUSE		To specify the diffuse reflectivity of the current material -  range of 0.0 to 1.0
         * Ks / AI_MATKEY_COLOR_SPECULAR	To specify the specular reflectivity of the current material -  range of 0.0 to 1.0
         * Ka / AI_MATKEY_COLOR_AMBIENT		To specify the ambient reflectivity of the current material -  range of 0.0 to 1.0
         * Ke / AI_MATKEY_COLOR_EMISSIVE	To specify the emission of the current material -  range of 0 to ?100?
         * illum	Illumination Properties that are turned on - range from 0 to 10
         * 		0	Color on and Ambient off
                1	Color on and Ambient on
                2	Highlight on
                3	Reflection on and Ray trace on
                4	Transparency: Glass on
                    Reflection: Ray trace on
                5	Reflection: Fresnel on and Ray trace on
                6	Transparency: Refraction on
                    Reflection: Fresnel off and Ray trace on
                7	Transparency: Refraction on
                    Reflection: Fresnel on and Ray trace on
                8	Reflection on and Ray trace off
                9	Transparency: Glass on
                    Reflection: Ray trace off
                10	Casts shadows onto invisible surfaces
         * Ns / AI_MATKEY_SHININESS (multiplied with 4 ?why?) Specifies the specular exponent for the current material (high value results in a tight, concentrated highlight)
         * 							- normally ranges from 0 to 1000
         * Ni / AI_MATKEY_REFRACTI	Specifies the optical density for the surface (index of refraction) 1.0 light doesnt bend, smaller 1.0 bizarre results
         * 							- ranges from 0.001 to 10
         * d / AI_MATKEY_COLOR_TRANSPARENT transparency of the Object (1.0 is opaque) - ranges of 0.0 to 1.0
         * Tr = 1 - d
         * Tf	To specify the transmission filter of the current material - ranges of 0.0 to 1.0
         */

        aiMaterial *mat = scene->mMaterials[i];

        //readout of mtl names
        aiString objName;
        mat->Get(AI_MATKEY_NAME, objName);

        //get the corresponding illum value
       // int illum = mtlLoader.getIllum(objName.C_Str());
        int illum = 0;

        //get the corresponding Kd value
        aiColor3D diffColor(0.f, 0.f, 0.f);
        mat->Get(AI_MATKEY_COLOR_DIFFUSE, diffColor);

        //get the corresponding Ks value
        aiColor3D specColor(0.f, 0.f, 0.f);
        mat->Get(AI_MATKEY_COLOR_SPECULAR, specColor);

        aiColor3D emitColor(0.f, 0.f, 0.f);
        mat->Get(AI_MATKEY_COLOR_EMISSIVE, emitColor);

        //get the corresponding Ns value
        float shininess = 0.0f;
        mat->Get(AI_MATKEY_SHININESS, shininess);

        //get the corresponding Ni value
        float ior = 1.0f;
        mat->Get(AI_MATKEY_REFRACTI, ior);

        //get the corresponding d value (transparanty)
        float alpha = 1.0f;
        mat->Get(AI_MATKEY_COLOR_TRANSPARENT, alpha);

        glm::vec3 diffuseColor(diffColor.r, diffColor.g, diffColor.b);
        glm::vec3 specularColor(specColor.r, specColor.g, specColor.b);
        glm::vec3 emittedColor(emitColor.r, emitColor.g, emitColor.b);

        aiString fileName;  // filename
        mat->Get(AI_MATKEY_TEXTURE(aiTextureType_DIFFUSE, 0), fileName);

        std::string s = path.substr(0, path.rfind("/")) + "/";
        s += fileName.data;

        //materials.push_back(new CVK::Material(glm::vec3(0.0,1.0,0.0), glm::vec3(0.0,0.0,1.0), 10));
        if(fileName.length > 0)
            materials.push_back(
                    std::make_shared<CVK::Material>(const_cast<char *> ( s.c_str()), 1.f, 1.f, specularColor, shininess, ior, alpha,
                                      illum));
            //should set kd and ks!!
        else
        {
            materials.push_back(
				std::make_shared<CVK::Material>(diffuseColor, diffuseColor, specularColor, specularColor, shininess, ior, alpha,
                                      illum));
            //materials.push_back(new CVK::Material(0.7f, diffuseColor, 0.3f, specularColor, shininess));
        }

        materials.back()->setEmittedColor(emittedColor);
    }

    // load all meshes in this file
    for(unsigned int i = 0; i < scene->mNumMeshes; i++)
    {
        aiMesh *mesh = scene->mMeshes[i];

        std::shared_ptr<CVK::Geometry> geometry = std::make_shared<CVK::Geometry>();

        // load geometry information of the current mesh
        for(unsigned int vCount = 0; vCount < mesh->mNumVertices; vCount++)
        {
            // vertices
            aiVector3D v = mesh->mVertices[vCount];
            geometry->getVertices()->push_back(glm::vec4(v.x, v.y, v.z, 1.0f));

            // normals
            if(mesh->HasNormals())
            {
                v = mesh->mNormals[vCount];
                geometry->getNormals()->push_back(glm::vec3(v.x, v.y, v.z));
            }

            // texture coordinates
            if(mesh->HasTextureCoords(0))
            {
                v = mesh->mTextureCoords[0][vCount];
                geometry->getUVs()->push_back(glm::vec2(v.x, v.y));
            }
        }

        for(unsigned int fCount = 0; fCount < mesh->mNumFaces; fCount++)
        {
            aiFace f = mesh->mFaces[fCount];
            // index numbers
            for(unsigned int iCount = 0; iCount < f.mNumIndices; iCount++)
            {
                geometry->getIndex()->push_back(f.mIndices[iCount]);
            }
        }

        geometry->setMaterialIndex(mesh->mMaterialIndex);
        geometry->createBuffers();

        // new child node with the geometry and a connection to material
		std::shared_ptr<CVK::Node> child = std::make_shared<CVK::Node>();
        child->setGeometry(geometry);
        child->setMaterial(materials[mesh->mMaterialIndex]);
        addChild(child);
    }
}

void CVK::Node::render()
{
    glm::mat4 modelMatrix = glm::mat4(1.0f);
    CVK::Node *parentNode = m_parent;
    while(parentNode != 0)
    {
        modelMatrix = *parentNode->getModelMatrix() * modelMatrix;
        parentNode = parentNode->getParent();
    }
    render(modelMatrix);
}

void CVK::Node::render(glm::mat4 modelMatrix)
{
    // accumulate model matrix
    modelMatrix = modelMatrix * *getModelMatrix();

    // update shader
    CVK::State::getInstance()->getShader()->update(this);
    CVK::State::getInstance()->getShader()->updateModelMatrix(modelMatrix);

    // render geometry
    if(hasGeometry())
        m_geometry->render();

    // render children
    for(unsigned int i = 0; i < m_children.size(); i++)
    {
        m_children[i]->render(modelMatrix);
    }
}

CVK::Node *CVK::Node::find(std::string name)
{
    if(m_name == name)
        return this;
    for(unsigned int i = 0; i < m_children.size(); i++)
    {
        CVK::Node *result = m_children[i]->find(name);
        if(result != 0)
            return result;
    }
    return 0;
}

void CVK::Node::setName(std::string name)
{
    m_name = name;
}

std::string *CVK::Node::getName()
{
    return &m_name;
}

void CVK::Node::setGeometry(std::shared_ptr<CVK::Geometry> geometry)
{
    m_geometry = geometry;
}

std::shared_ptr<CVK::Geometry> CVK::Node::getGeometry()
{
    return m_geometry;
}

bool CVK::Node::hasGeometry()
{
    return (m_geometry != NULL);
}

void CVK::Node::setMaterial(std::shared_ptr<CVK::Material> material)
{
    m_material = material;
}

std::shared_ptr<CVK::Material> CVK::Node::getMaterial()
{
    return m_material;
}

bool CVK::Node::hasMaterial()
{
    return (m_material != NULL);
}

void CVK::Node::setModelMatrix(glm::mat4 modelMatrix)
{
    m_modelMatrix = modelMatrix;
}

glm::mat4 *CVK::Node::getModelMatrix()
{
    return &m_modelMatrix;
}

void CVK::Node::setParent(CVK::Node *node)
{
    m_parent = node;
}

CVK::Node *CVK::Node::getParent()
{
    return m_parent;
}

void CVK::Node::addChild(std::shared_ptr<CVK::Node> node)
{
    m_children.push_back(node);
    node->setParent(this);
}

std::vector<std::shared_ptr<CVK::Node>> &CVK::Node::getChildren()
{
    return m_children;
}
