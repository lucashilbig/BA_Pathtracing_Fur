#include "SceneGraph.h"
#include "KIRK/Utils/jsonin.h"
#include "KIRK/Common/Light.h"
#include <experimental/filesystem>


bool KIRK::SceneGraph::endsWith(std::string const &fullString, std::string const &ending)
{
    if (fullString.length() >= ending.length())
    {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    }
    else
    {
        return false;
    }
}

std::shared_ptr<KIRK::SceneGraph> KIRK::SceneGraph::makeSceneGraph(const std::string &file_path)
{
	//Better this than crash around...
	if(!std::experimental::filesystem::exists(std::experimental::filesystem::path(file_path)))
	{
		LOG_ERRORs() << "File does not exist: " << file_path;

		LOG_INPUTs() << "Press [ENTER] to quit...";
		std::cin.ignore();
		exit(1);
	}

    if(file_path == "")
	{
        LOG_INFO("No file path specified. Creating a new empty scene..." );
        SceneGraph *sceneGraph = new SceneGraph;

        LOG_INFO("Adding default camera, light and environment...");
        sceneGraph->createDefaultCamera();
        sceneGraph->createDefaultLight();
        sceneGraph->createDefaultEnvironment();

        LOG_INFO("Finished creating the scene.");
        return std::shared_ptr<SceneGraph>(sceneGraph);
    }
	else if(endsWith(file_path, ".obj"))
    {
        LOG_INFO("Detected .obj file. Creating new scenegraph..." );
        SceneGraph *sceneGraph = new SceneGraph;
        sceneGraph->getRootNode()->attachChild(std::shared_ptr<SceneNode>(sceneGraph->importObj(file_path)));

        LOG_INFO("Imported file. Adding default camera, light and environment...");
        sceneGraph->createDefaultCamera();
        sceneGraph->createDefaultEnvironment();
        sceneGraph->createDefaultLight();

        LOG_INFO("Finished loading. Active camera is %.", 0);
        return std::shared_ptr<SceneGraph>(sceneGraph);
    }
    else if(endsWith(file_path, ".json"))
    {
		LOG_INFO("Detected .json file. Loading scenegraph...");
        return jsonio::makeSceneGraph(file_path);
    }
    else
    {
        LOG_ERROR("No supported file format. Exiting.");
        throw std::runtime_error("Can't construct SceneGraph.");
    }
}

KIRK::SceneGraph::~SceneGraph()
{
}

KIRK::SceneGraph::SceneGraph()
{
    initRoot();
}

void KIRK::SceneGraph::initRoot()
{
    m_root_node = std::shared_ptr<SceneNode>(new SceneNode);
    m_node_count = 1;
}

std::shared_ptr<KIRK::SceneNode> KIRK::SceneGraph::importObj(const std::string &file_path)
{
    LOG_INFO("Importing obj file...");
    Assimp::Importer importer;
    const aiScene* scene = importer.ReadFile(file_path, aiProcess_GenSmoothNormals
                                             | aiProcess_GenUVCoords);

    auto sub_root = std::shared_ptr<KIRK::SceneNode>(new SceneNode);
    m_node_count++;

	if (!scene)
	{
		LOG_ERROR("Scene failed to load: %", importer.GetErrorString());
		return sub_root;
	}

    LOG_DEBUG("----------------------------------------------" );
    LOG_INFO("Loading nodes..." );
	//Loading the raw file path to load the textures from.
	size_t found = file_path.find_last_of("\\/");
	std::string result_path = file_path.substr(0, found);

    assignAiNode(scene, scene->mRootNode, sub_root, result_path.c_str());
    LOG_DEBUG("Created % graph nodes.", m_node_count);
    LOG_DEBUG("----------------------------------------------" );

    return sub_root;
}

void KIRK::SceneGraph::assignAiNode(const aiScene* scene, aiNode* source_node, std::shared_ptr<KIRK::SceneNode> target_node, const std::string &scene_path)
{
    for(int i=0; i<source_node->mNumChildren; i++)
    {
        auto child = std::shared_ptr<KIRK::SceneNode>(new SceneNode);

        assignAiNode(scene, source_node->mChildren[i], child, scene_path);
        m_node_count++;

        //Set child-parent cross reference
        if(source_node->mChildren[i]->mNumMeshes>0)
		{
            std::shared_ptr<Mesh> mesh( createMeshFromAi(scene, source_node->mChildren[i], scene_path));
            if(!mesh)
			{
                LOG_ERROR("Found a non-triangulated Mesh! Skipping Node...");
                break;
            }

			child->setName("mesh " + std::to_string(m_node_count));
			//Set local mesh transformation to the child node, then you can simply reset it to a local center point.
			child->assignData(mesh, SceneNode::Type::MESH);
        }
        target_node->attachChild(std::shared_ptr<SceneNode>(child));
    }
}

KIRK::Mesh* KIRK::SceneGraph::createMeshFromAi(const aiScene* scene, aiNode* node, const std::string &scene_path)
{
    Mesh *mesh = new Mesh();
    int index_offset = 0;

    for(int m=0; m<node->mNumMeshes; m++)
	{
        aiMesh* source = scene->mMeshes[node->mMeshes[m]];
        aiMaterial* ai_material = scene->mMaterials[source->mMaterialIndex];

        // Let's initialize the materials coming from obj-files with just some default values.
        // Assimp is waaaaay too unpredictable to do anything useful there.
        // We can also more or less safely take the name from assimp.
        aiString name;
        ai_material->Get(AI_MATKEY_NAME, name);

        //Material does not exist yet. Create a new one and insert it into the map.
        KIRK::Material *kirk_material = new KIRK::Material(name.C_Str());

            aiColor3D color;
            ai_material->Get(AI_MATKEY_COLOR_DIFFUSE, color);
            kirk_material->m_diffuse.value = Color::RGBA{color.r, color.g, color.b, 1};

			aiString texture_path;
			if (ai_material->GetTextureCount(aiTextureType_DIFFUSE) > 0 && ai_material->GetTexture(aiTextureType_DIFFUSE, 0, &texture_path) == AI_SUCCESS) {
				kirk_material->m_diffuse.texture = std::make_shared<KIRK::Texture>((scene_path + ("/") + texture_path.C_Str()).c_str());
			}
			if (ai_material->GetTextureCount(aiTextureType_NORMALS) > 0 && ai_material->GetTexture(aiTextureType_NORMALS, 0, &texture_path) == AI_SUCCESS) {
				kirk_material->m_normal.texture = std::make_shared<KIRK::Texture>((scene_path + ("/") + texture_path.C_Str()).c_str());
			}

			ai_material->Get(AI_MATKEY_COLOR_EMISSIVE, color);
			kirk_material->m_emission.value = Color::RGBA{ color.r, color.g, color.b, 1 };

            kirk_material->name = name.C_Str();
            kirk_material->m_reflectivity.value = 0.2f;
            kirk_material->m_transparency.value = 0;
			kirk_material->m_roughness.value = 0.1f;
            kirk_material->m_ior = 1.56f;

		//Assign the material to the mesh
        mesh->m_materials.emplace_back(kirk_material);

		//Helper variable for when there are multiple sub-meshes which will be combined into one.
        index_offset = mesh->m_vertices.size();

        for(int f=0; f<source->mNumFaces; f++)
		{
            if(source->mFaces[f].mNumIndices<=2)
			{
                //We don't need those two-vertex-faces, let's throw them out.
                LOG_ERROR("Found weird face.");
                continue;
            }

            for(int vtx = 0; vtx < source->mFaces[f].mNumIndices-2; vtx++)
			{
                // We'll create a triangle fan out of all vertices contained in the face.
                // Kinda like a custom triangulation.

                unsigned int idx_a = source->mFaces[f].mIndices[0] + index_offset;
                unsigned int idx_b = source->mFaces[f].mIndices[vtx+1] + index_offset;
                unsigned int idx_c = source->mFaces[f].mIndices[vtx+2] + index_offset;

                Mesh::face face;
                face.material_index = mesh->m_materials.size()-1;
                face.vertex_index_a = idx_a;
                face.vertex_index_b = idx_b;
                face.vertex_index_c = idx_c;

                mesh->m_indices.push_back(idx_a);
                mesh->m_indices.push_back(idx_b);
                mesh->m_indices.push_back(idx_c);
                mesh->m_faces.push_back(face);
            }
        }

        for(int v=0; v<source->mNumVertices; v++)
		{
			//Getting the vertices from the previously loaded indices.
            Mesh::vertex vtx;
            vtx.position.x = source->mVertices[v].x;
            vtx.position.y = source->mVertices[v].y;
            vtx.position.z = source->mVertices[v].z;

            if(source->mNormals == NULL)
			{
                //Use default up-normal if assimp could not load one.
                vtx.normal = glm::vec3(0, 1, 0);
            }
			else
			{
                vtx.normal.x = source->mNormals[v].x;
                vtx.normal.y = source->mNormals[v].y;
                vtx.normal.z = source->mNormals[v].z;
            }

			//For now, we just use texCoord[0] here.
            if(source->HasTextureCoords(0))
			{
                vtx.texcoord.x = source->mTextureCoords[0][v].x;
                vtx.texcoord.y = source->mTextureCoords[0][v].y;
            }

			//Vertex ready for lift-off!!!
            mesh->m_vertices.push_back(vtx);
        }
    }

    LOG_DEBUG("Mesh loaded having");
    LOG_DEBUG("     % vertices.", mesh->m_vertices.size() );
    LOG_DEBUG("     % materials.", mesh->m_materials.size() );
    LOG_DEBUG("     % faces.", mesh->m_faces.size() );

    return mesh;
}

void KIRK::SceneGraph::createDefaultLight()
{
	//Just some default lights.
	Light *light_mid = new QuadLight({2.0f, 1.0f, 0.0f}, Color::WHITE, {-1.0f, -1.0f, 0.0f}, {1, 1}, 0, 0.01);
	Light *light_left = new PointLight({ -2.5f, 3.0f, 0.0f }, Color::RGBA(1.0f, 1.0f, 1.0f, 1.f), 0.8f, 0.002f);
	Light *light_right = new SpotLight({ 2.5f, 3.0f, 0.0f },  30*Color::RGBA(2.0f, 2.0f, 2.0f, 1.f), { -5.0f, -3.0f, 0.5f }, 0.5f, 11.25f, 10.f, 0, 0, 0.1);
	Light *sun = new SunLight(Color::WHITE, {10,-3,3}, 0.5);

    m_environment->setAmbientLight(Color::RGBA(0.2f, 0.2f, 0.2f, 1.f));

    std::shared_ptr<SceneNode> node_mid(new SceneNode());
	node_mid->setName("light_mid");
    node_mid -> assignData(std::shared_ptr<Light>(light_mid), SceneNode::Type::LIGHT);
	node_mid->transformWith(glm::rotate(glm::radians(45.0f), glm::vec3( -1.0f, -1.0f, 0.0f )));
    getRootNode()->attachChild(node_mid);

    std::shared_ptr<SceneNode> node_left(new SceneNode());
    node_left -> assignData(std::shared_ptr<Light>(light_left), SceneNode::Type::LIGHT);
	node_left->setName("node_left");
    getRootNode()->attachChild(node_left);

    std::shared_ptr<SceneNode> node_right(new SceneNode());
    node_right -> assignData(std::shared_ptr<Light>(light_right), SceneNode::Type::LIGHT);
	node_right->setName("node_right");
	node_right->transformWith(glm::translate(glm::vec3(2.0f, 0.0f, -2.0f)));
    getRootNode()->attachChild(node_right);

    std::shared_ptr<SceneNode> node_sun(new SceneNode());
    node_sun -> assignData(std::shared_ptr<Light>(sun), SceneNode::Type::LIGHT);
	node_sun->setName("node_sun");
    getRootNode()->attachChild(node_sun);
}

void KIRK::SceneGraph::createDefaultCamera()
{
    std::shared_ptr<Camera> camera(new Camera());

	//Just use some default position, look-at and up-vector.
    glm::vec3 position = glm::vec3(-3.f, 7.f, -9.f);
    glm::vec3 look_at = glm::normalize(-position);
    camera->setTransform(position, look_at, glm::vec3(0,1,0));

    std::shared_ptr<SceneNode> node(new SceneNode());
	node->setName("camera");
    node -> assignData(camera, SceneNode::Type::CAMERA);

    getRootNode()->attachChild(node);
    m_active_camera = camera;
}

void KIRK::SceneGraph::createDefaultEnvironment()
{
	//Use a light blue colored background as default-
    m_environment = std::unique_ptr<Environment>(new Environment());
    m_environment->setColor(Color::RGBA(0.7f, 0.9f, 1.f, 1.f));
}

void KIRK::SceneGraph::addFurFibersToAllMeshes(unsigned int num_fiber_verts, float fiber_radius)
{	
	for each (auto sceneNode in *this) {
		if (sceneNode->m_data_type == KIRK::SceneNode::MESH)
			std::dynamic_pointer_cast<KIRK::Mesh>(sceneNode->m_data_object)->addFurToFaces(num_fiber_verts, fiber_radius);
	}
}

////////////////////////////////////////////////////////////
///////////////
///////////////		GETTERS
///////////////
////////////////////////////////////////////////////////////

std::shared_ptr<KIRK::SceneNode> KIRK::SceneGraph::getRootNode() const
{
	return m_root_node;
}

std::shared_ptr<KIRK::Camera> KIRK::SceneGraph::getActiveCamera() const
{
    if (m_active_camera == nullptr)
        throw std::logic_error("There is no active camera!");
    return m_active_camera;
}

KIRK::Environment* KIRK::SceneGraph::getEnvironment() const
{
	return m_environment.get();
}

KIRK::SceneNodeIterator KIRK::SceneGraph::begin()
{
	return SceneNodeIterator(m_root_node);
}

KIRK::SceneNodeIterator KIRK::SceneGraph::end()
{
	return SceneNodeIterator(std::shared_ptr<SceneNode>(nullptr));
}









