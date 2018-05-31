#include "jsonin.h"

using json = nlohmann::json;

//well i needed to hack a bit :/
static std::string directory = "";

//===Helper Funktions===
float initfloat(json& input, std::string parameter, float def)
{
    try
    {
        return float(input.at(parameter));
    }
    catch (std::out_of_range)
    {
        std::string warning("No value given for " + parameter + " in " + input.value("name", ""));
        LOG_WARN(warning);
        LOG_WARN("Using default Parameter");
        return def;
    }
}

glm::vec2 initvec(json& input, std::string parameter, glm::vec2 def)
{
    try
    {
        auto val =input.at(parameter);
        return glm::vec2 (val[0], val[1]);
    }
    catch (std::out_of_range)
    {
        std::string warning("No value given for " + parameter + " in " + input.value("name", ""));
        LOG_WARN(warning);
        LOG_WARN("Using default Parameter");
        return def;
    }
}

glm::vec3 initvec(json& input, std::string parameter, glm::vec3 def)
{
    try
    {
        auto val =input.at(parameter);
        return glm::vec3 (val[0], val[1], val[2]);
    }
    catch (std::out_of_range)
    {
        std::string warning("No value given for " + parameter + " in " + input.value("name", ""));
        LOG_WARN(warning);
        LOG_WARN("Using default Parameter");
        return def;
    }
}

glm::vec4 initvec(json& input, std::string parameter, glm::vec4 def)
{
    try
    {
        auto val =input.at(parameter);
        return glm::vec4 (val[0], val[1], val[2], val[3]);
    }
    catch (std::out_of_range)
    {
        std::string warning("No value given for " + parameter + " in " + input.value("name", ""));
        LOG_WARN(warning);
        LOG_WARN("Using default Parameter");
        return def;
    }
}

glm::quat initquat(json& input, std::string parameter, glm::quat def)
{
    try
    {
        auto val =input.at(parameter);
        return glm::quat (val[0], val[1], val[2], val[3]);
    }
    catch (std::out_of_range)
    {
        std::string warning("No value given for " + parameter + " in " + input.value("name", ""));
        LOG_WARN(warning);
        LOG_WARN("Using default Parameter");
        return def;
    }
}
//=======================================================================

std::shared_ptr<KIRK::SceneGraph> jsonio::makeSceneGraph(std::string path)
{
    std::ifstream file(path, std::ios::in);

    if(!file.is_open())
    {
		LOG_ERROR("Can't open file: \"%\"", path);
        return nullptr;
    }

    std::size_t found = path.find_last_of("/\\");
    directory = path.erase(found, path.size());

	std::shared_ptr<KIRK::SceneGraph> scene(new KIRK::SceneGraph());

	// Let's use references to the according vectors and maps. Then we won't need any copy or move commands.
    std::vector<std::shared_ptr<KIRK::Camera>> cameras;
    std::vector<std::shared_ptr<KIRK::Light>> lights;
    std::vector<std::shared_ptr<KIRK::Mesh>> meshes;

    json jsonobject = json::parse(file);
    file.close();

    //This vector is localy needed to locate materials when the graph is build
    std::vector<std::shared_ptr<KIRK::Material>> mats;

    //This vector is needed to remember all nodes which are a child of the root
    //non-root children are linked to nodes referenced in this vec right after their creation,
    //but aren't referenced directly by the vec
    //at last all nodes referenced by the vec are linked to the root_node
    std::map<int, std::shared_ptr<KIRK::SceneNode>> nodes;

    for (json::iterator it = jsonobject.begin(); it != jsonobject.end(); ++it)
    {
      std::string key = (std::string)it.key();

      if (key == "Material")
      {
          for (json::iterator matit = it.value().begin(); matit != it.value().end(); ++matit)
          {
              mats << matit.value();
          }
      }

      else if (key == "Mesh")
      {
          for (json::iterator meshit = it.value().begin(); meshit != it.value().end(); ++meshit)
              meshes << meshit.value();
      }

      else if (key == "Node")
      {
          //iterates over the graph, but we do not work recursevly ... you know which word i mean
          for (json::iterator nodeit = it.value().begin(); nodeit != it.value().end(); ++nodeit)
          {
              if (nodeit.value()["object"] == "empty")
              {
                  std::shared_ptr<KIRK::SceneNode> empty_node (new KIRK::SceneNode());
                  try {
                      int parent_id = nodeit.value().at("parent_id");
                      nodes[parent_id]->attachChild(empty_node);
                      int id = nodeit.value().at("object_id");
                      nodes[id] = empty_node;

                      glm::vec3 trans_vec = initvec(nodeit.value(), "translation", glm::vec3(0.0,0.0,0.0));
                      glm::quat rot_quat = initquat(nodeit.value(), "rotation", glm::quat(1.0,0.0,0.0,0.0));
                      rot_quat = glm::normalize(rot_quat);
                      glm::vec3 scale_vec = initvec(nodeit.value(), "scale", glm::vec3(1.0,1.0,1.0));

                      glm::mat4 rotMatrix = glm::toMat4(rot_quat);
                      glm::mat4 transformation = rotMatrix * glm::scale((glm::translate(glm::mat4(1.0), trans_vec)), scale_vec);
                      empty_node->m_transform = transformation;
                  }
                  catch(std::out_of_range) {
                      int id = nodeit.value().at("object_id");
                      nodes[id] = empty_node;
                      scene->getRootNode()->attachChild(empty_node);

                      glm::vec3 trans_vec = initvec(nodeit.value(), "translation", glm::vec3(0.0,0.0,0.0));
                      glm::quat rot_quat = initquat(nodeit.value(), "rotation", glm::quat(1.0,0.0,0.0,0.0));
                      rot_quat = glm::normalize(rot_quat);
                      glm::vec3 scale_vec = initvec(nodeit.value(), "scale", glm::vec3(1.0,1.0,1.0));

                      glm::mat4 rotMatrix = glm::toMat4(rot_quat);
                      glm::mat4 transformation = rotMatrix * glm::scale((glm::translate(glm::mat4(1.0), trans_vec)), scale_vec);
                      empty_node->m_transform = transformation;
                  }
              }

              else if (nodeit.value()["object"] == "camera")
              {
                  cameras << nodeit.value();
                  std::shared_ptr<KIRK::SceneNode> cam_node (new KIRK::SceneNode);
                  cam_node->assignData(cameras.back(), KIRK::SceneNode::CAMERA);
                  try {
                      int parent_id = nodeit.value().at("parent_id");
                      nodes[parent_id]->attachChild(cam_node);
                      int id = nodeit.value().at("object_id");
                      nodes[id] = cam_node;
                      glm::vec3 trans_vec = initvec(nodeit.value(), "translation", glm::vec3(0.0,0.0,0.0));
                      glm::quat rot_quat = initquat(nodeit.value(), "rotation", glm::quat(1.0,0.0,0.0,0.0));
                      rot_quat = glm::normalize(rot_quat);
                      glm::vec3 scale_vec = initvec(nodeit.value(), "scale", glm::vec3(1.0,1.0,1.0));

                      glm::mat4 rotMatrix = glm::toMat4(rot_quat);
                      glm::mat4 transformation = rotMatrix * glm::scale((glm::translate(glm::mat4(1.0), trans_vec)), scale_vec);
                      cam_node->m_transform = transformation;
                      std::string name =  nodeit.value().value("name","no_name");
                      cam_node->m_name = name;
                  }
                  catch(std::out_of_range) {
                      int id = nodeit.value().at("object_id");
                      nodes[id] = cam_node;
                      glm::vec3 trans_vec = initvec(nodeit.value(), "translation", glm::vec3(0.0,0.0,0.0));
                      glm::quat rot_quat = initquat(nodeit.value(), "rotation", glm::quat(1.0,0.0,0.0,0.0));
                      rot_quat = glm::normalize(rot_quat);
                      glm::vec3 scale_vec = initvec(nodeit.value(), "scale", glm::vec3(1.0,1.0,1.0));

                      glm::mat4 rotMatrix = glm::toMat4(rot_quat);
                      glm::mat4 transformation = rotMatrix * glm::scale((glm::translate(glm::mat4(1.0), trans_vec)), scale_vec);
                      cam_node->m_transform = transformation;
                      scene->getRootNode()->attachChild(cam_node);
                      std::string name =  nodeit.value().value("name","no_name");
                      cam_node->m_name = name;
                  }
                  scene->setActiveCamera(cameras.back());
              }
              else if (nodeit.value()["object"] == "light")
              {
                  lights << nodeit.value();
                  std::shared_ptr<KIRK::SceneNode> light_node (new KIRK::SceneNode);
                  light_node->assignData(lights.back(), KIRK::SceneNode::LIGHT);
                  try
				  {
                      int parent_id = nodeit.value().at("parent_id");
                      nodes[parent_id]->attachChild(light_node);
                      int id = nodeit.value().at("object_id");
                      nodes[id] = light_node;
                      std::string name =  nodeit.value().value("name","no_name");
                      light_node->m_name = name;
                  }
                  catch(std::out_of_range)
				  {
                      int id = nodeit.value().at("object_id");
                      nodes[id] = light_node;
                      scene->getRootNode()->attachChild(light_node);
                      std::string name =  nodeit.value().value("name","no_name");
                      light_node->m_name = name;
                  }
              }
              else if ((nodeit.value()["object"] == "mesh"))
              {
                  std::shared_ptr<KIRK::SceneNode> mesh_node (new KIRK::SceneNode);
                  int mesh_id = nodeit.value().at("mesh_id");

                  mesh_node->assignData(meshes[mesh_id], KIRK::SceneNode::MESH);
                  try
				  {
                      int parent_id = nodeit.value().at("parent_id");
                      nodes[parent_id]->attachChild(mesh_node);
                      int id = nodeit.value().at("object_id");
                      nodes[id] = mesh_node;

                  }
                  catch(std::out_of_range)
				  {
                      int id = nodeit.value().at("object_id");
                      nodes[id] = mesh_node;
                      scene->getRootNode()->attachChild(mesh_node);
                  }

                  glm::vec3 trans_vec = initvec(nodeit.value(), "translation", glm::vec3(0.0,0.0,0.0));
                  glm::quat rot_quat = initquat(nodeit.value(), "rotation", glm::quat(1.0,0.0,0.0,0.0));
                  rot_quat = glm::normalize(rot_quat);
                  glm::vec3 scale_vec = initvec(nodeit.value(), "scale", glm::vec3(1.0,1.0,1.0));

                  glm::mat4 rotMatrix = glm::toMat4(rot_quat);
                  glm::mat4 transformation = rotMatrix * glm::scale((glm::translate(glm::mat4(1.0), trans_vec)), scale_vec);
                  mesh_node->m_transform = transformation;
                  std::string name =  nodeit.value().value("name","no_name");
                  mesh_node->m_name = name;


                  int mat_id;
                  try
                  {
                      mat_id = nodeit.value().at("material_id");
                      if (meshes[mesh_id]->m_materials.size() == 1)
                          meshes[mesh_id]->m_materials.erase(meshes[mesh_id]->m_materials.begin());
                      meshes[mesh_id]->m_materials.push_back(mats[mat_id]);
                  }
                  catch (std::domain_error)
                  {
                      json mat_ids =nodeit.value().at("material_id");
                      if (meshes[mesh_id]->m_materials.size() == 1)
                          meshes[mesh_id]->m_materials.erase(meshes[mesh_id]->m_materials.begin());
                      for (json::iterator i = mat_ids.begin(); i != mat_ids.end(); ++i)
                          meshes[mesh_id]->m_materials.push_back(mats[i.value()]);

                  }
              }
          }
      }
    }

    scene ->createDefaultEnvironment();

    try
    {
        json env = jsonobject.at("Environment");

        try
        {
            auto val = env.at("color");
            const KIRK::Color::RGBA color = KIRK::Color::RGBA{val[0],val[1],val[2],0.0};
            scene->getEnvironment() ->setColor(color);
        }
        catch (std::out_of_range)
        {
            if (env.at("kind") == "cube")
            {
                std::string realpath1 = directory + (env.value("posx", ""));
                const char* posx = realpath1.c_str();

                std::string realpath2 = directory + (env.value("posy", ""));
                const char* posy = realpath2.c_str();

                std::string realpath3 = directory + (env.value("posz", ""));
                const char* posz = realpath3.c_str();

                std::string realpath4 = directory + (env.value("negx", ""));
                const char* negx = realpath4.c_str();

                std::string realpath5 = directory + (env.value("negy", ""));
                const char* negy = realpath5.c_str();

                std::string realpath = directory + (env.value("negz", ""));
                const char* negz = realpath.c_str();

                scene->getEnvironment()->loadCubeMap(posx, posy, posz, negx, negy, negz);

            }
            else if (env.at("kind") == "sphere")
            {
                std::string relpath = "/";
                relpath.append(env.value("path", ""));
                std::string realpath = directory + relpath;
                const char* path = realpath.c_str();
				if (relpath != "/")
                    scene->getEnvironment()->loadSphereMap(path);
            }

        }
        try
        {
            auto val = env.at("light");
            scene->getEnvironment() ->setAmbientLight(KIRK::Color::RGBA{val[0],val[1],val[2],val[3]});
        }
        catch (std::out_of_range)
        {
        }


    }
    catch (std::out_of_range)
    {
    }

    if (scene->getRootNode()->m_children.size() == 1)
    {
        if (scene->getRootNode()->m_children[0]->m_data_type == KIRK::SceneNode::EMPTY)
        {
            std::shared_ptr<KIRK::SceneNode> new_root = scene->getRootNode()->m_children[0];
            scene->getRootNode() = new_root;
        }
    }

    return scene;
}

std::vector<std::shared_ptr<KIRK::Material>>& operator<<(std::vector<std::shared_ptr<KIRK::Material>>& materials, json input)
{
    try
    {
        std::string path = input.at("path");
        std::ifstream file(std::string(directory + path) , std::ios::in);

        if(!file.is_open())
        {
			LOG_WARN("Can't open file: \"%\"", path);
			LOG_WARN("Using default material.");
            input["name"] = "default";
        }
		else
        {
			input = json::parse(file);
			file.close();
        }
    }
    catch (std::out_of_range)
    { }
    std::string name = input.at("name");
    std::shared_ptr<KIRK::Material> material (new KIRK::Material(name));

    material->m_diffuse.value = initvec(input, "diffuse", material->m_diffuse.value);
	material->m_specular.value = initvec(input, "specular", material->m_specular.value);
	material->m_volume.value = initvec(input, "volume", material->m_volume.value);
	material->m_emission.value = initvec(input, "emission", material->m_emission.value);
	material->m_ior = initfloat(input, "ior", material->m_ior);
	material->m_transparency.value = initfloat(input, "transparency", material->m_transparency.value);
    material->m_reflectivity.value = initfloat(input, "reflectivity", material->m_reflectivity.value);
    material->m_roughness.value = initfloat(input, "roughness", material->m_roughness.value);

    //not tested
    std::string diff_path = input.value("diffuse_map", "");
    if (diff_path != "")
    {
        std::string path = (directory + diff_path);
        std::shared_ptr< KIRK::Texture> texture(new KIRK::Texture(path.c_str()));
        material->m_diffuse.texture = texture;
    }

    std::string spec_path = input.value("specular_map", "");
    if (spec_path != "")
    {
        std::string path = (directory + spec_path);
        std::shared_ptr< KIRK::Texture> texture(new KIRK::Texture(path.c_str()));
        material->m_specular.texture = texture;
    }


    std::string vol_path = input.value("volume_map", "");
    if (vol_path != "")
    {
        std::string path = (directory + vol_path);
        std::shared_ptr< KIRK::Texture> texture(new KIRK::Texture(path.c_str()));
        material->m_volume.texture = texture;
    }

    std::string emi_path =  input.value("emission_map", "");
    if (emi_path != "")
    {
        std::string path = (directory + emi_path);
        std::shared_ptr< KIRK::Texture> texture(new KIRK::Texture(path.c_str()));
        material->m_emission.texture = texture;
    }

    std::string norm_path =  input.value("normal_map", "");
    if (norm_path != "")
    {
        std::string path = (directory + norm_path);
        std::shared_ptr< KIRK::Texture> texture(new KIRK::Texture(path.c_str()));
        material->m_normal.texture = texture;
    }

    std::string bump_path =  input.value("bump_map", "");
    if (bump_path != "")
    {
        std::string path = (directory + bump_path);
        std::shared_ptr< KIRK::Texture> texture(new KIRK::Texture(path.c_str()));
        material->m_bump.texture = texture;
    }

    std::string tran_path =  input.value("transparency_map", "");
    if (tran_path != "")
    {
        std::string path = (directory + tran_path);
        std::shared_ptr< KIRK::Texture> texture(new KIRK::Texture(path.c_str()));
        material->m_transparency.texture = texture;
    }

    std::string roug_path =  input.value("roughness_map", "");
    if (roug_path != "")
    {
        std::string path = (directory + roug_path);
        std::shared_ptr< KIRK::Texture> texture(new KIRK::Texture(path.c_str()));
        material->m_roughness.texture = texture;
    }

    std::string bsdf_name =  input.value("bsdf", "");
    if (bsdf_name != "")
    {
        std::shared_ptr<KIRK::BSDF> bsdf =KIRK::BsdfFactory::getInstance().getBsdf(bsdf_name);
        material->m_bsdf = bsdf;
    }

    std::string shader_name =  input.value("shader", "");
    if (shader_name != "" && shader_name != "default")
    {
        std::shared_ptr<KIRK::Shader> shader =KIRK::ShaderFactory::getInstance().getShader(shader_name);
        material->m_shader = shader;
    }

    materials.push_back(material);

    return materials;
}

std::vector<std::shared_ptr<KIRK::Camera>>& operator<<(std::vector<std::shared_ptr<KIRK::Camera>>& cameras, json input)
{
    KIRK::Camera* camera (new KIRK::Camera());

    glm::vec3 position = initvec(input, "position", glm::vec3(0.0, 0.0, 0.0));
	glm::vec3 direction = initvec(input, "direction", glm::vec3(0));
	if (glm::length(direction) == 0) {
		direction = initvec(input, "look_at", glm::vec3(0.0, 0.0, 0.0)) - position;
	}
    glm::vec3 up_vector = initvec(input, "up_vector", glm::vec3(0.0, 0.0, 0.0));

    camera->setTransform(position, direction, up_vector);
	cameras.push_back(std::unique_ptr<KIRK::Camera>(camera));

    camera ->applyParameters();
    return cameras;
}

std::vector<std::shared_ptr<KIRK::Light>>& operator<<(std::vector<std::shared_ptr<KIRK::Light>>& lights, json input)
{

        KIRK::Light* light;

        if (input.at("kind")=="quad")
        {
            glm::vec3 position = initvec(input, "position",glm::vec3(0.0, 0.0, 0.0));
            glm::vec3 direction = initvec(input, "direction",glm::vec3(0.0, 0.0, 0.0));
            glm::vec2 size = initvec(input, "size",glm::vec2(0.0, 0.0));

            KIRK::Color::RGBA color = initvec(input, "color", glm::vec4(1.0,1.0,1.0,1.0));

            float constant = initfloat(input, "constant", 0.0);
            float lin_att = initfloat(input, "linear", 0.0);
            float quad_att = initfloat(input, "quadratic", 0.0);

            light = new KIRK::QuadLight(position, color, direction, size, constant, lin_att, quad_att);
        }

        else if (input.at("kind")=="point")
        {
            glm::vec3 position = initvec(input, "position",glm::vec3(0.0, 0.0, 0.0));
            KIRK::Color::RGBA color = initvec(input, "color", glm::vec4(1.0,1.0,1.0,1.0));
            float radius = initfloat(input, "radius", 0.0);
            float constant = initfloat(input, "constant", 0.0);
            float lin_att = initfloat(input, "linear", 0.0);
            float quad_att = initfloat(input, "quadratic", 0.0);
            light = new KIRK::PointLight(position, color, radius, constant, lin_att, quad_att);
        }

        else if (input.at("kind")=="spot")
        {
            glm::vec3 position = initvec(input, "position",glm::vec3(0.0, 0.0, 0.0));
            KIRK::Color::RGBA color = initvec(input, "color", glm::vec4(1.0,1.0,1.0,1.0));
            float radius = initfloat(input, "radius", 0.0);
            glm::vec3 direction = initvec(input, "direction",glm::vec3(0.0, 0.0, 0.0));
            float outer = initfloat(input, "outer", 0.0);

            float inner =  initfloat(input, "inner", 0.0);
            float constant = initfloat(input, "constant", 0.0);
            float lin_att = initfloat(input, "linear", 0.0);
            float quad_att = initfloat(input, "quadratic", 0.0);
            light = new KIRK::SpotLight(position, color, direction, radius, outer, inner, constant, lin_att, quad_att);
        }
		else if (input.at("kind")=="sun")
		{
            KIRK::Color::RGBA color = initvec(input, "color", glm::vec4(1.0,1.0,1.0,1.0));
            glm::vec3 direction = initvec(input, "direction",glm::vec3(0.0, 0.0, 0.0));
            float radius = initfloat(input, "radius", 0.0);
			light = new KIRK::SunLight(color, direction, radius);
		}

        lights.push_back(std::unique_ptr<KIRK::Light>(light));



    return lights;
}


std::vector<std::shared_ptr<KIRK::Mesh>>& operator<<(std::vector<std::shared_ptr<KIRK::Mesh>>& meshes, json input)
{
    //This code uses assimp to load the mesh from file and is copied from SceneGraph.cpp
    //No materials are imported from the mesh file.
    std::string realpath = directory + (input.value("path", ""));
    const char* path = realpath.c_str();
    std::size_t point = realpath.find_last_of(".");
    std::string end = realpath.substr(point, realpath.size());

    if (end == ".obj" || end == ".dae")
    {
        Assimp::Importer importer;

        std::string realpath = directory + (input.value("path", ""));
        const char* path = realpath.c_str();

		LOG_INFO("Importing mesh from %", path);
        const aiScene* scene = importer.ReadFile(path, aiProcess_GenSmoothNormals
                                                 | aiProcess_GenUVCoords);

        if (!scene)
            LOG_ERROR("Scene could not be loaded from %! Error: %", path, importer.GetErrorString());

        KIRK::Mesh* mesh = new KIRK::Mesh();

        int index_offset = mesh->m_vertices.size();

        aiNode* node = scene->mRootNode;

        for(int i=0; i<node->mNumChildren; i++)
        {
            aiNode* child_node = node->mChildren[i];

            for(int m=0; m<child_node->mNumMeshes; m++)
            {
                aiMesh* source = scene->mMeshes[child_node->mMeshes[m]];

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

                        KIRK::Mesh::face face;
                        face.vertex_index_a = idx_a;
                        face.vertex_index_b = idx_b;
                        face.vertex_index_c = idx_c;
                        if(!(end == ".dae"))
                            face.material_index = (source->mMaterialIndex)-1;
                        mesh->m_indices.push_back(idx_a);
                        mesh->m_indices.push_back(idx_b);
                        mesh->m_indices.push_back(idx_c);
                        mesh->m_faces.push_back(std::make_shared<KIRK::Mesh::face>(face));
                    }
                }

                for(int v=0; v<source->mNumVertices; v++)
                {
                    KIRK::Mesh::vertex vtx;
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

                    if(source->HasTextureCoords(0))
                    {
                        vtx.texcoord.x = source->mTextureCoords[0][v].x;
                        vtx.texcoord.y = source->mTextureCoords[0][v].y;
                    }

                    mesh->m_vertices.push_back(vtx);
                }
            }

        }

        LOG_DEBUG("Vertex count: %", mesh->m_vertices.size());

        meshes.push_back(std::unique_ptr<KIRK::Mesh>(mesh));
    }
    else if (end == ".b3df")
        meshes.push_back(KIRK::loadBinaryGeometry(path));
    return meshes;
}
