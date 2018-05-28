#include "jsonout.h"
#include <glm/gtx/matrix_decompose.hpp>

using json = nlohmann::json;

void graphToList(std::shared_ptr<KIRK::SceneNode> node, std::vector<std::shared_ptr<KIRK::SceneNode>>& outlist)
{
    //if(node->m_data_type != KIRK::SceneNode::EMPTY)
        outlist.push_back(node);
    for(std::shared_ptr<KIRK::SceneNode> child : node->m_children)
        graphToList(child, outlist);
}

void jsonio::writeScene(std::shared_ptr<KIRK::SceneGraph> scene, std::string path)
{
    if(!std::experimental::filesystem::exists(path))
        std::experimental::filesystem::create_directory(path);

    if(!std::experimental::filesystem::exists(path + "/meshes"))
        std::experimental::filesystem::create_directory(path + "/meshes");

    if(!std::experimental::filesystem::exists(path + "/textures"))
        std::experimental::filesystem::create_directory(path + "/textures");

    if(!std::experimental::filesystem::exists(path + "/materials"))
        std::experimental::filesystem::create_directory(path + "/materials");

    json materials;
    json meshes;
    json nodes;

    std::vector<std::shared_ptr<KIRK::SceneNode>> outlist;
    graphToList(scene->getRootNode(), outlist);
    int i = 0;
    int ci = 0;
    int mi = 0;
    int mati = 0;
    int li = 0;

    //This map assigns a node to his id, we need this for parent child relationships.
    std::map<KIRK::SceneNode*, int> indexmap;
    for (int i = 0; i < outlist.size(); i++)
        indexmap[outlist.at(i).get()] = i;

    for(std::shared_ptr<KIRK::SceneNode> node : outlist)
    {
        if(node->m_data_type == KIRK::SceneNode::EMPTY)
        {
            json out_empty;
            out_empty["object"] = "empty";
            out_empty["object_id"] = i;
            glm::mat4 transform = node->m_transform;
            glm::vec3 scale;
            glm::quat rotation;
            glm::vec3 translation;
            glm::vec3 skew;
            glm::vec4 perspective;
            glm::decompose(transform, scale, rotation, translation, skew, perspective);
            out_empty["translation"] = {translation.x, translation.y, translation.z};
            out_empty["rotation"] = {rotation.w, rotation.x, rotation.y, rotation.z};
            out_empty["scale"] = {scale.x, scale.y, scale.z};
            i++;
            KIRK::SceneNode* parent =  node->m_parent.lock().get();
            if (!(parent == scene->getRootNode().get() || parent == nullptr))
            {
                out_empty["parent_id"] = indexmap[parent];
            }
            nodes.push_back(out_empty);
        }

        else if(node->m_data_type == KIRK::SceneNode::CAMERA)
        {
            KIRK::Camera* camera = (KIRK::Camera*)node->m_data_object.get();
            json out_camera;
            out_camera["object"] = "camera";
            out_camera["object_id"] = i;
            out_camera["camera_id"] = ci;
            //This asumes for the moment that camera is child for root, since parent->child relationships aren't supported at the moment
            out_camera["position"] = {camera->getLocalPosition().x, camera->getLocalPosition().y, camera->getLocalPosition().z};
            out_camera["look_at"] = {camera->getLocalLookAt().x, camera->getLocalLookAt().y, camera->getLocalLookAt().z};
            out_camera["up_vector"] = {camera->getLocalUp().x, camera->getLocalUp().y, camera->getLocalUp().z};
            out_camera["name"] = node->m_name;

            KIRK::SceneNode* parent =  node->m_parent.lock().get();
            if (!(parent == scene->getRootNode().get() || parent == nullptr))
            {
                out_camera["parent_id"] = indexmap[parent];
            }
            nodes.push_back(out_camera);
            i++;
            ci++;

        }

        else if(node->m_data_type == KIRK::SceneNode::LIGHT)
        {
            KIRK::Light* light = (KIRK::Light*)node->m_data_object.get();
            json out_light;
            out_light["object"] = "light";
            out_light["object_id"] = i;
			out_light["color"] = { light->m_color.x, light->m_color.y, light->m_color.z, light->m_color.w };
            out_light["name"] = node->m_name;

			if (typeid(*light) != typeid(KIRK::PointLight)) {
				out_light["direction"] = { light->m_direction.x, light->m_direction.y, light->m_direction.z };
			}
			if (typeid(*light) != typeid(KIRK::QuadLight)) {
				out_light["radius"] = light->m_radius;
			}
			if (typeid(*light) != typeid(KIRK::SunLight)) {
				out_light["position"] = { light->m_position.x, light->m_position.y, light->m_position.z };
				out_light["constant"] = light->m_const;
				out_light["linear"] = light->m_lin_att;
				out_light["quadratic"] = light->m_lin_att;
			}
            out_light["light_id"] = li;

            if(typeid(*(light)) == typeid(KIRK::PointLight))
            {
                out_light["kind"] = "point";
            }
            else if(typeid(*(light)) == typeid(KIRK::QuadLight))
            {
                out_light["kind"] = "quad";
				KIRK::QuadLight* lite = dynamic_cast<KIRK::QuadLight*>(light);
				out_light["size"] = { lite->m_size[0], lite->m_size[1] };
            }
            else if(typeid(*(light)) == typeid(KIRK::SpotLight))
            {
                out_light["kind"] = "spot";
                out_light["outer"] = ((KIRK::SpotLight*)light)->m_outer_angle;
                out_light["inner"] = ((KIRK::SpotLight*)light)->m_inner_angle;
            }
            else if(typeid(*(light)) == typeid(KIRK::SunLight))
            {
                out_light["kind"] = "sun";
            }

            if (!(node->m_parent.lock().get() == scene->getRootNode().get() || node->m_parent.lock().get() == nullptr))
                out_light["parent_id"] = indexmap[node->m_parent.lock().get()];
            nodes.push_back(out_light);
            i++;
            li++;

        }

        else if (node->m_data_type == KIRK::SceneNode::MESH)
        {
            KIRK::Mesh* mesh = (KIRK::Mesh*)node->m_data_object.get();
            json out_mesh;
            out_mesh["mesh_id"] = mi;
            std::string mesh_path = path + "/meshes/";
            out_mesh["path"] =  "/meshes/mesh" + std::to_string(mi)+ ".b3df";
            KIRK::writeBinaryGeometry(mesh_path + "mesh" + std::to_string(mi)+ ".b3df",mesh, false);
            meshes.push_back(out_mesh);
            out_mesh["name"] = node->m_name;

            json out_meshnode;
            for (auto material : mesh->m_materials)
            {
                json out_material;
                out_material["material_id"] = mati;

                if(!std::experimental::filesystem::exists(path + "/materials/" + material->name + ".json"))
                {
                    json material_json;
                    material_json["name"] = material->name;
                    material_json["diffuse"] = {material->m_diffuse.value.x, material->m_diffuse.value.y, material->m_diffuse.value.z,
                            material->m_diffuse.value.w};
                    material_json["specular"] = {material->m_specular.value.x, material->m_specular.value.y, material->m_specular.value.z,
                            material->m_specular.value.w};
                    material_json["emission"] = {material->m_emission.value.x, material->m_emission.value.y, material->m_emission.value.z,
                            material->m_emission.value.w};
                    material_json["volume"] = { material->m_volume.value.x, material->m_volume.value.y, material->m_volume.value.z,
                            material->m_volume.value.w };
                    material_json["ior"] = material->m_ior;
                    material_json["transparency"] = material->m_transparency.value;
                    material_json["reflectivity"] = material->m_reflectivity.value;
                    material_json["roughness"] = material->m_roughness.value;
                    material_json["bsdf"] = material->m_bsdf->getName();
                    material_json["shader"] = material->m_shader->m_name;

                    if(material ->m_diffuse.texture != nullptr)
                    {
                        material_json["diffuse_map"] = "../textures/" + material->name + "_diffuse.png";
                        material->m_diffuse.texture->saveTo((path + "/textures/" + material->name + "_diffuse.png").c_str());
                    }
                    if(material ->m_bump.texture != nullptr)
                    {
                        material_json["bump_map"] = "../textures/" + material->name + "_bump.png";
                        material->m_bump.texture->saveTo((path + "/textures/" + material->name + "_bump.png").c_str());
                    }
                    if(material ->m_emission.texture != nullptr)
                    {
                        material_json["emission_map"] = "../textures/" + material->name + "_emission.png";
                        material->m_emission.texture->saveTo((path + "/textures/" + material->name + "_emission.png").c_str());
                    }
                    if(material ->m_roughness.texture != nullptr)
                    {
                        material_json["roughness_map"] = "../textures/" + material->name + "_roughness.png";
                        material->m_roughness.texture->saveTo((path + "/textures/" + material->name + "_roughness.png").c_str());
                    }
                    if(material ->m_specular.texture != nullptr)
                    {
                        material_json["specular_map"] = "../textures/" + material->name + "_specular.png";
                        material->m_specular.texture->saveTo((path + "/textures/" + material->name + "_specular.png").c_str());
                    }
                    if(material ->m_transparency.texture != nullptr)
                    {
                        material_json["transparency_map"] = "../textures/" + material->name + "_transparency.png";
                        material->m_transparency.texture->saveTo((path + "/textures/" + material->name + "_transparency.png").c_str());
                    }
                    if(material ->m_volume.texture != nullptr)
                    {
                        material_json["volume_map"] = "../textures/" + material->name + "_volume.png";
                        material->m_volume.texture->saveTo((path + "/textures/" + material->name + "_volume.png").c_str());
                    }

                    std::string out_path = path + "/materials/" + material->name + ".json";
                    std::ofstream o(out_path);
                    o << std::setw(4) << material_json << std::endl;
                }
                out_material["path"] = "/materials/" + material->name + ".json";
                materials.push_back(out_material);
                out_meshnode["material_id"].push_back(mati);
                mati++;
            }

            out_meshnode["object"] = "mesh";
            out_meshnode["object_id"] = i;
            out_meshnode["mesh_id"] = mi;
            glm::mat4 transform = node->m_transform;
            glm::vec3 scale;
            glm::quat rotation;
            glm::vec3 translation;
            glm::vec3 skew;
            glm::vec4 perspective;
            glm::decompose(transform, scale, rotation, translation, skew, perspective);
            out_meshnode["translation"] = {translation.x, translation.y, translation.z};
            out_meshnode["rotation"] = {rotation.w, rotation.x, rotation.y, rotation.z};
            out_meshnode["scale"] = {scale.x, scale.y, scale.z};


            if (!(node->m_parent.lock().get() == scene->getRootNode().get() || node->m_parent.lock().get() == nullptr))
                out_meshnode["parent_id"] = indexmap[node->m_parent.lock().get()];
            nodes.push_back(out_meshnode);

            mi++;
            i++;
        }
    }

    json out_env;
    if (scene->getEnvironment()->m_type == KIRK::Environment::COLOR)
    {
        out_env["color"] = {scene->getEnvironment()->getBackgroundColor().x,
                scene->getEnvironment()->getBackgroundColor().y,
                scene->getEnvironment()->getBackgroundColor().z,
                scene->getEnvironment()->getBackgroundColor().w};
    }
    else if (scene->getEnvironment()->m_type == KIRK::Environment::SPHERE_MAP)
    {
        out_env["kind"] = "sphere";
        out_env["path"] = "/textures/environment.png";
        scene->getEnvironment()->m_spheremap_texture->saveTo((path + "/textures/environment.png").c_str());
    }
    else
    {
        out_env["kind"] = "cube";
        out_env["posx"] = "/textures/environment_posx.png";
        out_env["posy"] = "/textures/environment_posy.png";
        out_env["posz"] = "/textures/environment_posz.png";
        out_env["negx"] = "/textures/environment_negx.png";
        out_env["negy"] = "/textures/environment_negy.png";
        out_env["negz"] = "/textures/environment_negz.png";
        scene->getEnvironment()->m_cubemap_textures[0]->saveTo((path + "/textures/environment_posx.png").c_str());
        scene->getEnvironment()->m_cubemap_textures[1]->saveTo((path + "/textures/environment_posy.png").c_str());
        scene->getEnvironment()->m_cubemap_textures[2]->saveTo((path + "/textures/environment_posz.png").c_str());
        scene->getEnvironment()->m_cubemap_textures[3]->saveTo((path + "/textures/environment_negx.png").c_str());
        scene->getEnvironment()->m_cubemap_textures[4]->saveTo((path + "/textures/environment_negy.png").c_str());
        scene->getEnvironment()->m_cubemap_textures[5]->saveTo((path + "/textures/environment_negz.png").c_str());
    }

    out_env["light"] = {scene->getEnvironment()->getAmbientLight().x,
            scene->getEnvironment()->getAmbientLight().y,
            scene->getEnvironment()->getAmbientLight().z,
            scene->getEnvironment()->getAmbientLight().w};


    json out_json;
    out_json["Material"] = materials;
    out_json["Mesh"]= meshes;
    out_json["Node"] = nodes;
    out_json["Environment"] = out_env;

    std::string out_path;
    out_path.append(path + "/scene.json");
    std::ofstream o(out_path);
    o << std::setw(4) << out_json << std::endl;

}
