#ifndef KIRK_JSONIO_H
#define KIRK_JSONIO_H

#include <string>
#include <map>
#include <fstream>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <vector>
#include <iterator>
#include <memory>
#include "json.h"
#include "KIRK/Common/SceneGraph.h"
#include "KIRK/Common/SceneNode.h"
#include "KIRK/Common/Material.h"
#include "KIRK/Common/Camera.h"
#include "KIRK/Common/Light.h"
#include "KIRK/Utils/BinaryModelUtils.h"
#include "KIRK/Common/Shading/BsdfFactory.h"
#include "KIRK/Common/Shading/ShaderFactory.h"
#include <assimp/scene.h>
#include <assimp/Importer.hpp>
#include <assimp/postprocess.h>


using json = nlohmann::json;

namespace jsonio
{
    /**
    * @brief Contructs a Scenegraph described in a .json file
    * @param path to the .json file
    * @return shared_ptr to a scenegraph
    */
    std::shared_ptr<KIRK::SceneGraph> makeSceneGraph(std::string path);
}

std::vector<std::shared_ptr<KIRK::Material>>& operator<<(std::vector<std::shared_ptr<KIRK::Material>>& materials, json input);

std::vector<std::shared_ptr<KIRK::Camera>>& operator<<(std::vector<std::shared_ptr<KIRK::Camera>>& cameras, json input);

std::vector<std::shared_ptr<KIRK::Light>>& operator<<(std::vector<std::shared_ptr<KIRK::Light>>& lights, json input);

std::vector<std::shared_ptr<KIRK::Mesh>>& operator<<(std::vector<std::shared_ptr<KIRK::Mesh>>& meshs, json input);

#endif  //KIRK_JSONIO_H
