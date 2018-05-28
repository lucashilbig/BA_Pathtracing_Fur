#ifndef KIRK_JSONOUT_H
#define KIRK_JSONOUT_H

#include <string>
#include <map>
#include <fstream>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <vector>
#include <iterator>
#include <memory>
#include <experimental/filesystem>
#include "json.h"
#include "KIRK/Common/SceneGraph.h"
#include "KIRK/Common/SceneNode.h"
#include "KIRK/Common/Material.h"
#include "KIRK/Common/Camera.h"
#include "KIRK/Common/Light.h"
#include "KIRK/Common/Shading/Shader.h"
#include "KIRK/Utils/BinaryModelUtils.h"


using json = nlohmann::json;

namespace jsonio
{
    /**
    * @brief Savesthe scene in a .json file
    * @param shared_ptr to the scenegraph
    * @param name of the .json file
    */
    void writeScene(std::shared_ptr<KIRK::SceneGraph> scene, std::string path);
}

#endif  //KIRK_JSONOUT_H
