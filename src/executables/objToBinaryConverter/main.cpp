#include <iostream>
#include "KIRK/Common/SceneGraph.h"
#include <iostream>
#include "KIRK/Common/SceneGraph.h"
#include "KIRK/Common/SceneNode.h"
#include "KIRK/Utils/BinaryModelUtils.h"
#include "KIRK/Utils/ArgParser.h"
#include <KIRK/Utils/jsonout.h>
#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		std::cout << "No directory given." << std::endl;
	}
	else if (argc < 3)
	{
		std::cout << "No destination path chosen." << std::endl;
	}
	else
	{
		std::string cur_path; //changes each iteration
		std::string dir_path = std::string(RESOURCES_PATH) + argv[1]; //path containing all the files for conversion
		std::string des_path = std::string(RESOURCES_PATH) + argv[2]; //destination path for written scenes
		fs::create_directory(des_path);
		int scene_counter = 0;
		for (auto& p : fs::recursive_directory_iterator(dir_path))
		{
			cur_path = p.path().string();
			std::cout << cur_path << std::endl;
			if (KIRK::SceneGraph::endsWith(cur_path, ".obj") || KIRK::SceneGraph::endsWith(cur_path, ".json"))
			{
				scene_counter++;
				std::shared_ptr<KIRK::SceneGraph> scene = KIRK::SceneGraph::makeSceneGraph(cur_path);
				jsonio::writeScene(scene, des_path + "/scene" + std::to_string(scene_counter));
			}
		}
	}

	return 0;
}
