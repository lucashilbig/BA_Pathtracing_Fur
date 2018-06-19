#include "CVK/CVK_Utils/CVK_ShaderUtils.h"
#include "CVK/CVK_Utils/CVK_ConverterUtils.h"

#include "KIRK/Common/SceneGraph.h"

#include "KIRK/Utils/Window.h"
#include "KIRK/Utils/ArgParser.h"
#include "KIRK/Utils/jsonout.h"
#include "KIRK/Utils/Gui/Gui.h"
#include "KIRK/Utils/Gui/GuiScene.h"
#include "KIRK/Utils/Gui/GuiFileDialog.h"
#include "KIRK/Utils/Gui/GuiMenu.h"

#include "KIRK/CPU/CPU_Raytracer/CPU_PathTracer.h"

#include "KIRK/CPU/CPU_Datastructures/CPU_BVH.h"
#include "KIRK/CPU/CPU_Datastructures/CPU_NoDataStructure.h"

#include <externals/ImGui/imgui_impl_glfw_gl3.h>

//////////////////////////////////////////////
//
//  Member Variables
//
//////////////////////////////////////////////

//booleans used for settings
bool m_is_raytracer_active = false;
bool m_show_render_texture = false;

//platform which is used to render image
static int m_platform_used = 0; //value 0 for OpenGL ; 1 for CPU 

//ID of datastructure which is used. Default: 1 BVH
static int m_cpu_ds_id = 1;

//create cpu raytracer objects
KIRK::CPU::CPU_Raytracer *m_active_cpu_raytracer;
auto cpu_pathtracer = new KIRK::CPU::PathTracer();

//textures we use for rendering the images
std::shared_ptr<KIRK::Texture> render_texture_cpu;

//wrapper for the camera
CVK::CVKCameraSynchronizer *camera_wrapper;

//////////////////////////////////////////////
//
//  Member Functions
//
//////////////////////////////////////////////

CVK::Trackball *getTrackballCamera(int width, int height)
{
	CVK::Perspective *projection = new CVK::Perspective(glm::radians(60.f), width / (float)height, 0.0001f, 50.f);
	CVK::Trackball *cam_trackball = new CVK::Trackball(width, height, projection);
	glm::vec3 vec = glm::vec3(0.0f, 0.0f, 0.0f);
	cam_trackball->setCenter(&vec);
	cam_trackball->setRadius(5);
	CVK::State::getInstance()->setCamera(cam_trackball);
	return cam_trackball;
}

CVK::FreeFlight *getFreeFlightCamera(int width, int height)
{
	CVK::Perspective *projection = new CVK::Perspective(glm::radians(60.f), width / static_cast<float>(height), 0.1f, 50.f);
	CVK::FreeFlight *cam_freeflight = new CVK::FreeFlight(width, height, projection);
	CVK::State::getInstance()->setCamera(cam_freeflight);
	return cam_freeflight;
}

void setActiveCPURaytracer(KIRK::CPU::CPU_Raytracer *active_cpu_raytracer)
{
	if (!m_is_raytracer_active && m_active_cpu_raytracer != active_cpu_raytracer)
	{
		//set active cpu raytracer
		m_active_cpu_raytracer = active_cpu_raytracer;
		m_platform_used = 1;

		//reset the pathtracer
		if (active_cpu_raytracer == cpu_pathtracer)
			static_cast<KIRK::CPU::PathTracer*>(m_active_cpu_raytracer)->reset();
	}
}



/** @brief Switches the used datastructure for all raytracers of the current used platform. 0 for NoDatastructure, 1 for BVH, 2 KDTree, 3 Octree.
* !! raytracers have to be initialised with a scene before calling this func !!
*@param datastructure_id ID of the Datastructure which should be used
*/
void setDatastructure(const int datastructure_id)
{
	if (!m_is_raytracer_active)
	{
		//CPU Raytracer is used
		if (m_platform_used == 1 && m_cpu_ds_id != datastructure_id)
		{
			m_cpu_ds_id = datastructure_id;
			switch (datastructure_id)
			{
			case 0://CPU::NoDatastructure
				LOG_INFO("Switching CPU Datastructure to: NoDatastructure ...");
				m_active_cpu_raytracer->getScene().setDataStructure(std::unique_ptr<KIRK::CPU::CPU_DataStructure>(std::make_unique<KIRK::CPU::NoDataStructure>()), true);
				LOG_INFO("... done.");
				break;
			case 1://CPU::BVH which is Bounding Volume Hierarchy
				LOG_INFO("Switching CPU Datastructure to: Bounding Volume Hierarchy ...");
				m_active_cpu_raytracer->getScene().setDataStructure(std::unique_ptr<KIRK::CPU::CPU_DataStructure>(std::make_unique<KIRK::CPU::BVH>()), true);
				LOG_INFO("... done.");
				break;
			case 2://CPU::KDTree
				LOG_INFO("Switching CPU Datastructure to: KD Tree ...");
				m_active_cpu_raytracer->getScene().setDataStructure(std::unique_ptr<KIRK::CPU::CPU_DataStructure>(std::make_unique<KIRK::CPU::KDTree>()), true);
				LOG_INFO("... done.");
				break;
			default://default datastructure is BVH
				LOG_INFO("Switching CPU Datastructure to: Bounding Volume Hierarchy ...");
				m_active_cpu_raytracer->getScene().setDataStructure(std::unique_ptr<KIRK::CPU::CPU_DataStructure>(std::make_unique<KIRK::CPU::BVH>()), true);
				LOG_INFO("... done.");
				break;
			}
		}
	}
}

void saveRenderTexture()
{

	if (m_platform_used == 1)//cpu mode is active
	{
		render_texture_cpu->glDownload();
		render_texture_cpu->saveTo(SCREENSHOTS_PATH "/render_out_cpu.png");
	}
	else//openGL mode is active
	{
		LOG_WARN("Can't save render texture in OpenGL mode! \n Switch to CPU or GPU mode to save the current CPU or GPU render texture.");
	}
}

void resizeCallback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
	camera_wrapper->update(window);
	render_texture_cpu->resize(width, height);
	cpu_pathtracer->windowResizeCallback();
}

void charCallback(GLFWwindow *window, unsigned int key)
{
	switch (key)
	{
	case 'r'://render new Image with raytracer
		if (m_platform_used == 1 || m_platform_used == 2)
		{
			LOG_INFO("\n Rendering new Image \n -------------------");
			m_is_raytracer_active = true;
		}
		else
		{
			LOG_WARN("You can't render an Raytracing image in OpenGL mode!");
		}
		break;
	case 't'://stop rendering
		m_is_raytracer_active = false;
		LOG_INFO("Stopped Rendering");
		break;
	case '0'://switch to OpenGL mode
		m_platform_used = 0;
		break;
	case '1'://switch to CPU mode
		m_platform_used = 1;
		break;
	default:
		break;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////
//
//  MAIN:
//  ARGUMENTS:
//		-s <scene_path>		Load a scene from file system with relative path  (Prioritized)
//		-r <scene_res>		Load a scene from RESOURCES_PATH
//		-w <img_width_px>	Set a render image width
//		-h <img_height_px>	Set a render image height
//
////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
	//////////////////////////////////////////////
	//
	//  Loading and creating a Scene
	//
	//////////////////////////////////////////////

	//Parse arguments
	auto arg_map = KIRK::ArgParser::toMap(argc, argv);

	//For a little more convenience, we can also input a file path for the obj to load and render.
	//Loading arguments is just as easy as to have a key char and a default value for when the argument is not set.
	//CAUTION: Cannot use any const here. So no const char* allowed.
	//Not important to have optimized const char*'s here for appending RESOURCES_PATH, all this is kind of temporary.
	std::string file = std::string(SCENES_PATH) + arg_map.get('r', std::string("/Fur_SkinPatch/scene_lowLight.json"));
	std::string relative_path = arg_map.get('s', std::string(""));
	unsigned int image_width = arg_map.get('w', 1280);
	unsigned int image_height = arg_map.get('h', 720);

	//If both, a relative and a resources path, are given, use the relative source path instead of the resources path.
	if (relative_path != "")
	{
		file = relative_path;
	}

	//Load up scene.
	std::shared_ptr<KIRK::SceneGraph> scene(KIRK::SceneGraph::makeSceneGraph(file));

	//Initializing environment map (cubeMap here)
	scene->getEnvironment()->loadCubeMap(RESOURCES_PATH "/cubeMap/posx.jpg",
		RESOURCES_PATH "/cubeMap/posy.jpg",
		RESOURCES_PATH "/cubeMap/posz.jpg",
		RESOURCES_PATH "/cubeMap/negx.jpg",
		RESOURCES_PATH "/cubeMap/negy.jpg",
		RESOURCES_PATH "/cubeMap/negz.jpg");

	//Apply fur on every triangle in the scene. For testing purpose.
	scene->addFurFibersToAllMeshes(15, 0.004f);

	//////////////////////////////////////////////
	//
	//  General stuff
	//
	//////////////////////////////////////////////

	//We need a render texture to render to
	render_texture_cpu = std::make_shared<KIRK::Texture>(image_width, image_height, int(KIRK::Texture::TEX_RGB_UB));

	//For the cvk camera on the one side to be up to date on the kirk camera's parameters and for the
	//kirk camera's position to be up to date with the user set cvk camera, we have to use this wrapper to keep all tidy
	camera_wrapper = new CVK::CVKCameraSynchronizer(scene, std::shared_ptr<CVK::Camera>(getFreeFlightCamera(image_width, image_height)));

	//Basic OpenGL initialization.	
	auto window = new KIRK::Window("Raytracing", image_width, image_height, 4, 3);
	window->setBackgroundColor(SkyBlue);

	//We use the screenFillShader to draw our texture to a screen filling quad.
	auto screen_shader = CVK::ShaderUtils::ScreenFillShader();
	auto phong_shader = std::make_shared<CVK::ShaderUtils::PhongToScreenShader>(CVK::SceneToCVK::exportScene(scene));

	//////////////////////////////////////////////
	//
	//  CPU Raytracer initialization
	//
	//////////////////////////////////////////////

	//Initialize the cpu_scene with BVH Datastructure
	auto cpu_scene = std::make_shared<KIRK::CPU::Scene>(scene, std::unique_ptr<KIRK::CPU::CPU_DataStructure>(std::make_unique<KIRK::CPU::BVH>()));

	//Initialize the CPU raytracers
	cpu_pathtracer->init(cpu_scene);

	//Set some cpu_pathtracer and bidirectional cpu_pathtracer settings.
	cpu_pathtracer->setSampleCount(20);

	//set default raytracer
	m_active_cpu_raytracer = cpu_pathtracer;


	//////////////////////////////////////////////
	//
	//  GUI initialization
	//
	//////////////////////////////////////////////

	//Initialize ImGUI
	auto gui = std::make_shared<KIRK::Gui>(*window);

	KIRK::GuiScene gui_scene(gui);
	gui_scene.setUpdateCallback([&](unsigned flags)
	{
		cpu_scene->updateFromSceneGraph(true);
		phong_shader->setRenderScene(CVK::SceneToCVK::exportScene(scene));
	});
	gui_scene.buildSceneGui(scene, *camera_wrapper);


	auto loadFile = [&](fs::path path)
	{
		LOG_INFO("Scene importing from %", path);
		scene = KIRK::SceneGraph::makeSceneGraph(path.string());//Memory usage verringert sich beim Aufruf dieser Zeile warum auch immer.

		//cpu scene and raytracer
		cpu_scene->setSceneGraph(scene);
		cpu_scene->updateFromSceneGraph(true);
		cpu_pathtracer->init(cpu_scene);
		//openGL stuff
		phong_shader->setRenderScene(CVK::SceneToCVK::exportScene(scene));
		screen_shader = CVK::ShaderUtils::ScreenFillShader();
		camera_wrapper->setCVKCamera(scene, std::shared_ptr<CVK::Camera>(getFreeFlightCamera(image_width, image_height)));
		gui_scene.setUpdateCallback([&](unsigned flags)
		{
			cpu_scene->updateFromSceneGraph(true);
			phong_shader->setRenderScene(CVK::SceneToCVK::exportScene(scene));
		});
		gui_scene.buildSceneGui(scene, *camera_wrapper);
	};


	//gui raytracer selection
	gui->make<KIRK::GuiWindow>("Render Mode", true)->make<KIRK::GuiNamedLambdaElement>("modes", [&]()
	{
		//------------ Raytracer Selection -----------------

		static int cpu_ds_id = 1;

		//select rendering platform
		ImGui::RadioButton("OpenGL", &m_platform_used, 0); ImGui::SameLine();
		ImGui::RadioButton("CPU", &m_platform_used, 1);

		//select Raytracer usage
		if (m_platform_used == 1)
		{
			//if the raytracer is active change the button color to green
			bool pop_color = false;
			if (m_is_raytracer_active)
			{
				ImGui::PushStyleColor(ImGuiCol_Button, ImColor::HSV(0.3f, 0.6f, 0.6f));
				ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImColor::HSV(0.3f, 0.7f, 0.7f));
				ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImColor::HSV(0.3f, 0.8f, 0.8f));
				pop_color = true;
			}
			if (ImGui::Button("Render Image", ImVec2(250.0f, 32.0f)))
			{
				LOG_INFO("\n       Rendering new Image \n       -------------------");
				m_is_raytracer_active = true;
			}
			if (pop_color) { ImGui::PopStyleColor(3); pop_color = false; }
			ImGui::SameLine();
			ImGui::PushStyleColor(ImGuiCol_Button, ImColor::HSV(0.0f, 0.6f, 0.6f));
			ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImColor::HSV(0.0f, 0.7f, 0.7f));
			ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImColor::HSV(0.0f, 0.8f, 0.8f));
			if (ImGui::Button("Stop Rendering", ImVec2(ImGui::GetContentRegionAvailWidth(), 32.0f)))
			{
				m_is_raytracer_active = false;
				LOG_INFO("Stopped Rendering");
			}
			ImGui::PopStyleColor(3);
			ImGui::Separator();
		}


		if (m_platform_used == 1) // CPU is used
		{
			//select cpu datastructure
			const char* cpu_datastructures[] = { "No Datastructure", "Bounding Volume Hierarchy", "KD Tree" };
			if (ImGui::Combo("Datastructure", &cpu_ds_id, cpu_datastructures, static_cast<int>(sizeof(cpu_datastructures) / sizeof(*cpu_datastructures)))) { setDatastructure(cpu_ds_id); }
			ImGui::Separator();
			ImGui::TextColored(ImVec4(0.f, 1.f, 0.f, 1.f), "Settings");

			//Draw and set up Pathtracer
			cpu_pathtracer->draw();
			setActiveCPURaytracer(cpu_pathtracer);
			m_show_render_texture = true;
		}
		else // OpenGL is used
		{
			m_is_raytracer_active = false;
			m_show_render_texture = false;
		}
	});

	auto menu_bar = gui->make<KIRK::GuiMainMenuBar>();
	menu_bar->addItem("File", "Export", [&]()
	{
		gui->make<KIRK::GuiFileDialog>("Export to...", "Export", "Cancel", SCENES_PATH, false, [&](const KIRK::GuiFileDialog &dialog, const fs::path &path)
		{
			LOG_INFO("Scene exporting to %", path);
			jsonio::writeScene(scene, path.string());
		});
	});
	menu_bar->addItem("File", "Import", [&]()
	{
		gui->make<KIRK::GuiFileDialog>("Import from...", "Import", "Cancel", SCENES_PATH, true, [&](const KIRK::GuiFileDialog &dialog, const fs::path &path)
		{
			loadFile(path);
		});
	});
	menu_bar->addSeparator("File");
	menu_bar->addItem("File", "Save Render", saveRenderTexture);
	menu_bar->addSeparator("File");
	menu_bar->addItem("File", "Exit", [&]()
	{
		window->close();
	});

	menu_bar->addItem("Camera", "Trackball", [&]()
	{
		camera_wrapper->setCVKCamera(scene, std::shared_ptr<CVK::Camera>(getTrackballCamera(image_width, image_height)));
	});

	menu_bar->addItem("Camera", "Free Flight", [&]()
	{
		camera_wrapper->setCVKCamera(scene, std::shared_ptr<CVK::Camera>(getFreeFlightCamera(image_width, image_height)));
	});

	menu_bar->enableViewMenu("View");

	auto fileDropCallback = [&](GLFWwindow *window, int count, const char ** files)
	{
		fs::path path(files[0]);
		loadFile(path);
	};

	//Setting our callbacks
	KIRK::Callbacks::addCharCallback(charCallback);
	KIRK::Callbacks::addFileDropCallback(fileDropCallback);
	glfwSetWindowSizeCallback(window->getGLFWWindow(), resizeCallback);

	//////////////////////////////////////////////
	//
	//  Rendering Loop
	//
	//////////////////////////////////////////////

	while (!window->shouldClose())
	{
		if (!gui->wantsMouseInput() && !m_is_raytracer_active && !m_show_render_texture)
			camera_wrapper->update(window->getGLFWWindow());

		if (m_is_raytracer_active)
		{
			//render the image to the render texture
			if (m_platform_used == 1)
			{
				//Render with the pathtracer to our texture
				m_active_cpu_raytracer->renderToTexture(render_texture_cpu);

				// upload only for cpu tracer to display the data on the gpu
				render_texture_cpu->glUpload();
			}

			//We want to see our image
			m_show_render_texture = true;
		}

		if (m_show_render_texture)
		{
			//If we want to show our texture on screen or after rendering, we just render a screen filling quad.
			if (m_platform_used == 1)//cpu mode is active
			{
				//fill the screen_shader with the cpu render texture
				screen_shader.setRenderTexture(render_texture_cpu->id());
			}
			//render a screen filling quad
			screen_shader.render();
		}
		else
		{
			//If we want to look around in CVK OpenGL mode, we render the scene in OpenGL.
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			phong_shader->render();
		}
		gui->render();
		window->endFrame();
	}

	//Don't leave your garbage around.
	//Clean up behind you please.
	delete window;
	glfwTerminate();
	delete cpu_pathtracer;

	return 0;
}
