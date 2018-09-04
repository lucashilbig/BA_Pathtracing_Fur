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
#include "KIRK/Utils/Clock.h"

#include "KIRK/CPU/CPU_Raytracer/Simple_CPU_Raytracer.h"
#include "KIRK/CPU/CPU_Raytracer/CPU_PathTracer.h"
#include "KIRK/CPU/CPU_Raytracer/CPU_BidirectionalPathTracer.h"
#include "KIRK/CPU/CPU_Raytracer/CPU_PhotonMapper.h"
#include "KIRK/GPU/GPU_Raytracer/Simple_Raytracer.h"
#include "KIRK/GPU/GPU_Raytracer/Pathtracer.h"

#include "KIRK/CPU/CPU_Datastructures/CPU_BVH.h"
#include "KIRK/GPU/GPU_Datastructures/BVH.h"
#include "KIRK/GPU/GPU_Datastructures/NoDatastructure.h"

#include <externals/ImGui/imgui_impl_glfw_gl3.h>
#include "KIRK/GPU/GPU_Raytracer/BsdfBuilder.h"

//////////////////////////////////////////////
//
//  Member Variables
//
//////////////////////////////////////////////

//booleans used for settings
bool m_is_raytracer_active = false;
bool m_show_render_texture = false;
bool m_use_clock = false;

//platform which is used to render image
static int m_platform_used = 0; //value 0 for OpenGL ; 1 for CPU ; 2 for GPU

//ID of datastructure which is used. Default: 1 BVH
static int m_gpu_ds_id = 1;
static int m_cpu_ds_id = 1;

//create cpu raytracer objects
KIRK::CPU::CPU_Raytracer *m_active_cpu_raytracer;
auto cpu_raytracer = new KIRK::CPU::SimpleCPURaytracer();
auto cpu_pathtracer = new KIRK::CPU::PathTracer();
auto cpu_bipathtracer = new KIRK::CPU::BidirectionalPathTracer();
auto cpu_photonmapper = new KIRK::CPU::CPU_PhotonMapper();

//create gpu raytracer objects
KIRK::GPU::Raytracer *m_active_gpu_raytracer;
KIRK::GPU::Simple_Raytracer *gpu_raytracer;
KIRK::GPU::Pathtracer *gpu_pathtracer;

//textures we use for rendering the images
std::shared_ptr<KIRK::Texture> render_texture_cpu;
std::shared_ptr<KIRK::Texture> render_texture_gpu;

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

		//reset the bi-/pathtracer
		if (active_cpu_raytracer == cpu_pathtracer)
			static_cast<KIRK::CPU::PathTracer*>(m_active_cpu_raytracer)->reset();
		else if (active_cpu_raytracer == cpu_bipathtracer)
			static_cast<KIRK::CPU::BidirectionalPathTracer*>(m_active_cpu_raytracer)->reset();
	}
}

void setActiveGPURaytracer(KIRK::GPU::Raytracer *active_gpu_raytracer)
{
	if (!m_is_raytracer_active && m_active_gpu_raytracer != active_gpu_raytracer)
	{
		//set active gpu raytracer
		m_active_gpu_raytracer = active_gpu_raytracer;
		m_platform_used = 2;
		m_active_gpu_raytracer->updateDataFromScene();

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
		//GPU Raytracer is used
		if (m_platform_used == 2 && m_gpu_ds_id != datastructure_id)
		{
			m_gpu_ds_id = datastructure_id;
			switch (datastructure_id)
			{
			case 0://GPU::NoDatastructure
				LOG_INFO("Switching GPU Datastructure to: NoDatastructure ...");
				m_active_gpu_raytracer->getScene().setDataStructure(std::unique_ptr<KIRK::GPU::Datastructure>(std::make_unique<KIRK::GPU::NoDatastructure>()), true);
				LOG_INFO("... done.");
				break;
			case 1://GPU::BVH which is Bounding Volume Hierarchy
				LOG_INFO("Switching GPU Datastructure to: Bounding Volume Hierarchy ...");
				m_active_gpu_raytracer->getScene().setDataStructure(std::unique_ptr<KIRK::GPU::Datastructure>(std::make_unique<KIRK::GPU::BVH>()), true);
				LOG_INFO("... done.");
				break;
			default://default datastructure is BVH
				LOG_INFO("Switching GPU Datastructure to: Bounding Volume Hierarchy ...");
				m_active_gpu_raytracer->getScene().setDataStructure(std::unique_ptr<KIRK::GPU::Datastructure>(std::make_unique<KIRK::GPU::BVH>()), true);
				LOG_INFO("... done.");
				break;
			}
		}
		//CPU Raytracer is used
		else if (m_platform_used == 1 && m_cpu_ds_id != datastructure_id)
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
			case 3://CPU::Octree
				LOG_INFO("Switching CPU Datastructure to: Octree ...");
				m_active_cpu_raytracer->getScene().setDataStructure(std::unique_ptr<KIRK::CPU::CPU_DataStructure>(std::make_unique<KIRK::CPU::Octree>
					(m_active_cpu_raytracer->getScene().getBounds()[0], m_active_cpu_raytracer->getScene().getBounds()[1], 5)), true);
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
	else if (m_platform_used == 2)//gpu mode is active
	{
		render_texture_gpu->glDownload();
		render_texture_gpu->saveTo(SCREENSHOTS_PATH "/render_out_gpu.png");
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
	render_texture_gpu->resize(width, height);
	if (gpu_raytracer)
		gpu_raytracer->windowResizeCallback();
	if (gpu_pathtracer)
		gpu_pathtracer->windowResizeCallback();
	cpu_pathtracer->windowResizeCallback();
	cpu_bipathtracer->windowResizeCallback();
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
	case '2'://switch to GPU mode
		m_platform_used = 2;
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
	std::string file = std::string(SCENES_PATH) + arg_map.get('r', std::string("/Lucy/scene.json"));
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

	//////////////////////////////////////////////
	//
	//  General stuff
	//
	//////////////////////////////////////////////

	//We need a render texture to render to
	render_texture_cpu = std::make_shared<KIRK::Texture>(image_width, image_height, int(KIRK::Texture::TEX_RGB_UB));
	render_texture_gpu = std::make_shared<KIRK::Texture>(image_width, image_height, int(KIRK::Texture::TEX_RGBA_UB));

	//For the cvk camera on the one side to be up to date on the kirk camera's parameters and for the
	//kirk camera's position to be up to date with the user set cvk camera, we have to use this wrapper to keep all tidy
	camera_wrapper = new CVK::CVKCameraSynchronizer(scene, std::shared_ptr<CVK::Camera>(getFreeFlightCamera(image_width, image_height)));

	//Basic OpenGL initialization.
	LOG_INPUT("Can you use OpenGL 4.3 with ARB_bindless_texture? (y/n):");
	char selection;
	bool inferiorPC = true;
	std::cin >> selection;
	switch (selection)
	{
	case 'y':
	case 'Y':
		inferiorPC = false;
	default:;
	}


	auto window = inferiorPC ? new KIRK::Window("Raytracing", image_width, image_height, 3, 3) : new KIRK::Window("Raytracing", image_width, image_height, 4, 3);
	window->setBackgroundColor(SkyBlue);

	//Check if device got OpenGL version 4.3 or higher (Needed for ComputeShader)
	bool compute_compatible = window->getVersionCode() >= 403;

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
	cpu_raytracer->init(cpu_scene);
	cpu_pathtracer->init(cpu_scene);
	cpu_bipathtracer->init(cpu_scene);
	cpu_photonmapper->init(cpu_scene);

	//Set some cpu_raytracer settings.
	cpu_raytracer->setFlags(KIRK::CPU::SimpleCPURaytracer::RTFLAG_USE_SOFT_SHADOWS |
		KIRK::CPU::SimpleCPURaytracer::RTFLAG_USE_PROGRESSBAR);

	//Set some cpu_pathtracer and bidirectional cpu_pathtracer settings.
	cpu_pathtracer->setSampleCount(20);
	cpu_bipathtracer->setSampleCount(20);

	//set default raytracer
	m_active_cpu_raytracer = cpu_raytracer;

	//////////////////////////////////////////////
	//
	//  GPU Raytracer initialization
	//
	//////////////////////////////////////////////

	std::shared_ptr<KIRK::GPU::Scene> gpu_scene;

	//make shure the device is capable of using compute shader (openGL version >= 4.3)
	if (compute_compatible)
	{
		gpu_raytracer = new KIRK::GPU::Simple_Raytracer();
		gpu_pathtracer = new KIRK::GPU::Pathtracer(false);

		//Initialize the gpu_scene with BVH datastructure
		gpu_scene = std::make_shared<KIRK::GPU::Scene>(scene, std::unique_ptr<KIRK::GPU::Datastructure>(std::make_unique<KIRK::GPU::BVH>()));

		//Initialize the gpu_raytracer and gpu_pathtracer with gpu_scene
		gpu_raytracer->init(gpu_scene);
		gpu_pathtracer->init(gpu_scene);
		//set default raytracer
		m_active_gpu_raytracer = gpu_raytracer;
	}

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
		gpu_scene->updateFromSceneGraph(true);
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
		cpu_raytracer->init(cpu_scene);
		cpu_pathtracer->init(cpu_scene);
		cpu_bipathtracer->init(cpu_scene);
		cpu_photonmapper->init(cpu_scene);
		//gpu scene and raytracer
		if (compute_compatible)
		{
			gpu_scene = std::make_shared<KIRK::GPU::Scene>(scene, std::unique_ptr<KIRK::GPU::Datastructure>(std::make_unique<KIRK::GPU::BVH>()));
			gpu_raytracer->init(gpu_scene);
			gpu_pathtracer->init(gpu_scene);
			m_gpu_ds_id = 1;
		}
		//openGL stuff
		phong_shader->setRenderScene(CVK::SceneToCVK::exportScene(scene));
		screen_shader = CVK::ShaderUtils::ScreenFillShader();
		camera_wrapper->setCVKCamera(scene, std::shared_ptr<CVK::Camera>(getFreeFlightCamera(image_width, image_height)));
		gui_scene.setUpdateCallback([&](unsigned flags)
		{
			gpu_scene->updateFromSceneGraph(true);
			cpu_scene->updateFromSceneGraph(true);
			phong_shader->setRenderScene(CVK::SceneToCVK::exportScene(scene));
		});
		gui_scene.buildSceneGui(scene, *camera_wrapper);
	};


	//gui raytracer selection
	gui->make<KIRK::GuiWindow>("Render Mode", true)->make<KIRK::GuiNamedLambdaElement>("modes", [&]()
	{
		//------------ Raytracer Selection -----------------

		static int gpu_ds_id = 1;
		static int cpu_ds_id = 1;
		static int raytracer = 0;

		//select rendering platform
		ImGui::RadioButton("OpenGL", &m_platform_used, 0); ImGui::SameLine();
		ImGui::RadioButton("CPU", &m_platform_used, 1);
		if (compute_compatible) { ImGui::SameLine(); ImGui::RadioButton("GPU", &m_platform_used, 2); }

		//select Raytracer usage
		if (m_platform_used == 1 || m_platform_used == 2)
		{
			ImGui::RadioButton("Raytracer", &raytracer, 0); ImGui::SameLine();
			ImGui::RadioButton("Pathtracer", &raytracer, 1);
			if (m_platform_used == 1)
			{
				ImGui::SameLine();
				ImGui::RadioButton("Bidirectional PT", &raytracer, 2);
				ImGui::RadioButton("Photonmapper", &raytracer, 3);
			}
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
			ImGui::Checkbox("Calculate Render Time", &m_use_clock);
			ImGui::Separator();
		}

		if (m_platform_used == 2) // GPU is used
		{
			//select gpu datastructure
			const char* gpu_datastructures[] = { "No Datastructure", "Bounding Volume Hierarchy" };
			if (ImGui::Combo("Datastructure", &gpu_ds_id, gpu_datastructures, static_cast<int>(sizeof(gpu_datastructures) / sizeof(*gpu_datastructures)))) { setDatastructure(gpu_ds_id); }
			ImGui::Separator();
			ImGui::TextColored(ImVec4(0.f, 1.f, 0.f, 1.f), "Settings");

			//set active raytracer
			switch (raytracer)
			{
			case 0://simple Raytracer
				gpu_raytracer->draw();
				setActiveGPURaytracer(gpu_raytracer);
				m_show_render_texture = true;
				break;
			case 1:// Pathtracer
				gpu_pathtracer->draw();
				gpu_pathtracer->m_tonemapper.get()->draw();
				setActiveGPURaytracer(gpu_pathtracer);
				m_show_render_texture = true;
				break;
			default: // default is simple Raytracer
				gpu_raytracer->draw();
				setActiveGPURaytracer(gpu_raytracer);
				m_show_render_texture = true;
				break;
			}
		}
		else if (m_platform_used == 1) // CPU is used
		{
			//select gpu datastructure
			const char* cpu_datastructures[] = { "No Datastructure", "Bounding Volume Hierarchy", "KD Tree", "Octree"};
			if (ImGui::Combo("Datastructure", &cpu_ds_id, cpu_datastructures, static_cast<int>(sizeof(cpu_datastructures) / sizeof(*cpu_datastructures)))) { setDatastructure(cpu_ds_id); }
			ImGui::Separator();
			ImGui::TextColored(ImVec4(0.f, 1.f, 0.f, 1.f), "Settings");

			switch (raytracer)
			{
			case 0://simple Raytracer
				cpu_raytracer->draw();
				setActiveCPURaytracer(cpu_raytracer);
				m_show_render_texture = true;
				break;
			case 1:// Pathtracer
				cpu_pathtracer->draw();
				setActiveCPURaytracer(cpu_pathtracer);
				m_show_render_texture = true;
				break;
			case 2:// Bidirectional Pathtracer
				cpu_bipathtracer->draw();
				setActiveCPURaytracer(cpu_bipathtracer);
				m_show_render_texture = true;
				break;
			case 3:// Photonmapper
				   //cpu_photonmapper->draw();
				setActiveCPURaytracer(cpu_photonmapper);
				m_show_render_texture = true;
			default: // default is simple Raytracer
				cpu_raytracer->draw();
				setActiveCPURaytracer(cpu_raytracer);
				m_show_render_texture = true;
				break;
			}
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
			//render the image to the render texture, when GPU is used
			if (m_platform_used == 2)
			{
				if (m_use_clock)
				{
					START_GPU_CLOCK;
					m_active_gpu_raytracer->renderToTexture(render_texture_gpu);
					LOG_INFO("Render on GPU took % milliseconds", END_GPU_CLOCK / 1000000.0);
				}
				else
				{
					m_active_gpu_raytracer->renderToTexture(render_texture_gpu);
				}
				// uploading the texture here would be unnecessary slow because the data is already on the gpu.
				// also since the cpu data contains garbage, uploading it would overwrite the rendered image with garbage
			}
			else if (m_platform_used == 1) //or CPU is used
			{
				if (m_use_clock && (m_active_cpu_raytracer == cpu_raytracer || m_active_cpu_raytracer == cpu_photonmapper))
				{
					KIRK::Clock<> clock;
					m_active_cpu_raytracer->renderToTexture(render_texture_cpu);
					LOG_INFO("Render on CPU took % milliseconds", clock.getElapsedTime<KIRK::Milliseconds>());
				}
				else
				{
					m_active_cpu_raytracer->renderToTexture(render_texture_cpu);
				}
				// upload only for cpu tracer to display the data on the gpu
				render_texture_cpu->glUpload();
			}

			//We don't want to keep rendering.
			if (m_active_cpu_raytracer == cpu_raytracer && m_platform_used == 1)
				m_is_raytracer_active = false;
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
			else if (m_platform_used == 2)//gpu mode is active
			{
				//fill the screen_shader with the gpu render texture
				screen_shader.setRenderTexture(render_texture_gpu->id());
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
	delete cpu_raytracer;
	delete cpu_bipathtracer;
	delete cpu_photonmapper;
	if (gpu_raytracer)
		delete gpu_raytracer;
	if (gpu_pathtracer)
		delete gpu_pathtracer;

	return 0;
}
