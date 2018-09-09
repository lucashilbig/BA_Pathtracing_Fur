#include "CVK/CVK_Utils/CVK_ShaderUtils.h"
#include "CVK/CVK_Utils/CVK_ConverterUtils.h"
#include "CVK/CVK_2/CVK_LineList.h"

#include "KIRK/Common/SceneGraph.h"

#include "KIRK/Utils/Window.h"
#include "KIRK/Utils/ArgParser.h"
#include "KIRK/Utils/jsonout.h"
#include "KIRK/Utils/Gui/Gui.h"
#include "KIRK/Utils/Gui/GuiScene.h"
#include "KIRK/Utils/Gui/GuiFileDialog.h"
#include "KIRK/Utils/Gui/GuiMenu.h"

#include "KIRK/CPU/CPU_Raytracer/CPU_PathTracer.h"
#include "KIRK/CPU/CPU_Raytracer/Simple_CPU_Raytracer.h"

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
auto cpu_raytracer = new KIRK::CPU::SimpleCPURaytracer();
auto cpu_pathtracer = new KIRK::CPU::PathTracer();

//textures we use for rendering the images
std::shared_ptr<KIRK::Texture> render_texture_cpu;

//wrapper for the camera
CVK::CVKCameraSynchronizer *camera_wrapper;

//Intersection Rays for debugging
std::vector<KIRK::Ray> in_rays;
std::vector<KIRK::Ray> normal_rays;
std::vector<KIRK::Ray> out_rays;
int m_ray_index = 0;//Index of the current rays
int m_out_ray_index = 0;//Index of the current output rays since we have 2(TT-Path) or 3(TRT-Path) for one inpute ray

//phong shader
std::shared_ptr<CVK::ShaderUtils::PhongToScreenShader> phong_shader;

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
	case '+'://View next intersection point
		m_ray_index++;
		m_out_ray_index += 3;
		if (in_rays.size() > 0 && in_rays.size() > m_ray_index - 1) {
			//create parent node
			std::shared_ptr<CVK::Node> parent = std::make_shared<CVK::Node>("Parent_Node");;

			//create CVK::LineList for Input/Normal/Output rays
			std::shared_ptr<CVK::LineList> in_lines = std::make_shared<CVK::LineList>();
			std::shared_ptr<CVK::LineList> normal_lines = std::make_shared<CVK::LineList>();
			std::shared_ptr<CVK::LineList> out_lines = std::make_shared<CVK::LineList>();
			in_lines->init_LineList();
			normal_lines->init_LineList();
			out_lines->init_LineList();

			//Input Rays in Green
			in_lines->add_Line(in_rays[m_ray_index].m_origin, in_rays[m_ray_index].followDistance(glm::length(in_rays[m_ray_index].m_origin - normal_rays[m_out_ray_index].m_origin)), glm::vec3(0.f, 1.f, 0.f));
			//Normal Rays in Blue, Light Blue, White
			normal_lines->add_Line(normal_rays[m_out_ray_index].m_origin, normal_rays[m_out_ray_index].followDistance(0.3f), glm::vec3(0.f, 0.f, 1.f));
			normal_lines->add_Line(normal_rays[m_out_ray_index + 1].m_origin, normal_rays[m_out_ray_index + 1].followDistance(0.3f), glm::vec3(0.f, 1.f, 1.f));
			normal_lines->add_Line(normal_rays[m_out_ray_index + 2].m_origin, normal_rays[m_out_ray_index + 2].followDistance(0.3f), glm::vec3(1.f, 1.f, 1.f));
			//Output Rays in RED. First Refraction ray RED, second ray in YELLOW, third ray in ORANGE
			out_lines->add_Line(out_rays[m_out_ray_index].m_origin, out_rays[m_out_ray_index].followDistance(glm::length(out_rays[m_out_ray_index].m_origin - out_rays[m_out_ray_index + 1].m_origin)), glm::vec3(1.f, 0.f, 0.f));
			out_lines->add_Line(out_rays[m_out_ray_index + 1].m_origin, out_rays[m_out_ray_index + 1].followDistance(glm::length(out_rays[m_out_ray_index+1].m_origin - out_rays[m_out_ray_index + 2].m_origin)), glm::vec3(1.f, 1.f, 0.f));
			out_lines->add_Line(out_rays[m_out_ray_index + 2].m_origin, out_rays[m_out_ray_index + 2].followDistance(2), glm::vec3(1.f, 0.65f, 0.f));

			//finish LineLists
			in_lines->finish_LineList();
			normal_lines->finish_LineList();
			out_lines->finish_LineList();
			//add LineLists to Nodes
			std::shared_ptr<CVK::Node> in_node = std::make_shared<CVK::Node>("Input_Rays_Node");
			in_node->setGeometry(in_lines);
			in_node->setModelMatrix(glm::mat4(1.f));
			std::shared_ptr<CVK::Node> normal_node = std::make_shared<CVK::Node>("Normal_Rays_Node");
			normal_node->setGeometry(normal_lines);
			normal_node->setModelMatrix(glm::mat4(1.f));
			std::shared_ptr<CVK::Node> out_node = std::make_shared<CVK::Node>("Ouput_Rays_Node");
			out_node->setGeometry(out_lines);
			out_node->setModelMatrix(glm::mat4(1.f));

			//add nodes to parent
			parent->addChild(in_node);
			parent->addChild(normal_node);
			parent->addChild(out_node);
			
			//assign new scene to phong shader
			phong_shader->setRenderScene(parent);
		}
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
	std::string file = std::string(SCENES_PATH) + arg_map.get('r', std::string("/Fur_SmallSkinPatch/scene.json"));
	std::string relative_path = arg_map.get('s', std::string(""));
	unsigned int image_width = arg_map.get('w', 1280);
	unsigned int image_height = arg_map.get('h', 720);

	if (relative_path != "")
	{
		file = relative_path;
	}

	//Load up empty default scene. Contains Camera, Lights and Environment
	std::shared_ptr<KIRK::SceneGraph> scene(KIRK::SceneGraph::makeSceneGraph(file));


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
	phong_shader = std::make_shared<CVK::ShaderUtils::PhongToScreenShader>(CVK::SceneToCVK::exportScene(scene, true));

	//////////////////////////////////////////////
	//
	//  CPU Raytracer initialization
	//
	//////////////////////////////////////////////

	//Initialize the cpu_scene with BVH Datastructure
	auto cpu_scene = std::make_shared<KIRK::CPU::Scene>(scene, std::unique_ptr<KIRK::CPU::CPU_DataStructure>(std::make_unique<KIRK::CPU::BVH>()), true);
	//clear scene
	cpu_scene->m_materials.clear();
	cpu_scene->m_scene_objects.clear();

	//create object as cylinder or triangle
	bool useCylinder = true;
	KIRK::Cylinder *obj;

	if (useCylinder) {
		//create cylinder object at coordinate origin radius -= (radius / ((float)i + 5));
		float radius = 0.004f;
		obj = new KIRK::Cylinder(glm::vec3(0.f), glm::vec3(0.f, 1.f, 0.f), radius, radius - (radius / 6.f), &glm::mat4(1.f));
		//Material for fur fibers with marschnerHairBSDF and -Shader
		std::shared_ptr<KIRK::Material> mat = std::make_shared<KIRK::Material>("Cylinder_Mat");
		mat->m_diffuse.value = KIRK::Color::RGBA(0.545f, 0.353f, 0.169f, 1.0f);//Brown color 0.545f, 0.353f, 0.169f
		mat->m_ior = 1.55f;//suggested value from marschner hair paper
		mat->m_current_bsdf = 10;//Intersection Test shade function
		cpu_scene->m_materials.push_back(mat);
		obj->setMaterial(mat.get());
		//add object to scene
		cpu_scene->m_scene_objects.push_back(obj);
	}
	else {//use triangles
		//create fur fiber at position 0,0,0 with height 1
		KIRK::Mesh::furFiber fiber;
		fiber.fiber_positions.push_back(glm::vec3(0));
		fiber.fiber_positions.push_back(glm::vec3(0, 1, 0));
		fiber.fiber_radius.push_back(0.2f);
		fiber.fiber_radius.push_back(0.2f);
		//create triangles from the current fur fiber
		std::vector<KIRK::Triangle *> triangles = std::move(cpu_scene->fiberToTriangles(fiber, glm::mat4(1.f), 5));
		//Material for fur fibers with marschnerHairBSDF and -Shader
		std::shared_ptr<KIRK::Material> mat = std::make_shared<KIRK::Material>("Cylinder_Mat");
		mat->m_diffuse.value = KIRK::Color::RGBA(0.545f, 0.353f, 0.169f, 1.0f);//Brown color 0.545f, 0.353f, 0.169f
		mat->m_ior = 1.55f;//suggested value from marschner hair paper
		mat->m_current_bsdf = 10;//Intersection Test shade function
		cpu_scene->m_materials.clear();
		cpu_scene->m_materials.push_back(mat);
		//push triangles to m_scene_objects
		for (int t = 0; t < triangles.size(); t++)
		{
			triangles[t]->setMaterial(mat.get());
			cpu_scene->m_scene_objects.push_back(triangles[t]);
		}
	}



	cpu_scene->setDataStructure(std::unique_ptr<KIRK::CPU::CPU_DataStructure>(std::make_unique<KIRK::CPU::BVH>()));
	cpu_scene->updateFromSceneGraph(true);


	//Initialize the CPU raytracers
	cpu_raytracer->init(cpu_scene);
	cpu_pathtracer->init(cpu_scene);

	//Set some cpu_raytracer settings.
	cpu_raytracer->setFlags(KIRK::CPU::SimpleCPURaytracer::RTFLAG_USE_SOFT_SHADOWS);

	//Set some cpu_pathtracer and bidirectional cpu_pathtracer settings.
	cpu_pathtracer->setSampleCount(100);
	cpu_pathtracer->setDepth(5);

	//set default raytracer
	m_active_cpu_raytracer = cpu_raytracer;

	//////////////////////////////////////////////
	//
	//  GUI initialization
	//
	//////////////////////////////////////////////

	//Initialize ImGUI
	auto gui = std::make_shared<KIRK::Gui>(*window);

	auto menu_bar = gui->make<KIRK::GuiMainMenuBar>();

	menu_bar->addItem("Camera", "Trackball", [&]()
	{
		camera_wrapper->setCVKCamera(scene, std::shared_ptr<CVK::Camera>(getTrackballCamera(image_width, image_height)));
	});

	menu_bar->addItem("Camera", "Free Flight", [&]()
	{
		camera_wrapper->setCVKCamera(scene, std::shared_ptr<CVK::Camera>(getFreeFlightCamera(image_width, image_height)));
	});


	//Setting our callbacks
	KIRK::Callbacks::addCharCallback(charCallback);
	glfwSetWindowSizeCallback(window->getGLFWWindow(), resizeCallback);

	//////////////////////////////////////////////
	//
	//  Create CVK Scene with intersection Rays
	//
	//////////////////////////////////////////////

	//Render with the raytracer to our texture. This will store intersection rays to m_intersec_rays in the cpu_raytracer
	cpu_raytracer->renderToTexture(render_texture_cpu);

	//get rays from raytracer
	in_rays = cpu_raytracer->getInRays();
	normal_rays = cpu_raytracer->getNormalRays();
	out_rays = cpu_raytracer->getOutRays();

	//Display rays with openGL
	if (in_rays.size() > 0) {
		//create parent node
		std::shared_ptr<CVK::Node> parent = std::make_shared<CVK::Node>("Parent_Node");;

		//create CVK::LineList for Input/Normal/Output rays
		std::shared_ptr<CVK::LineList> in_lines = std::make_shared<CVK::LineList>();
		std::shared_ptr<CVK::LineList> normal_lines = std::make_shared<CVK::LineList>();
		std::shared_ptr<CVK::LineList> out_lines = std::make_shared<CVK::LineList>();
		in_lines->init_LineList();
		normal_lines->init_LineList();
		out_lines->init_LineList();

		//Input Rays in GREEN
		in_lines->add_Line(in_rays[m_ray_index].m_origin, in_rays[m_ray_index].followDistance(glm::length(in_rays[m_ray_index].m_origin - normal_rays[m_out_ray_index].m_origin)), glm::vec3(0.f, 1.f, 0.f));
		//Normal Rays in BLUE
		for(int i = 0; i < normal_rays.size(); i+=2)
			normal_lines->add_Line(normal_rays[i].m_origin, normal_rays[i].followDistance(0.3f), glm::vec3(0.f, 0.f, 1.f));

		//Normal Rays in Blue, Light Blue, White
		//normal_lines->add_Line(normal_rays[m_out_ray_index].m_origin, normal_rays[m_out_ray_index].followDistance(0.3f), glm::vec3(0.f, 0.f, 1.f));
		normal_lines->add_Line(normal_rays[m_out_ray_index + 1].m_origin, normal_rays[m_out_ray_index + 1].followDistance(0.3f), glm::vec3(0.f, 1.f, 1.f));
		normal_lines->add_Line(normal_rays[m_out_ray_index + 2].m_origin, normal_rays[m_out_ray_index + 2].followDistance(0.3f), glm::vec3(1.f, 1.f, 1.f));
		//Output Rays in RED. First Refraction ray RED, second ray in YELLOW, third ray in ORANGE
		out_lines->add_Line(out_rays[m_out_ray_index].m_origin, out_rays[m_out_ray_index].followDistance(glm::length(out_rays[m_out_ray_index].m_origin - out_rays[m_out_ray_index + 1].m_origin)), glm::vec3(1.f, 0.f, 0.f));
		out_lines->add_Line(out_rays[m_out_ray_index + 1].m_origin, out_rays[m_out_ray_index + 1].followDistance(glm::length(out_rays[m_out_ray_index + 1].m_origin - out_rays[m_out_ray_index + 2].m_origin)), glm::vec3(1.f, 1.f, 0.f));
		out_lines->add_Line(out_rays[m_out_ray_index + 2].m_origin, out_rays[m_out_ray_index + 2].followDistance(2), glm::vec3(1.f, 0.65f, 0.f));
		
		//finish LineLists
		in_lines->finish_LineList();
		normal_lines->finish_LineList();
		out_lines->finish_LineList();
		//add LineLists to Nodes
		std::shared_ptr<CVK::Node> in_node = std::make_shared<CVK::Node>("Input_Rays_Node");
		in_node->setGeometry(in_lines);
		in_node->setModelMatrix(glm::mat4(1.f));
		std::shared_ptr<CVK::Node> normal_node = std::make_shared<CVK::Node>("Normal_Rays_Node");
		normal_node->setGeometry(normal_lines);
		normal_node->setModelMatrix(glm::mat4(1.f));
		std::shared_ptr<CVK::Node> out_node = std::make_shared<CVK::Node>("Ouput_Rays_Node");
		out_node->setGeometry(out_lines);
		out_node->setModelMatrix(glm::mat4(1.f));

		//add nodes to parent
		parent->addChild(in_node);
		parent->addChild(normal_node);
		parent->addChild(out_node);
		
		//assign new scene to phong shader
		phong_shader->setRenderScene(parent);
	}

	//Display Cylinder Axis
	/*{
		//create parent node
		std::shared_ptr<CVK::Node> parent = std::make_shared<CVK::Node>("Parent_Node");;

		//create CVK::LineList for Input/Normal/Output rays
		std::shared_ptr<CVK::LineList> axis_lines = std::make_shared<CVK::LineList>();
		axis_lines->init_LineList();

		//Add Axis lines. U-axis RED, V-axis GREEN, W-axis BLUE
		axis_lines->add_Line(obj->getU(), obj->getU() * -2, glm::vec3(1.f, 0.f, 0.f));
		axis_lines->add_Line(obj->getV(), obj->getV() * -2, glm::vec3(0.f, 1.f, 0.f));
		axis_lines->add_Line(obj->getW(), obj->getW() * -2, glm::vec3(0.f, 0.f, 1.f));	
		
		//finish LineLists
		axis_lines->finish_LineList();
		//add LineLists to Nodes
		std::shared_ptr<CVK::Node> in_node = std::make_shared<CVK::Node>("Input_Rays_Node");
		in_node->setGeometry(axis_lines);
		in_node->setModelMatrix(glm::mat4(1.f));

		//add nodes to parent
		parent->addChild(in_node);

		//assign new scene to phong shader
		phong_shader->setRenderScene(parent);
	}*/

	//////////////////////////////////////////////
	//
	//  Rendering Loop
	//
	//////////////////////////////////////////////

	while (!window->shouldClose())
	{
		if (!gui->wantsMouseInput() && !m_is_raytracer_active && !m_show_render_texture)
			camera_wrapper->update(window->getGLFWWindow());

		//Render scene with OpenGL
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		phong_shader->render();

		gui->render();
		window->endFrame();
	}

	//Don't leave your garbage around.
	//Clean up behind you please.
	delete window;
	glfwTerminate();
	delete cpu_pathtracer;
	delete cpu_raytracer;

	return 0;
}
