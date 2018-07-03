#include "CVK_ConverterUtils.h"
#include "../CVK_2/CVK_Trackball.h"

CVK::CVKCameraSynchronizer::CVKCameraSynchronizer(std::weak_ptr<KIRK::SceneGraph> scene, std::shared_ptr<CVK::Camera> source)
{
	m_kirk_camera = scene.lock()->getActiveCamera().get();
	m_cvk_camera = source;

	int w, h;
	m_cvk_camera->getWidthHeight(&w, &h);
	m_kirk_camera->setResolution(glm::vec2(w, h));

	if (typeid(*(m_cvk_camera.get())) == typeid(CVK::FreeFlight))
	{
		((CVK::FreeFlight*)m_cvk_camera.get())->setUpvector(m_kirk_camera->getLocalUp());
		((CVK::FreeFlight*)m_cvk_camera.get())->setPosition(m_kirk_camera->getLocalPosition());
		((CVK::FreeFlight*)m_cvk_camera.get())->m_direction = glm::normalize(m_kirk_camera->getLocalLookAt());
	}

	else if (typeid(*(m_cvk_camera.get())) == typeid(CVK::Trackball))
	{
		((CVK::Trackball*)m_cvk_camera.get())->setUpvector((m_kirk_camera->getLocalUp()));
		((CVK::Trackball*)m_cvk_camera.get())->setPosition(m_kirk_camera->getLocalPosition());
		((CVK::Trackball*)m_cvk_camera.get())->m_direction = glm::normalize(m_kirk_camera->getLocalLookAt());
	}

}

void CVK::CVKCameraSynchronizer::update(GLFWwindow* window)
{
	//Adapt projection of CVK::Camera to the KIRK::Camera settings.
	m_cvk_camera->getProjection()->setProjMatrix(glm::mat4(glm::perspective(m_kirk_camera->getFieldOfView(), m_kirk_camera->getAspect(), 0.01f, 1000.f)));
	int w, h;
	glfwGetWindowSize(window, &w, &h);
	if (w == 0 || h == 0) {
		LOG_WARNs() << "Tried to resize on a zero-size.";
		return;
	}
	m_cvk_camera->setWidthHeight(w, h);
	m_cvk_camera->update(window);

	//Adapt view matrix
	glm::mat4 view = m_kirk_camera->calculateTransform() * m_cvk_camera->getView();
	m_cvk_camera->setView(view);

	//Apply cvk transformation settings to KIRK::Camera
	//Retrieve CVK::Camera transform
	glm::vec3 position, xvec, yvec, zvec;
	m_cvk_camera->getView(&xvec, &yvec, &zvec, &position);

	//Retrieve camera image size
	int width, height;
	m_cvk_camera->getWidthHeight(&width, &height);

	//Transform the KIRK::Camera to the retrieved coordinate system.
	m_kirk_camera->setResolution(glm::vec2(width, height));
	m_kirk_camera->setTransform(position, -zvec, yvec);
	m_kirk_camera->applyParameters();
}

std::shared_ptr<CVK::Camera> CVK::CVKCameraSynchronizer::getCVKCamera() const
{
	return m_cvk_camera;
}

void CVK::CVKCameraSynchronizer::setCVKCamera(std::weak_ptr<KIRK::SceneGraph> scene, std::shared_ptr<CVK::Camera> source)
{
	m_kirk_camera = scene.lock()->getActiveCamera().get();
	m_cvk_camera = source;

	int w, h;
	m_cvk_camera->getWidthHeight(&w, &h);
	m_kirk_camera->setResolution(glm::vec2(w, h));


	if (typeid(*(m_cvk_camera.get())) == typeid(CVK::FreeFlight))
	{
		((CVK::FreeFlight*)m_cvk_camera.get())->setUpvector(m_kirk_camera->getLocalUp());
		((CVK::FreeFlight*)m_cvk_camera.get())->setPosition(m_kirk_camera->getLocalPosition());
		((CVK::FreeFlight*)m_cvk_camera.get())->m_direction = glm::normalize(m_kirk_camera->getLocalLookAt());
	}
	else if (typeid(*(m_cvk_camera.get())) == typeid(CVK::Trackball))
	{
		((CVK::Trackball*)m_cvk_camera.get())->setUpvector(m_kirk_camera->getLocalUp());
		((CVK::Trackball*)m_cvk_camera.get())->setPosition(m_kirk_camera->getLocalPosition());
		((CVK::Trackball*)m_cvk_camera.get())->setRadius(glm::length(m_kirk_camera->getLocalPosition()));
		//((CVK::Trackball*)m_cvk_camera.get())->m_direction= glm::normalize(m_kirk_camera->getLocalLookAt());
	}
}

namespace CVK
{

	std::shared_ptr<CVK::Node> CVK::SceneToCVK::exportScene(std::shared_ptr<KIRK::SceneGraph> sceneGraph, bool showLightGeometry)
	{
		LOG_DEBUG("Exporting scene to CVK::Node...");

		return toCVKNode(sceneGraph->getRootNode(), showLightGeometry)[0];
	}

	std::vector<std::shared_ptr<CVK::Node>> CVK::SceneToCVK::toCVKNode(const std::shared_ptr<KIRK::SceneNode> sceneNode, bool showLightGeometry)
	{
		std::vector<std::shared_ptr<CVK::Node>> m_nodes;

		//////////////////////////////////////////////////////////////////////////////////////////
		///////		NODE IS LIGHT OR CAMERA
		//////////////////////////////////////////////////////////////////////////////////////////
		if (sceneNode->m_data_type == KIRK::SceneNode::LIGHT || sceneNode->m_data_type == KIRK::SceneNode::CAMERA)
		{
			//Those nodes are leaves, but not necessary for CVK. We should discard them.
			if (sceneNode->m_data_type == KIRK::SceneNode::LIGHT)
			{
				std::shared_ptr<KIRK::Light> light = std::dynamic_pointer_cast<KIRK::Light>(sceneNode->m_data_object);
				CVK::Light *cvk_light;

				if (auto sun = std::dynamic_pointer_cast<KIRK::SunLight>(light))
				{
					KIRK::SunLight t(*sun);
					t.transform(sceneNode->m_data_object->calculateTransform());
					cvk_light = new CVK::Light(glm::vec4(t.m_position, 0), glm::vec3(t.m_color), t.m_direction);
				}
				else if (auto spot = std::dynamic_pointer_cast<KIRK::SpotLight>(light))
				{
					KIRK::SpotLight t(*spot);
					t.transform(sceneNode->m_data_object->calculateTransform());
					cvk_light = new CVK::Light(glm::vec4(t.m_position, 1), glm::vec3(t.m_color), t.m_direction, glm::abs(t.m_outer_angle - t.m_inner_angle), glm::radians(t.m_outer_angle));
				}
				else if (auto quad = std::dynamic_pointer_cast<KIRK::QuadLight>(light))
				{
					KIRK::QuadLight t(*quad);
					t.transform(sceneNode->m_data_object->calculateTransform());
					cvk_light = new CVK::Light(glm::vec4(t.m_position, 1), glm::vec3(t.m_color), t.m_direction, 1, glm::radians(90.f));
				}
				else
				{
					KIRK::PointLight t(*std::dynamic_pointer_cast<KIRK::PointLight>(light));
					cvk_light = new CVK::Light(glm::vec4(t.m_position, 1), glm::vec3(t.m_color), t.m_direction);
				}

				//We also put in our light sources, all as being point lights.
				CVK::State::getInstance()->addLight(cvk_light);

				//Render sphere at the lights location for debugging purpose
				if (showLightGeometry) {
					//transformation for the geometry
					glm::mat4 trans = sceneNode->m_data_object->calculateTransform();
					//Node and Geometry for the Light position
					std::shared_ptr<CVK::Node> light_node = std::make_shared<CVK::Node>("Light_Geometry_Node");
					std::shared_ptr<CVK::Sphere> light_geometry = std::make_shared<CVK::Sphere>(light->m_position, 0.3f);//sphere at the lights position with radius 0.3
					//set geometry and transform for node
					light_node->setGeometry(light_geometry);
					light_node->setModelMatrix(trans);
					m_nodes.push_back(light_node);
					//Node and Geometry for the Light direction
					std::shared_ptr<CVK::Node> light_node_dir = std::make_shared<CVK::Node>("Light_Geometry_Node");
					std::shared_ptr<CVK::Cone> light_geometry_dir = std::make_shared<CVK::Cone>(light->m_position, light->m_direction, 0.05f, 0.0f , 5);//cone for lights direction
					 //set geometry and transform for node
					light_node_dir->setGeometry(light_geometry_dir);
					light_node_dir->setModelMatrix(trans);
					m_nodes.push_back(light_node_dir);
				}
			}

			//Return an vector of nodes (With nodes for light geometry or empty when data object is a camera).
			return m_nodes;
		}

		//////////////////////////////////////////////////////////////////////////////////////////
		///////		NODE IS MESH
		//////////////////////////////////////////////////////////////////////////////////////////
		if (sceneNode->m_data_type == KIRK::SceneNode::MESH)
		{
			//If we have a mesh node, there is no chance, we have children, as mesh nodes are always leaves.
			std::shared_ptr<KIRK::Mesh> kirk_mesh = std::dynamic_pointer_cast<KIRK::Mesh>(sceneNode->m_data_object);

			//Use already separated meshes, because CVK can only use one material per mesh
			std::vector<KIRK::Mesh> separate_meshes = kirk_mesh->separateByMaterial();

			for (int j = 0; j < separate_meshes.size(); j++)
			{
				KIRK::Mesh &mesh = separate_meshes[j];

				//Naming conventions... yeah.
				std::shared_ptr<CVK::Node> mesh_node = std::make_shared<CVK::Node>("Mesh_Node");
				std::shared_ptr<CVK::Geometry> node_mesh = std::make_shared<CVK::Geometry>();

				///////////////////////////////////////////////////////////
				//
				// A CVK::Geometry consists of indices, vertices,
				// normals and uv-coords
				// We can let it calculate its tangents by itself,
				// the other stuff has to be set manually.
				//
				///////////////////////////////////////////////////////////

				//First we put all vertices, normals and uvs
				//For that we need to reserve the needed sizes in the vectors.
				node_mesh->getVertices()->reserve(mesh.m_vertices.size());
				node_mesh->getNormals()->reserve(mesh.m_vertices.size());
				node_mesh->getUVs()->reserve(mesh.m_vertices.size());

				//Then we can add everything.
				for (KIRK::Mesh::vertex vertex : mesh.m_vertices)
				{
					node_mesh->getVertices()->emplace_back(glm::vec4(vertex.position, 1));
					node_mesh->getNormals()->emplace_back(glm::vec3(glm::vec4(vertex.normal, 0)));
					node_mesh->getUVs()->emplace_back(vertex.texcoord);
				}

				//Then we put in all our indices.
				//Reserve again.
				node_mesh->getIndex()->reserve(mesh.m_indices.size());
				for (unsigned int i = 0; i < mesh.m_indices.size(); i++)
				{
					//Add vertices according to loaded indices.
					node_mesh->getIndex()->emplace_back(mesh.m_indices[i]);
				}

				//////////////////////////////////////////////////////////////////////////////////
				//This is where OpenGL is needed. So if it crashes somewhere here, you
				//have to ensure that you have initialized a valid OpenGL context before.
				node_mesh->createBuffers();
				//////////////////////////////////////////////////////////////////////////////////

				//Take the default material.
				node_mesh->setMaterialIndex(0);
				std::shared_ptr<CVK::Material> material = std::make_shared<CVK::Material>(glm::vec3(mesh.m_materials[0]->m_diffuse.value),
					glm::vec3(mesh.m_materials[0]->m_specular.value),
					mesh.m_materials[0]->shininess(),
					mesh.m_materials[0]->m_reflectivity.value);

				// Set Diffuse map if available.
				if (mesh.m_materials[0]->m_diffuse.texture)
				{
					mesh.m_materials[0]->m_diffuse.texture->glUpload();
					material->setTexture(CVK::TextureType::COLOR_TEXTURE, mesh.m_materials[0]->m_diffuse.texture->id());
				}

				//Set Normal map is available
				if (mesh.m_materials[0]->m_normal.texture)
				{
					mesh.m_materials[0]->m_normal.texture->glUpload();
					material->setTexture(CVK::TextureType::NORMAL_TEXTURE, mesh.m_materials[0]->m_normal.texture->id());
				}

				// Finally assign mesh, transformation and material
				mesh_node->setGeometry(node_mesh);
				glm::mat4 transform = sceneNode->m_data_object->calculateTransform();
				mesh_node->setModelMatrix(transform);
				mesh_node->setMaterial(material);

				//Push node with mesh geometry
				m_nodes.push_back(mesh_node);

				//////////////////////////////////////////////////////////////////////////////////
				//
				// Add CVK::Cone Geometries to the nodes for every fur fiber(instanced) in the KIRK::Mesh
				//
				//////////////////////////////////////////////////////////////////////////////////



				//If we have fur fiber information in KIRK::Mesh, we add a Cone Geometry for every fiber to CVK scene.
				if (!mesh.m_furFibers.empty()) {
					//create different colored materials for every cone that exists in a single fur fiber
					std::vector<std::shared_ptr<CVK::Material>> fiber_materials;
					fiber_materials.push_back(std::make_shared<CVK::Material>(glm::vec3(1.f, 0.f, 0.f), glm::vec3(1.0f), 0.01f, 1.5f));// red
					fiber_materials.push_back(std::make_shared<CVK::Material>(glm::vec3(0.f, 1.f, 0.f), glm::vec3(1.0f), 0.01f, 1.5f));// green
					fiber_materials.push_back(std::make_shared<CVK::Material>(glm::vec3(0.f, 0.f, 1.f), glm::vec3(1.0f), 0.01f, 1.5f));// blue
					fiber_materials.push_back(std::make_shared<CVK::Material>(glm::vec3(1.f, 1.f, 0.f), glm::vec3(1.0f), 0.01f, 1.5f));// yellow
					fiber_materials.push_back(std::make_shared<CVK::Material>(glm::vec3(1.f, 0.f, 1.f), glm::vec3(1.0f), 0.01f, 1.5f));// magenta
					fiber_materials.push_back(std::make_shared<CVK::Material>(glm::vec3(0.f, 1.f, 1.f), glm::vec3(1.0f), 0.01f, 1.5f));// cyan
					fiber_materials.push_back(std::make_shared<CVK::Material>(glm::vec3(1.f, 1.f, 1.f), glm::vec3(1.0f), 0.01f, 1.5f));// white
					fiber_materials.push_back(std::make_shared<CVK::Material>(glm::vec3(0.f, 0.f, 0.f), glm::vec3(1.0f), 0.01f, 1.5f));// black
					fiber_materials.push_back(std::make_shared<CVK::Material>(glm::vec3(138.0f / 255.0f, 7.0f / 255.0f, 7.0f / 255.0f), glm::vec3(1.0f), 0.01f, 1.5f));// blood red
					fiber_materials.push_back(std::make_shared<CVK::Material>(glm::vec3(0.7f, 0.9f, 1.f), glm::vec3(1.0f), 0.01f, 1.5f));// sky blue

					//Parent Node for all fiber_ConeGeometry_Nodes
					std::shared_ptr<CVK::Node> fiber_node = std::make_shared<CVK::Node>("Fiber_Node");

					//Iterate once over all fiber positions of one fiber and create a node with geometry for every cone
					for (int j = 0; j < mesh.m_furFibers[0].fiber_positions.size() - 1; j++) {
						//get base and apex position for cone
						glm::vec3 basepos = mesh.m_furFibers[0].fiber_positions[j];
						glm::vec3 apexpos = mesh.m_furFibers[0].fiber_positions[j + 1];
						float baseradius = mesh.m_furFibers[0].fiber_radius[j];
						//move base position a bit to hide cone edges in the fiber struct
						basepos -= 0.008f * (apexpos - basepos);
						//lower base radius a bit so it doesnt stick out of the previous cylinder of the fiber						
						baseradius -= (j > 3) ? 0.1f * baseradius : 0.05f * baseradius;//higher multiplier at the top of the fiber, because the values are smaller there
						//create CVK::Cone geometry with fiber_positions and fiber_radius
						std::shared_ptr<CVK::Cone> fur_geometry = std::make_shared<CVK::Cone>(basepos, apexpos,
							baseradius, mesh.m_furFibers[0].fiber_radius[j + 1], 5);
						//create node for the geometry object
						std::shared_ptr<CVK::Node> cone_node = std::make_shared<CVK::Node>("Fiber_ConeGeometry_Node");
						//Assign Cone geometry to node
						cone_node->setGeometry(fur_geometry);
						cone_node->setMaterial(fiber_materials[(j < fiber_materials.size()) ? j : 0]);
						//Add cone_node to parent node
						fiber_node->addChild(cone_node);
					}
					//Now iterate over every fur fiber in the Mesh
					for (int i = 0; i < mesh.m_furFibers.size(); i++) {
						//Instancing node for fur fibers. They get the node with the geometry as child and only change the transformation matrix for the new fiber position
						std::shared_ptr<CVK::Node> fiber_inst_node = std::make_shared<CVK::Node>("Fiber_Instancing_Node");
						//add the parent node(with geometries stored) as child
						fiber_inst_node->addChild(fiber_node);
						//calculate translation vector. Since our geometry is created at the position of the first fiber in our mesh 
						//we have to use the difference between the first fiber and the current fiber. Also we only need to use the start positions of the fiber.
						glm::vec3 trans_vec = mesh.m_furFibers[i].fiber_positions[0] - mesh.m_furFibers[0].fiber_positions[0];
						//Now calculate the actual transformation matrix with the meshes matrix as base and our translation vector on top of that
						glm::mat4 trans_mat = glm::translate(transform, trans_vec);
						//set transformation matrix for current fiber
						fiber_inst_node->setModelMatrix(trans_mat);
						//Push node to all other scene nodes
						m_nodes.push_back(fiber_inst_node);						
					}
				}
			}
		}

		//If we just have an empty node, we convert all its children and add them to it.
		m_nodes.push_back(std::make_shared<CVK::Node>("Sub_Root"));

		for (std::shared_ptr<KIRK::SceneNode> child : sceneNode->m_children)
		{
			//First create all children.
			std::vector<std::shared_ptr<CVK::Node>> child_nodes = toCVKNode(child, showLightGeometry);

			if (child_nodes.empty())
				continue;

			for (auto &&node : child_nodes)
			{
				m_nodes[0]->addChild(node);
			}
		}

		return m_nodes;
	}
}
