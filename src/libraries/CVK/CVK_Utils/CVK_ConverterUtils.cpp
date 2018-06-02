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
		((CVK::FreeFlight*)m_cvk_camera.get())->m_direction= glm::normalize(m_kirk_camera->getLocalLookAt());
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

	std::shared_ptr<CVK::Node> CVK::SceneToCVK::exportScene(std::shared_ptr<KIRK::SceneGraph> sceneGraph)
	{
		LOG_DEBUG("Exporting scene to CVK::Node...");

		return toCVKNode(sceneGraph->getRootNode())[0];
	}

	std::vector<std::shared_ptr<CVK::Node>> CVK::SceneToCVK::toCVKNode(const std::shared_ptr<KIRK::SceneNode> sceneNode)
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
			}

			//Return an empty vector of nodes.
			return std::vector<std::shared_ptr<CVK::Node>>();
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

				//If we have fur fiber information in KIRK::Mesh, we add a LineList Geometry to CVK scene. This contains every fiber in the mesh as a line
				if (!mesh.m_furFaces.empty()) {
					//Node and geometry as LineList
					std::shared_ptr<CVK::Node> fur_node = std::make_shared<CVK::Node>("LineList_Node");
					std::shared_ptr<CVK::LineList> fur_geometry = std::make_shared<CVK::LineList>();
					//Init the LineList once
					fur_geometry->init_LineList();

					//Iterate over every face in the mesh
					for (int i = 0; i < mesh.m_furFaces.size(); i++) {
						//Iterate over the fiber positions in the face
						for (int j = 0; j < mesh.m_furFaces[i].fiber_positions.size() - 1; j++) {
							//Add the Lines of the current fur fiber to the LineList
							fur_geometry->add_Line(mesh.m_furFaces[i].fiber_positions[j], mesh.m_furFaces[i].fiber_positions[j+1], glm::vec3(255.0f, 0.0f, 0.0f));
						}
					}

					//Finish the LineList
					fur_geometry->finish_LineList();
					//Assign LineList to node
					fur_node->setGeometry(fur_geometry);
					fur_node->setModelMatrix(transform);
					//Push node with LineList to m_nodes
					m_nodes.push_back(fur_node);
				}
			}
		}

		//If we just have an empty node, we convert all its children and add them to it.
		m_nodes.push_back(std::make_shared<CVK::Node>("Sub_Root"));

		for (std::shared_ptr<KIRK::SceneNode> child : sceneNode->m_children)
		{
			//First create all children.
			std::vector<std::shared_ptr<CVK::Node>> child_nodes = toCVKNode(child);

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
