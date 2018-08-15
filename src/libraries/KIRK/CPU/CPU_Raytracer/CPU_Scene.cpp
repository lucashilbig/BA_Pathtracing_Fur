
#include "CPU_Scene.h"
#include <stdexcept>

KIRK::CPU::Scene::Scene(std::weak_ptr<KIRK::SceneGraph> sceneGraph, std::unique_ptr<CPU_DataStructure> datastructure)
{
	setDataStructure(std::move(datastructure), false);
	setSceneGraph(sceneGraph);
}

KIRK::CPU::Scene::~Scene()
{
	for (unsigned int i = 0; i < m_scene_objects.size(); i++)
		delete m_scene_objects[i];
	m_scene_objects.clear();
}

void KIRK::CPU::Scene::setSceneGraph(std::weak_ptr<KIRK::SceneGraph> sceneGraph)
{
	if (!sceneGraph.expired())
	{
		m_sceneGraph = sceneGraph;
	}

	if (m_sceneGraph.expired())
	{
		throw std::invalid_argument("Empty SceneGraph");
	}
	m_camera = std::unique_ptr<KIRK::Camera>(new KIRK::Camera(*(sceneGraph.lock()->getActiveCamera())));
	m_environment = std::unique_ptr<KIRK::Environment>(new KIRK::Environment(*sceneGraph.lock()->getEnvironment()));

	flattenNode(m_sceneGraph.lock()->getRootNode(), glm::mat4(1.0f));

	computeBounds();
	buildDatastructure();
}

void KIRK::CPU::Scene::setDataStructure(std::unique_ptr<CPU_DataStructure> datastructure, bool buildNew)
{
	if (!datastructure)
		throw std::invalid_argument("Don't give me a null pointer for the datastructure!");
	m_datastructure = std::move(datastructure);
	if (buildNew)
		buildDatastructure();
}

void KIRK::CPU::Scene::updateFromSceneGraph(bool updateAll)
{
	if (m_sceneGraph.expired())
	{
		throw std::invalid_argument("Empty SceneGraph");
	}

	m_camera = std::unique_ptr<KIRK::Camera>(new KIRK::Camera(*(m_sceneGraph.lock()->getActiveCamera().get())));
	m_environment = std::unique_ptr<KIRK::Environment>(new KIRK::Environment(*m_sceneGraph.lock()->getEnvironment()));

	if (updateAll)
	{
		flattenNode(m_sceneGraph.lock()->getRootNode(), glm::mat4(1.0f));

		computeBounds();
		buildDatastructure();
	}
}

void KIRK::CPU::Scene::flattenNode(std::shared_ptr<KIRK::SceneNode> sceneNode, glm::mat4 base_transform)
{
	//for every child of the RootNode "sceneNode" do ...
	for (std::shared_ptr<KIRK::SceneNode> child : sceneNode->m_children)
	{	//If the child is a MeshTypeObject do ...
		if (child->m_data_type == KIRK::SceneNode::MESH)
		{	//Copy the meshes dataObject

			KIRK::Mesh *mesh = (KIRK::Mesh*)child->m_data_object.get();
			//Go over every face
			for (int f = 0; f < mesh->m_faces.size(); f++)
			{
				KIRK::Mesh::face face = mesh->m_faces[f];
				//Converte the meshstructure to the CV::TriangleStructure with transformation
				Triangle *tri = new Triangle(
					glm::vec3(base_transform * child->m_transform * glm::vec4(mesh->m_vertices[face.vertex_index_a].position, 1)),
					glm::vec3(base_transform * child->m_transform * glm::vec4(mesh->m_vertices[face.vertex_index_b].position, 1)),
					glm::vec3(base_transform * child->m_transform * glm::vec4(mesh->m_vertices[face.vertex_index_c].position, 1)),
					glm::normalize(glm::vec3(glm::transpose(glm::inverse(base_transform * child->m_transform)) * glm::vec4(mesh->m_vertices[face.vertex_index_a].normal, 0))),
					glm::normalize(glm::vec3(glm::transpose(glm::inverse(base_transform * child->m_transform)) * glm::vec4(mesh->m_vertices[face.vertex_index_b].normal, 0))),
					glm::normalize(glm::vec3(glm::transpose(glm::inverse(base_transform * child->m_transform)) * glm::vec4(mesh->m_vertices[face.vertex_index_c].normal, 0))),
					mesh->m_vertices[face.vertex_index_a].texcoord,
					mesh->m_vertices[face.vertex_index_b].texcoord,
					mesh->m_vertices[face.vertex_index_c].texcoord);

				m_materials.push_back(mesh->m_materials[face.material_index]);
				tri->setMaterial(mesh->m_materials[face.material_index].get());

				m_scene_objects.push_back(tri);
			}						

			//If we have fur fibers in our mesh, we convert them to cylinders or triangles
			if (!mesh->m_furFibers.empty())
			{
				if (m_fiberAsCylinder) {
					//////
					// FIBER TO CYLINDER
					//////

					//random value generator
					std::random_device rd;
					std::mt19937 mt(rd());
					std::uniform_real_distribution<float> dist(5.0f, std::nextafter(10.0f, DBL_MAX));//std::nextafter so we get the [5,10] interval instead of [5,10)

					//transformation of fur fiber
					glm::mat4 transform = base_transform * child->m_transform;
					//Material for fur fibers
					std::shared_ptr<KIRK::Material> mat = std::make_shared<KIRK::Material>("Fiber_Mat");
					mat->m_diffuse.value = KIRK::Color::RGBA(0.545f, 0.353f, 0.169f, 1.0f);//Brown color
					mat->m_transparency.value = 0.4f;
					mat->m_reflectivity.value = 0.4f;
					mat->m_ior = 1.55f;//suggested value from marschner hair paper
					mat->m_alpha_shift.value = glm::radians(-1.0f * dist(mt));//suggested value from marschner hair paper between -10 and -5 degrees
					mat->m_beta_width.value = glm::radians(dist(mt));//suggested value from marschner hair paper between 5 and 10 degrees
					mat->m_bsdf = std::make_shared<BSDF>("MarschnerHairBSDF", MarschnerHairBSDF::localSample, MarschnerHairBSDF::evaluateLight);//Use marschner hair bsdf in case we use cylinder objects
					mat->m_current_bsdf = 6;
					mat->m_shader =  std::make_shared<Shader>("MarschnerHairShader");//Use marschner hair shader together with bsdf
					m_materials.push_back(mat);

					//Iterate over every fur fiber
					for (int i = 0; i < mesh->m_furFibers.size(); i++)
					{
						KIRK::Mesh::furFiber *fiber = &mesh->m_furFibers[i];
						//Iterate over every cone in the fur fiber
						for (int c = 0; c < fiber->fiber_positions.size() - 1; c++)
						{
							//get base and apex position for cylinder
							glm::vec3 basepos = fiber->fiber_positions[c];
							glm::vec3 apexpos = fiber->fiber_positions[c + 1];
							float baseradius = fiber->fiber_radius[c];
							//move base position a bit to hide cone edges in the fiber struct
							basepos -= 0.008f * (apexpos - basepos);
							//lower base radius a bit so it doesnt stick out of the previous cylinder of the fiber						
							baseradius -= (c > 3) ? 0.1f * baseradius : 0.05f * baseradius;//higher multiplier at the top of the fiber, because the values are smaller there
							//create raytracing cylinder object
							KIRK::Cylinder *obj = new KIRK::Cylinder(basepos, apexpos,
								baseradius, fiber->fiber_radius[c + 1], &transform);
							//set cylinder material
							obj->setMaterial(mat.get());
							//add object to scene
							m_scene_objects.push_back(obj);
						}
					}
				}
				else {
					//////
					// FIBER TO TRIANGLE
					//////
					//Iterate over every fur fiber
					for (int i = 0; i < mesh->m_furFibers.size(); i++)
					{
						//create triangles from the current fur fiber
						std::vector<Triangle *> triangles = std::move(fiberToTriangles(mesh->m_furFibers[i], base_transform * child->m_transform, 5));
						//push triangles to m_scene_objects
						for (int t = 0; t < triangles.size(); t++)
						{
							m_scene_objects.push_back(triangles[t]);
						}
					}
				}
			}
		}
		else if (child->m_data_type == KIRK::SceneNode::LIGHT)
		{
			if (typeid(*(child->m_data_object.get())) == typeid(KIRK::PointLight))
			{
				m_lights.push_back(std::unique_ptr<KIRK::Light>
					(new KIRK::PointLight(*((KIRK::PointLight*)child->m_data_object.get()))));
				m_lights.back()->transform(child->m_data_object->calculateTransform());
			}
			else if (typeid(*(child->m_data_object.get())) == typeid(KIRK::QuadLight))
			{
				m_lights.push_back(std::unique_ptr<KIRK::Light>
					(new KIRK::QuadLight(*((KIRK::QuadLight*)child->m_data_object.get()))));
				m_lights.back()->transform(child->m_data_object->calculateTransform());
			}
			else if (typeid(*(child->m_data_object.get())) == typeid(KIRK::SpotLight))
			{
				m_lights.push_back(std::unique_ptr<KIRK::Light>
					(new KIRK::SpotLight(*((KIRK::SpotLight*)child->m_data_object.get()))));
				m_lights.back()->transform(child->m_data_object->calculateTransform());
			}
			else if (typeid(*(child->m_data_object.get())) == typeid(KIRK::SunLight))
			{
				m_lights.push_back(std::unique_ptr<KIRK::Light>
					(new KIRK::SunLight(*((KIRK::SunLight*)child->m_data_object.get()))));
				m_lights.back()->transform(child->m_data_object->calculateTransform());
			}
		}

		if (child->m_children.size() > 0)
			flattenNode(child, base_transform * child->m_transform);
	}
}

void KIRK::CPU::Scene::buildDatastructure()
{
	if (m_datastructure) m_datastructure->addBaseDataStructure(this);
}

const std::vector<KIRK::Object *> &KIRK::CPU::Scene::getSceneObjects() const
{
	return m_scene_objects;
}

const glm::vec3 *KIRK::CPU::Scene::getBounds() const
{
	return m_bound;
}

void KIRK::CPU::Scene::computeBounds()
{
	glm::vec3 *obj_bound;
	m_bound[0] = glm::vec3(FLT_MAX);
	m_bound[1] = glm::vec3(-FLT_MAX);
	for (unsigned int index = 0; index < m_scene_objects.size(); index++)
	{
		obj_bound = m_scene_objects[index]->getBounds();
		for (unsigned int coord = 0; coord <= 2; coord++)
		{
			if (obj_bound[0][coord] < m_bound[0][coord])
				m_bound[0][coord] = obj_bound[0][coord];
			if (obj_bound[1][coord] > m_bound[1][coord])
				m_bound[1][coord] = obj_bound[1][coord];
		}
	}
}

std::vector<KIRK::Triangle *> KIRK::CPU::Scene::fiberToTriangles(KIRK::Mesh::furFiber fiber, glm::mat4 mesh_transform, unsigned int resolution)
{
	// iniatialize the variable we are going to use
	std::vector<Triangle *> triangles;
	std::vector<glm::vec4> m_vertices, m_normals;
	std::vector<glm::vec2> m_uvs;

	float u, v;
	float radius, phi;
	glm::vec3 q;
	int i, j;
	glm::vec3 n1, n2, n;
	glm::vec3 m_v, m_u, m_w;
	float m_height, m_slope;
	//Material for hair fibers
	std::shared_ptr<KIRK::Material> mat = std::make_shared<KIRK::Material>("Fiber_Mat");
	mat->m_diffuse.value = KIRK::Color::RGBA(0.545f, 0.353f, 0.169f, 1.0f);//Brown color
	mat->m_transparency.value = 0.4f;
	mat->m_reflectivity.value = 0.4f;
	m_materials.push_back(mat);

	//Iterate over every cone in the fur fiber
	for (int c = 0; c < fiber.fiber_positions.size() - 1; c++)
	{
		//initialise cone specific variables
		glm::vec3 m_basepoint = fiber.fiber_positions[c];
		glm::vec3 m_apexpoint = fiber.fiber_positions[c + 1];
		int offset = 0;
		m_v = m_apexpoint - m_basepoint;
		m_height = glm::length(m_v);
		m_v = glm::normalize(m_v);
		//clear vectors so the indexes for triangle creation are correct
		m_vertices.clear();
		m_normals.clear();
		m_uvs.clear();

		/* find two axes which are at right angles to cone_v */
		glm::vec3 tmp(0.f, 1.f, 0.f);
		if (1.f - fabs(glm::dot(tmp, m_v)) < KIRK::cRayEpsilon)
			tmp = glm::vec3(0.f, 0.f, 1.f);

		m_u = glm::normalize(glm::cross(m_v, tmp));
		m_w = glm::normalize(glm::cross(m_u, m_v));

		m_slope = (fiber.fiber_radius[c] - fiber.fiber_radius[c + 1]) / m_height;

		// Envelope
		for (j = 0; j <= resolution; j++)  //radius
			for (i = 0; i <= resolution; i++) //phi
			{
				u = i / (float)resolution;
				phi = 2 * glm::pi<float>() * u;
				v = j / (float)resolution;
				v = m_height * v;

				radius = fiber.fiber_radius[c] - m_slope * v;
				q = m_basepoint + radius * sinf(phi) * m_u + v * m_v + radius * cosf(phi) * m_w;

				float t = glm::dot(q, m_v) - glm::dot(m_basepoint, m_v);
				glm::vec3 q_1 = q - t * m_v;
				glm::vec3 n = glm::normalize(q_1 - m_basepoint);
				n = glm::normalize(n + m_slope * m_v);

				m_vertices.push_back(glm::vec4(q, 1.0f));
				m_normals.push_back(glm::vec4(n, 0.0f));
				m_uvs.push_back(glm::vec2(u, v / m_height));
			}


		// create triangles from above calculated vertices 
		for (j = 0; j < resolution; j++)
		{
			for (i = 0; i < resolution; i++)
			{
				// 1. Triangle
				Triangle *tri1 = new Triangle(
					glm::vec3(mesh_transform * m_vertices[offset + i + resolution + 1]),
					glm::vec3(mesh_transform * m_vertices[offset + i]),
					glm::vec3(mesh_transform * m_vertices[offset + i + 1]),
					glm::normalize(glm::vec3(glm::transpose(glm::inverse(mesh_transform)) * m_normals[offset + i + resolution + 1])),
					glm::normalize(glm::vec3(glm::transpose(glm::inverse(mesh_transform)) * m_normals[offset + i])),
					glm::normalize(glm::vec3(glm::transpose(glm::inverse(mesh_transform)) * m_normals[offset + i + 1])),
					m_uvs[offset + i + resolution + 1],
					m_uvs[offset + i],
					m_uvs[offset + i + 1]);
				tri1->setMaterial(mat.get());

				// 2. Triangle
				Triangle *tri2 = new Triangle(
					glm::vec3(mesh_transform * m_vertices[offset + i + 1]),
					glm::vec3(mesh_transform * m_vertices[offset + i + resolution + 1 + 1]),
					glm::vec3(mesh_transform * m_vertices[offset + i + resolution + 1]),
					glm::normalize(glm::vec3(glm::transpose(glm::inverse(mesh_transform)) * m_normals[offset + i + 1])),
					glm::normalize(glm::vec3(glm::transpose(glm::inverse(mesh_transform)) * m_normals[offset + i + resolution + 1 + 1])),
					glm::normalize(glm::vec3(glm::transpose(glm::inverse(mesh_transform)) * m_normals[offset + i + resolution + 1])),
					m_uvs[offset + i + 1],
					m_uvs[offset + i + resolution + 1 + 1],
					m_uvs[offset + i + resolution + 1]);
				tri2->setMaterial(mat.get());

				//Add triangles to return vector
				triangles.push_back(tri1);
				triangles.push_back(tri2);
			}
			offset += resolution + 1;
		}
	}
	return triangles;
}


