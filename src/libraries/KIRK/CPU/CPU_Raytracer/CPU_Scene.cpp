
#include "CPU_Scene.h"
#include <stdexcept>

KIRK::CPU::Scene::Scene(std::weak_ptr<KIRK::SceneGraph> sceneGraph, std::unique_ptr<CPU_DataStructure> datastructure)
{
    setDataStructure(std::move(datastructure), false);
	setSceneGraph(sceneGraph);
}

KIRK::CPU::Scene::~Scene()
{
    for(unsigned int i = 0; i < m_scene_objects.size(); i++)
    delete m_scene_objects[i];
    m_scene_objects.clear();
}

void KIRK::CPU::Scene::setSceneGraph(std::weak_ptr<KIRK::SceneGraph> sceneGraph)
{
    if(!sceneGraph.expired())
    {
        m_sceneGraph = sceneGraph;
    }

	if(m_sceneGraph.expired())
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
    if(!datastructure)
        throw std::invalid_argument("Don't give me a null pointer for the datastructure!");
	m_datastructure = std::move(datastructure);
    if( buildNew)
        buildDatastructure();
}

void KIRK::CPU::Scene::updateFromSceneGraph(bool updateAll)
{
    if(m_sceneGraph.expired())
    {
        throw std::invalid_argument("Empty SceneGraph");
    }

    m_camera = std::unique_ptr<KIRK::Camera>(new KIRK::Camera (*(m_sceneGraph.lock()->getActiveCamera().get())));
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
	for(std::shared_ptr<KIRK::SceneNode> child : sceneNode->m_children)
	{	//If the child is a MeshTypeObject do ...
		if(child->m_data_type == KIRK::SceneNode::MESH)
		{	//Copy the meshes dataObject

			KIRK::Mesh *mesh = (KIRK::Mesh*)child->m_data_object.get();
			//Go over every face
			for(int f = 0; f < mesh->m_faces.size(); f++)
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

		}
        else if (child->m_data_type == KIRK::SceneNode::LIGHT)
        {
            if(typeid(*(child->m_data_object.get())) == typeid(KIRK::PointLight))
            {
                m_lights.push_back(std::unique_ptr<KIRK::Light>
                                   (new KIRK::PointLight(*((KIRK::PointLight*)child->m_data_object.get()))));
                m_lights.back()->transform(child->m_data_object->calculateTransform());
            }
            else if(typeid(*(child->m_data_object.get())) == typeid(KIRK::QuadLight))
            {
                m_lights.push_back(std::unique_ptr<KIRK::Light>
                                   (new KIRK::QuadLight(*((KIRK::QuadLight*)child->m_data_object.get()))));
                m_lights.back()->transform(child->m_data_object->calculateTransform());
            }
            else if(typeid(*(child->m_data_object.get())) == typeid(KIRK::SpotLight))
            {
                m_lights.push_back(std::unique_ptr<KIRK::Light>
                                   (new KIRK::SpotLight(*((KIRK::SpotLight*)child->m_data_object.get()))));
                m_lights.back()->transform(child->m_data_object->calculateTransform());
            }
            else if(typeid(*(child->m_data_object.get())) == typeid(KIRK::SunLight))
            {
                m_lights.push_back(std::unique_ptr<KIRK::Light>
                                   (new KIRK::SunLight(*((KIRK::SunLight*)child->m_data_object.get()))));
                m_lights.back()->transform(child->m_data_object->calculateTransform());
            }
        }

        if (child ->m_children.size() > 0)
            flattenNode(child,  base_transform * child->m_transform);
	}
}

void KIRK::CPU::Scene::buildDatastructure()
{
	if(m_datastructure) m_datastructure->addBaseDataStructure(this);
}

const std::vector<KIRK::Triangle *> &KIRK::CPU::Scene::getSceneObjects() const
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
	for(unsigned int index = 0; index < m_scene_objects.size(); index++)
	{
		obj_bound = m_scene_objects[index]->getBounds();
		for(unsigned int coord = 0; coord <= 2; coord++)
		{
			if(obj_bound[0][coord] < m_bound[0][coord])
				m_bound[0][coord] = obj_bound[0][coord];
			if(obj_bound[1][coord] > m_bound[1][coord])
				m_bound[1][coord] = obj_bound[1][coord];
		}
	}
}


