#include "GuiScene.h"

namespace KIRK
{
GuiScene::GuiScene(std::shared_ptr<KIRK::Gui> gui)
	: m_gui(gui), m_flags(UPDATE_MATERIALS)
{
}

GuiScene::~GuiScene()
{
}

void GuiScene::buildSceneGui(std::shared_ptr<SceneGraph> graph, CVK::CVKCameraSynchronizer &sync)
{
	std::map<std::string, std::shared_ptr<Material>> materials;

	bool mat_open = false;
	if (auto mat_win = m_gui->get<GuiWindow>("Materials")) {
		mat_open = mat_win->isEnabled();
		mat_win->detachChildren();
	}
	else {
		m_gui->make<GuiWindow>("Materials", mat_open);
	}

	bool graph_open = false;
	if(auto graph_win = m_gui->get<GuiWindow>("Scene"))
	{
		graph_open = graph_win->isEnabled();
		graph_win->detachChildren();
	} else
	{
		m_gui->make<GuiWindow>("Scene", graph_open);
	}

	m_gui->get<GuiWindow>("Scene")->make<GuiNamedLambdaElement>("ThisTotallyWorks...Now", [&, graph]()
	{
		std::vector<const char*> texts;
		for(int i=0; i<m_cameras.size(); ++i)
		{
			texts.push_back(m_cameras[i]->m_parent_node.lock()->m_name.c_str());
		}

		int last = m_current_camera;
		ImGui::Combo("Active Cam", &m_current_camera, texts.data(), texts.size());

		if(last != m_current_camera)
		{
			graph->setActiveCamera(m_cameras[m_current_camera]);
			sync.setCVKCamera(graph, sync.getCVKCamera());
			sync.update(m_gui->getWindow());
		}

		if(ImGui::Button("New Camera"))
		{
			graph->createDefaultCamera();
			m_cameras.clear();
			for (auto &&node : *graph)
			{
				if (node->m_data_type == KIRK::SceneNode::CAMERA)
				{
					m_cameras.push_back(std::dynamic_pointer_cast<KIRK::Camera>(node->m_data_object));
				}
			}
		}
	});
	m_gui->get<GuiWindow>("Scene")->add(graph->getRootNode().get());

	bool lights_open = false;
	if (auto lights_win = m_gui->get<GuiWindow>("Lights"))
	{
		lights_open = lights_win->isEnabled();
		lights_win->detachChildren();
	}
	else
	{
		m_gui->make<GuiWindow>("Lights", lights_open);
	}

	for (auto &&node : *graph)
	{
		if(node->m_data_type == KIRK::SceneNode::MESH)
		{
			auto mesh = std::dynamic_pointer_cast<KIRK::Mesh>(node->m_data_object);
			for (auto &&material : mesh->m_materials) {
				materials.emplace(material->name, material);
			}
		} else if (node->m_data_type == KIRK::SceneNode::LIGHT)
		{
			auto light = std::dynamic_pointer_cast<KIRK::Light>(node->m_data_object);
			m_gui->get<GuiWindow>("Lights")->add(light.get());
		} else if(node->m_data_type == KIRK::SceneNode::CAMERA)
		{
			m_cameras.push_back(std::dynamic_pointer_cast<KIRK::Camera>(node->m_data_object));
		}
	}

	for (auto &&mat_pair : materials) {
		m_gui->get<GuiWindow>("Materials")->add(mat_pair.second.get());
	}

	if (m_editable_flags && m_update_callback) {
		auto upd_win = m_gui->get<GuiWindow>("Scene Updater");
		bool upd_open = false;
		if (upd_win) {
			upd_open = upd_win->isEnabled();
			upd_win->detachChildren();
		}
		else {
			upd_win = m_gui->make<GuiWindow>("Scene Updater", upd_open);
		}

		upd_win->make<GuiNamedLambdaElement>("Scene Update Controls", [&]()
		{
			if ((m_editable_flags & SceneUpdateFlags::UPDATE_MATERIALS) == SceneUpdateFlags::UPDATE_MATERIALS)
				ImGui::CheckboxFlags("Update Materials", &m_flags, SceneUpdateFlags::UPDATE_MATERIALS);

			if ((m_editable_flags & SceneUpdateFlags::UPDATE_GRAPH) == SceneUpdateFlags::UPDATE_GRAPH)
				ImGui::CheckboxFlags("Update Graph", &m_flags, SceneUpdateFlags::UPDATE_GRAPH);

			if (ImGui::Button("Update", ImVec2(ImGui::GetContentRegionAvailWidth(), 32)))
			{
				if (m_update_callback)
					m_update_callback(m_flags);
			}
		});
	}
}

void GuiScene::setUpdateCallback(std::function<void(unsigned)> callback, unsigned not_editable)
{
	m_update_callback = callback;
	m_editable_flags = ~not_editable;
}
}