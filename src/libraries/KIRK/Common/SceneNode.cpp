#include "SceneNode.h"

///////////////////////////////////////////////////////
///////////
///////////		SceneNode
///////////
///////////////////////////////////////////////////////

int KIRK::SceneNode::counter = 0;

void KIRK::SceneNode::assignData(std::shared_ptr<NodeDataObject> data, const Type data_type)
{
	m_data_object = data;
	m_data_type = data_type;
	m_data_object->m_parent_node = this->shared_from_this();
}

void KIRK::SceneNode::removeChild(const std::shared_ptr<SceneNode> child)
{
	child->m_parent.reset();
	m_children.erase(std::remove(m_children.begin(), m_children.end(), child), m_children.end());
}

void KIRK::SceneNode::attachChild(const std::shared_ptr<SceneNode> child)
{
	//For safety, don't attach the same Node twice.
	removeChild(child);
	m_children.emplace_back(child);

	//make this shared_ptr to a weak_ptr.
	child->m_parent = std::weak_ptr<SceneNode>(shared_from_this());
}

void KIRK::SceneNode::transformWith(glm::mat4 transform)
{
	m_transform = transform * m_transform;
}


///////////////////////////////////////////////////////
///////////
///////////		NodeDataObject
///////////
///////////////////////////////////////////////////////

glm::mat4 KIRK::NodeDataObject::calculateTransform()
{
	//Go up through node tree to calculate actual world transform
	glm::mat4 transform = glm::mat4(1.f);
	// for all nodes do stuff
	for (std::shared_ptr<SceneNode> parent = m_parent_node.lock(); parent; parent = parent->m_parent.lock())
	{
		transform = parent->m_transform * transform;
	}

	return transform;
}

void KIRK::SceneNode::onGui()
{
	ImGui::PushID(this->m_name.c_str());
	if (ImGui::TreeNode(m_name.c_str()))
	{
		ImGui::SameLine();
		ImGui::TextColored(ImVec4(0.8f, 0, 0, 1), "<------- THIS IS SELECTED");
		for (auto &child : m_children)
		{
			child->draw();
		}

		if (auto mesh = std::dynamic_pointer_cast<KIRK::Mesh>(m_data_object))
		{
			ImGui::TextColored(ImVec4(0.8f, 0, 0, 1), "Properties");
			if (ImGui::TreeNode(("Materials of " + m_name).c_str())) {
				auto materials = mesh->m_materials;
				for (auto material : materials)
				{
					material->draw();
				}
				ImGui::TreePop();
			}
		}

		if (auto light = std::dynamic_pointer_cast<KIRK::Light>(m_data_object))
		{
			ImGui::TextColored(ImVec4(0.8f, 0, 0, 1), "Properties");
			light->draw();
		}

		ImGui::TreePop();
	}
	ImGui::PopID();
}
