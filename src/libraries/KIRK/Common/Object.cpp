#include "KIRK/Common/Object.h"

int KIRK::Object::getLongestAxis() { return m_lA; }

void KIRK::Object::setMaterial(KIRK::Material *material) {
	m_material = material;
}
KIRK::Material *KIRK::Object::getMaterial() {
	return m_material;
}

glm::vec3 *KIRK::Object::getBounds() {
	return m_bound;
}
glm::vec3 KIRK::Object::getMinBound() const {
	return m_bound[0];
}
glm::vec3 KIRK::Object::getMaxBound() const {
	return m_bound[1];
}
glm::vec3 KIRK::Object::getCentroid() const {
	return m_centroid;
}