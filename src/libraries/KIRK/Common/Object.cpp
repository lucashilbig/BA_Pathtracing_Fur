#include "KIRK/Common/Object.h"

void KIRK::Object::setU(glm::vec3 u){
	m_u = u;
}
void KIRK::Object::setV(glm::vec3 v) {
	m_v = v;
}
void KIRK::Object::setW(glm::vec3 w) {
	m_w = w;
}

glm::vec3 KIRK::Object::getU() const {
	return m_u;
}
glm::vec3 KIRK::Object::getV() const {
	return m_v;
}
glm::vec3 KIRK::Object::getW() const {
	return m_w;
}

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