#include "Cylinder.h"



KIRK::Cylinder::Cylinder(glm::vec3 basePoint, glm::vec3 apexPoint, float baseRadius, float apexRadius, glm::mat4 transform_objectToWolrd, glm::mat4 transform_worldToObject)
{	
	//set member
	m_basepoint = basePoint;
	m_apexpoint = apexPoint;
	m_baseradius = baseRadius;
	m_apexradius = apexRadius;
	trans_o2w = transform_objectToWolrd;
	trans_w2o = transform_worldToObject;
	//calculate zmin, zmax in local space
	glm::vec4 local_bp = trans_w2o * glm::vec4(m_basepoint, 1.0f);//base point in local object space
	glm::vec4 local_ap = trans_w2o * glm::vec4(m_apexpoint, 1.0f);//apex point in local object space
	m_zmin = glm::min(local_bp.z, local_ap.z);
	m_zmax = glm::max(local_bp.z, local_ap.z);

	//calculate BoundingBox
	float max_radius = glm::max(m_baseradius, m_apexradius);
	m_bound[0] = glm::vec3(-max_radius, -max_radius, m_zmin) - glm::vec3(cCylinderEpsilon);
	m_bound[1] = glm::vec3(max_radius, max_radius, m_zmax) + glm::vec3(cCylinderEpsilon);

	//calculate longest axis of the bounding box
	glm::vec3 boundDiff = m_bound[1] - m_bound[0];
	m_lA = 0;
	float longestAxisDiff = boundDiff.x;
	if (boundDiff.y > longestAxisDiff)
	{
		longestAxisDiff = boundDiff.y;
		m_lA = 1;
	}
	if (boundDiff.z > longestAxisDiff)
	{
		longestAxisDiff = boundDiff.z;
		m_lA = 2;
	}

	//calculate centroid
	glm::vec3 bpap = glm::normalize(m_basepoint - m_apexpoint);
	m_centroid = m_basepoint + (0.5 * bpap);

}

KIRK::Cylinder::Cylinder(glm::vec3 basePoint, glm::vec3 apexPoint, float baseRadius, float apexRadius)
{
	//set member
	m_basepoint = basePoint;
	m_apexpoint = apexPoint;
	m_baseradius = baseRadius;
	m_apexradius = apexRadius;

	//calculate the transformation matrix for world to object transformation
	glm::vec3 bpap = glm::normalize(m_basepoint - m_apexpoint);//Vector between base and apex point
	glm::vec3 rot_vec = glm::normalize(glm::cross(glm::vec3(0.0f, 0.0f, 1.0f), bpap));//Rotation Vector(Axis) for local object space
	float rot_angle = glm::dot(glm::vec3(0.0f,0.0f,1.0f), bpap);//Rotation angle for new z-axis of the local object space
	glm::mat4 w2o = glm::rotate(glm::mat4(1.0f), rot_angle, rot_vec);//rotation of the transformation matrix around the calculated rotation vector and angle
	w2o = glm::translate(w2o, m_basepoint);//translation of the world to object trans matrix
	//set the transformation matrices
	trans_w2o = w2o;
	trans_o2w = glm::inverse(trans_w2o);

	//calculate zmin, zmax in local space
	glm::vec4 local_bp = trans_w2o * glm::vec4(m_basepoint, 1.0f);//base point in local object space
	glm::vec4 local_ap = trans_w2o * glm::vec4(m_apexpoint, 1.0f);//apex point in local object space
	m_zmin = glm::min(local_bp.z, local_ap.z);
	m_zmax = glm::max(local_bp.z, local_ap.z);

	//calculate BoundingBox
	float max_radius = glm::max(m_baseradius, m_apexradius);
	m_bound[0] = glm::vec3(-max_radius, -max_radius, m_zmin) - glm::vec3(cCylinderEpsilon);
	m_bound[1] = glm::vec3(max_radius, max_radius, m_zmax) + glm::vec3(cCylinderEpsilon);

	//calculate longest axis of the bounding box
	glm::vec3 boundDiff = m_bound[1] - m_bound[0];
	m_lA = 0;
	float longestAxisDiff = boundDiff.x;
	if (boundDiff.y > longestAxisDiff)
	{
		longestAxisDiff = boundDiff.y;
		m_lA = 1;
	}
	if (boundDiff.z > longestAxisDiff)
	{
		longestAxisDiff = boundDiff.z;
		m_lA = 2;
	}

	//calculate centroid
	m_centroid = m_basepoint + (0.5 * bpap);
}

