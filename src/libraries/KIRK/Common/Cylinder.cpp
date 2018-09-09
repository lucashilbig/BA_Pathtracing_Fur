#include "Cylinder.h"



KIRK::Cylinder::Cylinder(glm::vec3 basePoint, glm::vec3 apexPoint, float baseRadius, float apexRadius, glm::mat4 *modelMatrix)
{
	glm::mat3 M_ti = glm::mat3(glm::transpose(glm::inverse(*modelMatrix)));//Transposed Inverse of Modelmatrix
	//set member
	m_basepoint = glm::vec3(*modelMatrix * glm::vec4(basePoint, 1.f));
	m_apexpoint = glm::vec3(*modelMatrix * glm::vec4(apexPoint, 1.f));
	m_baseradius = baseRadius;
	m_apexradius = apexRadius;

	//pre-calculated axis of the objects local space
	m_v = apexPoint - basePoint;//v-axis(vector through the center of the cylinder(apexpoint - basepoint))
	m_height = glm::length(m_v);
	m_v = glm::normalize(m_v);

	/* find two axes which are at right angles to cone_v */
	glm::vec3 tmp(0.f, 1.f, 0.f);
	if (1.f - fabs(glm::dot(tmp, m_v)) < KIRK::cRayEpsilon)
		tmp = glm::vec3(0.f, 0.f, 1.f);

	m_u = glm::normalize(glm::cross(m_v, tmp));
	m_w = glm::normalize(glm::cross(m_u, m_v));
	//multiply with transposed inverse of modelmatrix for final local axis
	m_u = glm::normalize(M_ti * m_u);
	m_v = glm::normalize(M_ti * m_v);
	m_w = glm::normalize(M_ti * m_w);

	//calculate slope and cylinders bottom and top parameters
	m_slope = (m_baseradius - m_apexradius) / m_height;
	m_base_d = glm::dot(m_basepoint, m_v);

	m_min_d = glm::dot(m_v, m_basepoint);
	m_max_d = glm::dot(m_v, m_apexpoint);

	if (m_max_d < m_min_d)
	{
		float ftmp = m_max_d;
		m_max_d = m_min_d;
		m_min_d = ftmp;
	}

	//calculate BoundingBox and its sizeIndicator
	computeBounds();
	m_sizeIndicator = m_bound[1] - m_bound[0];

	//calculate centroid, which is located at the middle of the cylinder axis(v-axis)
	m_centroid = m_basepoint + (0.5 * m_v);

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

}

KIRK::Cylinder::~Cylinder()
{
}

bool KIRK::Cylinder::closestIntersection(Intersection *hit, float tMin, float tMax)
{
	//ray direction normalized
	float a, b, c, d;
	int nroots;
	bool enter;
	Ray *ray = &hit->m_ray;
	glm::vec3 Q;

	// Transformation of ray in local coordinate system

	glm::vec3 P = ray->m_origin - m_basepoint;
	glm::vec3 dir = ray->m_direction; //d not necessarily normalized!

	P = glm::vec3(glm::dot(P, m_u), glm::dot(P, m_v), glm::dot(P, m_w));
	glm::vec3 D(glm::dot(dir, m_u), glm::dot(dir, m_v), glm::dot(dir, m_w));

	a = 1.f - D.y*D.y * (1.f + m_slope * m_slope);

	b = P.x * D.x + P.z * D.z
		+ m_baseradius * m_slope * D.y
		- m_slope * m_slope * P.y * D.y;

	c = m_baseradius - m_slope * P.y;
	c = P.x * P.x + P.z * P.z - c * c;

	float disc = b * b - a * c;

	if (disc < 0.0)
		return false;

	disc = sqrtf(disc);
	float t1 = (-b - disc) / a; //d not necessarily normalized!
	float t2 = (-b + disc) / a; //d not necessarily normalized!

	if ((t2 < tMin) || (t1 > tMax))
		return false;
	if (t1 < KIRK::cRayEpsilon)
		nroots = 1;
	else
		nroots = 2;

	// ensure that the points are between the two bounding planes...

	switch (nroots) {
	case 1:
		if ((t2 > tMax) || (t2 < tMin))
			return false;
		Q = ray->followDistance(t2);
		d = glm::dot(m_v, Q);
		if (d >= m_min_d && d <= m_max_d)
		{
			enter = false;
			hit->update(this, t2, enter);
			return true;
		}
		else
			return false;
		break;
	case 2:
		if ((t1 < tMin) && (t2 > tMax))
			return false;
		Q = ray->followDistance(t1);
		d = glm::dot(m_v, Q);
		if (d >= m_min_d && d <= m_max_d)
		{
			enter = true;
			hit->update(this, t1, enter);
			return true;
		}
		else {
			Q = ray->followDistance(t2);
			d = glm::dot(m_v, Q);
			if (d >= m_min_d && d <= m_max_d)
			{
				enter = false;
				hit->update(this, t2, enter);
				return true;
			}
		}
		return false;
	}
	return false;
}

bool KIRK::Cylinder::isIntersection(Ray* ray, float tMax)
{
	//ray direction NOT normalized
	float a, b, c, d;
	int nroots;
	glm::vec3 Q;

	// Transformation of ray in local coordinate system

	glm::vec3 P = ray->m_origin - m_basepoint;
	glm::vec3 dir = ray->m_direction; //d not necessarily normalized!

	P = glm::vec3(glm::dot(P, m_u), glm::dot(P, m_v), glm::dot(P, m_w));
	glm::vec3 D(glm::dot(dir, m_u), glm::dot(dir, m_v), glm::dot(dir, m_w));

	a = D.x*D.x + D.z*D.z - m_slope * m_slope*D.y*D.y;

	b = P.x * D.x + P.z * D.z
		+ m_baseradius * m_slope * D.y
		- m_slope * m_slope * P.y * D.y;

	c = m_baseradius - m_slope * P.y;
	c = P.x * P.x + P.z * P.z - c * c;

	float disc = b * b - a * c;

	if (disc < 0.0)
		return false;

	disc = sqrtf(disc);
	float t1 = (-b - disc) / a; //d not necessarily normalized!
	float t2 = (-b + disc) / a; //d not necessarily normalized!

	if ((t2 < 0.f) || (t1 > tMax))
		return false;
	if (t1 < KIRK::cRayEpsilon)
		nroots = 1;
	else
		nroots = 2;

	// ensure that the points are between the two bounding planes...

	switch (nroots) {
	case 1:
		if ((t2 > tMax) || (t2 < 0.f))
			return false;
		Q = ray->followDistance(t2);
		d = glm::dot(m_v, Q);
		if (d >= m_min_d && d <= m_max_d)
			return true;
		else
			return false;
		break;
	case 2:
		if ((t1 < 0.f) && (t2 > tMax))
			return false;
		Q = ray->followDistance(t1);
		d = glm::dot(m_v, Q);
		if (d >= m_min_d && d <= m_max_d)
			return true;
		else
		{
			Q = ray->followDistance(t2);
			d = glm::dot(m_v, Q);
			if (d >= m_min_d && d <= m_max_d)
				return true;
		}
		return false;
	}
	return false;
}

void KIRK::Cylinder::calcNormal(Intersection *hit)
{
	glm::vec3 Q = hit->m_location;
	float t = glm::dot(Q, m_v) - m_base_d;
	glm::vec3 q_1 = Q - t * m_v;
	glm::vec3 n = glm::normalize(q_1 - m_basepoint);
	hit->m_normal = glm::normalize(n + m_slope * m_v);
}

void KIRK::Cylinder::calcTcoord(Intersection *hit)
{
	float phi;

	glm::vec3 Q = hit->m_location - m_basepoint;
	float u = glm::dot(Q, m_u);
	float v = glm::dot(Q, m_v);
	float w = glm::dot(Q, m_w);
	float r = m_baseradius - m_slope * v;

	float tmp = w / r;
	if (tmp < -1.f)
		tmp = -1.f;
	else if (tmp > 1.f)
		tmp = 1.f;

	if (u < 0)
		phi = 2.0f * glm::pi<float>() - acos(tmp);
	else
		phi = acos(tmp);
	hit->m_texcoord = glm::vec2(phi / 2.f / glm::pi<float>(), v / m_height);
}

bool KIRK::Cylinder::isInAABB(glm::vec3 *bbox)
{
	glm::vec3 P[8];
	float min = FLT_MAX, max = -FLT_MAX;

	//Transform cube points in cone coordiante system and 
	//project cube points to v-axis 
	glm::mat3 cone_matrix(glm::transpose(glm::mat3(m_u, m_v, m_w)));
	P[0] = cone_matrix * (glm::vec3(bbox[0].x, bbox[1].y, bbox[1].z) - m_basepoint);
	P[1] = cone_matrix * (glm::vec3(bbox[0].x, bbox[0].y, bbox[1].z) - m_basepoint);
	P[2] = cone_matrix * (glm::vec3(bbox[1].x, bbox[0].y, bbox[1].z) - m_basepoint);
	P[3] = cone_matrix * (glm::vec3(bbox[1].x, bbox[1].y, bbox[1].z) - m_basepoint);
	P[4] = cone_matrix * (glm::vec3(bbox[1].x, bbox[1].y, bbox[0].z) - m_basepoint);
	P[5] = cone_matrix * (glm::vec3(bbox[1].x, bbox[0].y, bbox[0].z) - m_basepoint);
	P[6] = cone_matrix * (glm::vec3(bbox[0].x, bbox[0].y, bbox[0].z) - m_basepoint);
	P[7] = cone_matrix * (glm::vec3(bbox[0].x, bbox[1].y, bbox[0].z) - m_basepoint);

	for (unsigned int i = 0; i < 8; i++)
	{
		if (P[i].y < min) min = P[i].y;
		if (P[i].y > max) max = P[i].y;
	}
	if ((min > m_height) || (max < 0.f))
		return false;

	// compare radius of P[i] with the radius of the cone
	int n = 0;
	for (unsigned int i = 0; i < 8; i++)
	{
		float r = m_baseradius - m_slope * P[i].y;
		float r_p = sqrt(P[i].x*P[i].x + P[i].z*P[i].z);
		if ((r < r_p) || (P[i].y > m_height) || (P[i].y < 0.f))
			n++;
	}
	if (n == 8)
		return false;
	//Special case (box surrounds cone) ignored, since voxels are assumed to be smaller
	if (n > 0)
		return true;

	return false;
}


void KIRK::Cylinder::computeBounds()
{
	glm::vec3 lbound[2];
	glm::vec3 P[8];

	float radius = (m_baseradius > m_apexradius) ? m_baseradius + 1e-6f : m_apexradius + 1e-6f;
	lbound[0] = glm::vec3(-radius, 0, -radius);
	lbound[1] = glm::vec3(radius, m_height, radius);
	glm::mat3 cone_matrix = glm::mat3(m_u, m_v, m_w);

	P[0] = cone_matrix * glm::vec3(lbound[0].x, lbound[1].y, lbound[1].z) + m_basepoint;
	P[1] = cone_matrix * glm::vec3(lbound[0].x, lbound[0].y, lbound[1].z) + m_basepoint;
	P[2] = cone_matrix * glm::vec3(lbound[1].x, lbound[0].y, lbound[1].z) + m_basepoint;
	P[3] = cone_matrix * glm::vec3(lbound[1].x, lbound[1].y, lbound[1].z) + m_basepoint;
	P[4] = cone_matrix * glm::vec3(lbound[1].x, lbound[1].y, lbound[0].z) + m_basepoint;
	P[5] = cone_matrix * glm::vec3(lbound[1].x, lbound[0].y, lbound[0].z) + m_basepoint;
	P[6] = cone_matrix * glm::vec3(lbound[0].x, lbound[0].y, lbound[0].z) + m_basepoint;
	P[7] = cone_matrix * glm::vec3(lbound[0].x, lbound[1].y, lbound[0].z) + m_basepoint;

	m_bound[0] = glm::vec3(FLT_MAX);
	m_bound[1] = glm::vec3(-FLT_MAX);
	for (unsigned int i = 0; i < 8; i++)
	{
		if (P[i].x < m_bound[0].x) m_bound[0].x = P[i].x;
		if (P[i].x > m_bound[1].x) m_bound[1].x = P[i].x;
		if (P[i].y < m_bound[0].y) m_bound[0].y = P[i].y;
		if (P[i].y > m_bound[1].y) m_bound[1].y = P[i].y;
		if (P[i].z < m_bound[0].z) m_bound[0].z = P[i].z;
		if (P[i].z > m_bound[1].z) m_bound[1].z = P[i].z;
	}
}
