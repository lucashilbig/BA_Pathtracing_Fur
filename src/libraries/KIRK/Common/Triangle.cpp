#include "Triangle.h"

KIRK::Triangle::Triangle(glm::vec3 a, glm::vec3 b, glm::vec3 c, glm::vec3 na, glm::vec3 nb, glm::vec3 nc,
                           glm::vec2 tca, glm::vec2 tcb, glm::vec2 tcc) : i0(0), i1(0), i2(0) 
{
	glm::mat4 *ModelMatrix = new glm::mat4(1.0f);

    glm::mat3 M_ti = glm::mat3(glm::transpose(glm::inverse(*ModelMatrix)));

    a = glm::vec3(*ModelMatrix * glm::vec4(a, 1.f));
    b = glm::vec3(*ModelMatrix * glm::vec4(b, 1.f));
    c = glm::vec3(*ModelMatrix * glm::vec4(c, 1.f));

    m_bound[0] = glm::min(glm::min(a, b), c) - glm::vec3(cRayEpsilon);
    m_bound[1] = glm::max(glm::max(a, b), c) + glm::vec3(cRayEpsilon);

    glm::vec3 boundDiff = m_bound[1] - m_bound[0];
    m_lA = 0;
    float longestAxisDiff = boundDiff.x;
    if(boundDiff.y > longestAxisDiff)
    {
        longestAxisDiff = boundDiff.y;
        m_lA = 1;
    }
    if(boundDiff.z > longestAxisDiff)
    {
        longestAxisDiff = boundDiff.z;
        m_lA = 2;
    }

    m_A = a;
    m_B = b;
    m_C = c;
    m_na = glm::normalize(M_ti * na);
    m_nb = glm::normalize(M_ti * nb);
    m_nc = glm::normalize(M_ti * nc);
    m_tca = tca;
    m_tcb = tcb;
    m_tcc = tcc;

    if(a[m_lA] <= b[m_lA] && b[m_lA] <= c[m_lA])
    {
        m_A = a;
        m_B = b;
        m_C = c;
        m_na = glm::normalize(M_ti * na);
        m_nb = glm::normalize(M_ti * nb);
        m_nc = glm::normalize(M_ti * nc);
        m_tca = tca;
        m_tcb = tcb;
        m_tcc = tcc;
    }
    if(b[m_lA] <= a[m_lA] && a[m_lA] <= c[m_lA])
    {
        m_A = b;
        m_B = a;
        m_C = c;
        m_na = glm::normalize(M_ti * nb);
        m_nb = glm::normalize(M_ti * na);
        m_nc = glm::normalize(M_ti * nc);
        m_tca = tcb;
        m_tcb = tca;
        m_tcc = tcc;
    }
    if(a[m_lA] <= c[m_lA] && c[m_lA] <= b[m_lA])
    {
        m_A = a;
        m_B = c;
        m_C = b;
        m_na = glm::normalize(M_ti * na);
        m_nb = glm::normalize(M_ti * nc);
        m_nc = glm::normalize(M_ti * nb);
        m_tca = tca;
        m_tcb = tcc;
        m_tcc = tcb;
    }
    if(c[m_lA] <= a[m_lA] && a[m_lA] <= b[m_lA])
    {
        m_A = c;
        m_B = a;
        m_C = b;
        m_na = glm::normalize(M_ti * nc);
        m_nb = glm::normalize(M_ti * na);
        m_nc = glm::normalize(M_ti * nb);
        m_tca = tcc;
        m_tcb = tca;
        m_tcc = tcb;
    }
    if(b[m_lA] <= c[m_lA] && c[m_lA] <= a[m_lA])
    {
        m_A = b;
        m_B = c;
        m_C = a;
        m_na = glm::normalize(M_ti * nb);
        m_nb = glm::normalize(M_ti * nc);
        m_nc = glm::normalize(M_ti * na);
        m_tca = tcb;
        m_tcb = tcc;
        m_tcc = tca;
    }
    if(c[m_lA] <= b[m_lA] && b[m_lA] <= a[m_lA])
    {
        m_A = c;
        m_B = b;
        m_C = a;
        m_na = glm::normalize(M_ti * nc);
        m_nb = glm::normalize(M_ti * nb);
        m_nc = glm::normalize(M_ti * na);
        m_tca = tcc;
        m_tcb = tcb;
        m_tcc = tca;
    }

    m_ab = m_B - m_A;
    m_ac = m_C - m_A;
    m_bc = m_C - m_B;

    if(m_ab[m_lA] == 0.0f)
        m_ab[m_lA] = 0.0001f;
    if(m_ac[m_lA] == 0.0f)
        m_ac[m_lA] = 0.0001f;
    if(m_bc[m_lA] == 0.0f)
        m_bc[m_lA] = 0.0001f;

    m_Normal = glm::normalize((m_na + m_nb + m_nc) / 3.0f);

    m_sizeIndicator = m_bound[1] - m_bound[0];
    m_centroid = (m_A + m_B + m_C) / 3.f;
}

KIRK::Triangle::~Triangle()
{
}

bool KIRK::Triangle::testAxis(const glm::vec3 e, const glm::vec3 v[], const glm::vec3 boxhalfsize, const int i, const int j, const int vi1, const int vi2, const int f, float &min, float &max, float &p0, float &p2, float &rad) {
	p0 = f * (e[i] * v[vi1][j] - e[j] * v[vi1][i]);
	p2 = f * (e[i] * v[vi2][j] - e[j] * v[vi2][i]);
	if (p0 < p2)
	{
		min = p0; max = p2;
	}
	else
	{
		min = p2; max = p0;
	}
	rad = fabsf(e[i]) * boxhalfsize[j] + fabsf(e[j]) * boxhalfsize[i];
	if (min > rad || max < -rad)
		return false;
	return true;
}

bool KIRK::Triangle::closestIntersection(Intersection *hit, float tMin, float tMax)
{
    glm::vec3 dir = hit->m_ray.m_direction;        // dir has length 1

    glm::vec3 d_v = glm::cross(dir, m_ac);
    float det = glm::dot(d_v, m_ab);
    if(glm::abs(det) < cTriangleEpsilon)
        return false;

    float invDet = 1.f / det;

    glm::vec3 P = hit->m_ray.m_origin;
    glm::vec3 w = P - m_A;

    float u = glm::dot(d_v, w) * invDet;
    if(u < 0.f || u > 1.f)
        return false;

    glm::vec3 w_u = glm::cross(w, m_ab);

    float v = glm::dot(w_u, dir) * invDet;
    if(v < 0.f || u + v > 1.f)
        return false;

    float t = glm::dot(w_u, m_ac) * invDet;
    if((t < tMin) || (t > tMax))
        return false;

    bool enter = glm::dot(m_Normal, dir) < 0 ? true : false;
    hit->update(this, t, enter, glm::vec3(1 - u - v, u, v));

    return true;
}

bool KIRK::Triangle::closestIntersectionAsPlane(Intersection *hit, float tMin, float tMax)
{
    glm::vec3 dir = hit->m_ray.m_direction;        // dir has length 1

    glm::vec3 d_v = glm::cross(dir, m_ac);
    float det = glm::dot(d_v, m_ab);
    if(glm::abs(det) < cTriangleEpsilon)
        return false;

    float invDet = 1.f / det;

    glm::vec3 P = hit->m_ray.m_origin;
    glm::vec3 w = P - m_A;

    float u = glm::dot(d_v, w) * invDet;
    glm::vec3 w_u = glm::cross(w, m_ab);
    float v = glm::dot(w_u, dir) * invDet;
    float t = glm::dot(w_u, m_ac) * invDet;
    if((t < tMin) || (t > tMax))
        return false;

    bool enter = glm::dot(m_Normal, dir) < 0 ? true : false;
    hit->update(this, t, enter, glm::vec3(1 - u - v, u, v));

    return true;
}

bool KIRK::Triangle::isIntersection(Ray *ray, float tMax)
{
    glm::vec3 dir = ray->m_direction;

    glm::vec3 d_v = glm::cross(dir, m_ac);
    float det = glm::dot(d_v, m_ab);
    if(glm::abs(det) < cTriangleEpsilon)
        return false;

    float invDet = 1.f / det;

    glm::vec3 P = ray->m_origin;
    glm::vec3 w = P - m_A;

    float u = glm::dot(d_v, w) * invDet;
    if(u < 0.f || u > 1.f)
        return false;

    glm::vec3 w_u = glm::cross(w, m_ab);

    float v = glm::dot(w_u, dir) * invDet;
    if(v < 0.f || u + v > 1.f)
        return false;

    float t = glm::dot(w_u, m_ac) * invDet;
    if((t < 0.f) || (t > tMax))
        return false;

    return true;
}

void KIRK::Triangle::calcNormal(Intersection *hit)
{
    glm::vec3 bc = hit->m_barycentric_coord;
    hit->m_normal = glm::normalize(bc.x * m_na + bc.y * m_nb + bc.z * m_nc);
}

void KIRK::Triangle::calcTcoord(Intersection *hit)
{
    glm::vec3 bc = hit->m_barycentric_coord;
    hit->m_texcoord = bc.x * m_tca + bc.y * m_tcb + bc.z * m_tcc;
}

void KIRK::Triangle::computeBounds()
{
    glm::vec3 Vertex[2];
    Vertex[0] = m_A + m_ab;
    Vertex[1] = m_A + m_ac;

    m_bound[0] = m_A;
    m_bound[1] = m_A;

    for(int v = 0; v <= 1; v++)
    {
        for(unsigned int coord = 0; coord <= 2; coord++)
        {
            if(Vertex[v][coord] < m_bound[0][coord])
                m_bound[0][coord] = Vertex[v][coord];
            if(Vertex[v][coord] > m_bound[1][coord])
                m_bound[1][coord] = Vertex[v][coord];
        }
    }
    m_bound[0] -= glm::vec3(cRayEpsilon, cRayEpsilon, cRayEpsilon);
    m_bound[1] += glm::vec3(cRayEpsilon, cRayEpsilon, cRayEpsilon);
}

bool KIRK::Triangle::planeBoxOverlap(const glm::vec3 &n, const glm::vec3 &v0, const glm::vec3 &halfsize) const
{
    glm::vec3 vMin, vMax;
    float v;
    for(int i = 0; i <= 2; i++)
    {
        v = v0[i];
        if(n[i] > 0.0f)
        {
            vMin[i] = -halfsize[i] - v;
            vMax[i] = halfsize[i] - v;
        } else
        {
            vMin[i] = halfsize[i] - v;
            vMax[i] = -halfsize[i] - v;
        }
    }

    if(glm::dot(n, vMin) > 0.0f)
        return false;
    if(glm::dot(n, vMax) >= 0.0f)
        return true;

    return false;
}

bool KIRK::Triangle::isInAABB(glm::vec3 *bbox)
{
    glm::vec3 diff = (bbox[1] - bbox[0]);
    glm::vec3 boxhalfsize = diff * 0.5f;
    glm::vec3 mid = bbox[0] + diff * 0.5f;
    glm::vec3 v[3] = {(m_A - mid), (m_B - mid), (m_C - mid)};

    glm::vec3 e0 = m_ab;
    glm::vec3 e1 = m_bc;
    glm::vec3 e2 = -m_ac;

    float min, max, p0, p2, rad;

	if (!testAxis(e0, v, boxhalfsize, 2, 1, 0, 2, 1, min, max, p0, p2, rad)) return false;
	if (!testAxis(e0, v, boxhalfsize, 2, 0, 0, 2, -1, min, max, p0, p2, rad)) return false;
	if (!testAxis(e0, v, boxhalfsize, 1, 0, 1, 2, 1, min, max, p0, p2, rad)) return false;

	if (!testAxis(e1, v, boxhalfsize, 2, 1, 0, 2, 1, min, max, p0, p2, rad)) return false;
	if (!testAxis(e1, v, boxhalfsize, 2, 0, 0, 2, -1, min, max, p0, p2, rad)) return false;
	if (!testAxis(e1, v, boxhalfsize, 1, 0, 0, 1, 1, min, max, p0, p2, rad)) return false;

	if (!testAxis(e2, v, boxhalfsize, 2, 1, 0, 1, 1, min, max, p0, p2, rad)) return false;
	if (!testAxis(e2, v, boxhalfsize, 2, 0, 0, 1, -1, min, max, p0, p2, rad)) return false;
	if (!testAxis(e2, v, boxhalfsize, 1, 0, 1, 2, 1, min, max, p0, p2, rad)) return false;

	for (int i = 0; i < 2; i++) {
		min = max = v[0][i];
		if (v[1][i] < min) min = v[1][i];
		if (v[2][i] < min) min = v[2][i];
		if (v[1][i] > max) max = v[1][i];
		if (v[2][i] > max) max = v[2][i];
		if (min>boxhalfsize[i] || max<-boxhalfsize[i])
			return false;
	}

    glm::vec3 n = glm::cross(e0, e1);

    if(!planeBoxOverlap(n, v[0], boxhalfsize))
        return false;

    return true;
}
