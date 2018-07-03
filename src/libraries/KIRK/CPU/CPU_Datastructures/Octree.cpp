#include "KIRK/Utils/Log.h"
#include "Octree.h"

int KIRK::CPU::Octree::s_maxRecursionDepth = 0;
glm::vec3 *KIRK::CPU::Octree::s_NodeRadius = new glm::vec3[12];

KIRK::CPU::Octree::Octree(glm::vec3 minBound, glm::vec3 maxBound, int maxRecursionDepth)
        : CPU_DataStructure(minBound, maxBound)
{
    s_maxRecursionDepth = maxRecursionDepth;

    s_NodeRadius[0] = m_size / 2.f; //radius is half of size
    for(int i = 0; i < maxRecursionDepth; i++)
        s_NodeRadius[i + 1] = s_NodeRadius[i] / 2.f;

    m_halfSize = m_size / 2.0f;
    m_mid = m_minBound + m_halfSize;
    m_recursionDepth = 0;
}

KIRK::CPU::Octree::Octree(int depth, glm::vec3 mid) : CPU_DataStructure(mid - s_NodeRadius[depth],
                                                                     mid + s_NodeRadius[depth])
{

    m_halfSize = m_size / 2.0f;
    m_mid = mid;
    m_recursionDepth = depth;
}

KIRK::CPU::Octree::~Octree()
{
}

void KIRK::CPU::Octree::printDebugInfo()
{
	LOG_DEBUG("Octree with depth: %", s_maxRecursionDepth);
	LOG_DEBUG("Memory size: % bytes.", getSizeInBytes());
	LOG_DEBUG("# of Nodes: % , # of Leafs %", numTotalChildren(), numTotalLeafs());
}

int KIRK::CPU::Octree::getSizeInBytes()
{
    int result = sizeof(*this);
    result += m_candidateList.size() * sizeof(KIRK::Triangle *);
    result += m_children.size() * sizeof(Octree *);
    for(Octree *c : m_children)
        result += c->getSizeInBytes();
    return result;
}


void KIRK::CPU::Octree::addBaseDataStructure(KIRK::CPU::Scene *scene)
{
    for(KIRK::Object *obj : scene->getSceneObjects())
        addObject(obj);
}


void KIRK::CPU::Octree::subdivide()
{
    m_children.reserve(8);

    int nextDepth = m_recursionDepth + 1;

    glm::vec3 offset = s_NodeRadius[nextDepth];

    addChild(new Octree(nextDepth, (m_mid + glm::vec3(-offset.x, -offset.y, -offset.z))));
    addChild(new Octree(nextDepth, (m_mid + glm::vec3(-offset.x, -offset.y, offset.z))));
    addChild(new Octree(nextDepth, (m_mid + glm::vec3(-offset.x, offset.y, -offset.z))));
    addChild(new Octree(nextDepth, (m_mid + glm::vec3(-offset.x, offset.y, offset.z))));
    addChild(new Octree(nextDepth, (m_mid + glm::vec3(offset.x, -offset.y, -offset.z))));
    addChild(new Octree(nextDepth, (m_mid + glm::vec3(offset.x, -offset.y, offset.z))));
    addChild(new Octree(nextDepth, (m_mid + glm::vec3(offset.x, offset.y, -offset.z))));
    addChild(new Octree(nextDepth, (m_mid + glm::vec3(offset.x, offset.y, offset.z))));
}

void KIRK::CPU::Octree::addObject(KIRK::Object *obj)
{
    glm::vec3 bounds[2] = {m_minBound, m_maxBound};
    if(!(obj->isInAABB(bounds)))
        return;

    if(m_recursionDepth < s_maxRecursionDepth)
    {
        if(!hasChildren())
            subdivide();
        for(Octree *oct : m_children)
            oct->addObject(obj);
    } else
    {
        addCandidate(obj);
    }
}

bool KIRK::CPU::Octree::closestIntersection(KIRK::Intersection *hit)
{
    glm::vec3 direction = hit->m_ray.m_direction;
    glm::vec3 origin = hit->m_ray.m_origin;
    glm::vec3 tnext[2];

    unsigned char directionBits = 0;

    if(direction.x < 0.f)
    {
        origin.x = m_minBound.x + m_maxBound.x - origin.x;
        direction.x = -(direction.x);
        directionBits |= 4;
    } else if(direction.x == 0.f)
        direction.x = cRayEpsilon;

    if(direction.y < 0.f)
    {
        origin.y = m_minBound.y + m_maxBound.y - origin.y;
        direction.y = -(direction.y);
        directionBits |= 2;
    } else if(direction.y == 0.f)
        direction.y = cRayEpsilon;

    if(direction.z < 0.f)
    {
        origin.z = m_minBound.z + m_maxBound.z - origin.z;
        direction.z = -(direction.z);
        directionBits |= 1;
    } else if(direction.z == 0.f)
        direction.z = cRayEpsilon;

    tnext[0] = (m_minBound - origin) / direction;
    tnext[1] = (m_maxBound - origin) / direction;

    float tmin = glm::compMax(tnext[0]);
    float tmax = glm::compMin(tnext[1]);

    bool intersectionFound = false;

    if(tmin < tmax)
        intersectionFound |= traverseNode(&tnext[0], &tnext[1], directionBits, hit);

    return intersectionFound;
}


bool KIRK::CPU::Octree::isIntersection(KIRK::Ray *ray, float tMax)
{

    glm::vec3 direction = ray->m_direction;
    glm::vec3 origin = ray->m_origin;
    glm::vec3 tnext[2];

    unsigned char directionBits = 0;

    if(direction.x < 0.f)
    {
        origin.x = m_minBound.x + m_maxBound.x - origin.x;
        direction.x = -(direction.x);
        directionBits |= 4;
    } else if(direction.x == 0.f)
        direction.x = cRayEpsilon;

    if(direction.y < 0.f)
    {
        origin.y = m_minBound.y + m_maxBound.y - origin.y;
        direction.y = -(direction.y);
        directionBits |= 2;
    } else if(direction.y == 0.f)
        direction.y = cRayEpsilon;

    if(direction.z < 0.f)
    {
        origin.z = m_minBound.z + m_maxBound.z - origin.z;
        direction.z = -(direction.z);
        directionBits |= 1;
    } else if(direction.z == 0.f)
        direction.z = cRayEpsilon;

    tnext[0] = (m_minBound - origin) / direction;
    tnext[1] = (m_maxBound - origin) / direction;

    float tmin = glm::compMax(tnext[0]);

    if((tmin < tMax))
        return traverseNode(&tnext[0], &tnext[1], directionBits, tMax, ray);
    return false;
}


bool KIRK::CPU::Octree::traverseNode(glm::vec3 *tmin, glm::vec3 *tmax, unsigned char directionBits, KIRK::Intersection *hit)
{
    if((tmax->x < 0.f) || (tmax->y < 0.f) || (tmax->z < 0.f))
        return false;

    if(m_hasCandidates)
    {
        float tMax = glm::compMin(*tmax);
        return Container::closestIntersectionWithCandidates(hit, 0.f, tMax);
    }

    if(!hasChildren())
        return false;

    glm::vec3 tmid = 0.5f * (*tmin + *tmax);

    unsigned char currentNode = getFirstNode(tmin, &tmid);

    do
    {
        glm::vec3 new_tmin, new_tmax;
        switch(currentNode)
        {
            case 0:
                new_tmin = glm::vec3(tmin->x, tmin->y, tmin->z);
                new_tmax = glm::vec3(tmid.x, tmid.y, tmid.z);
                if(m_children[directionBits]->traverseNode(&new_tmin, &new_tmax, directionBits, hit))
                    return true;
                currentNode = getNextNode(&new_tmax, 4, 2, 1);
                break;
            case 1:
                new_tmin = glm::vec3(tmin->x, tmin->y, tmid.z);
                new_tmax = glm::vec3(tmid.x, tmid.y, tmax->z);
                if(m_children[1 ^ directionBits]->traverseNode(&new_tmin, &new_tmax, directionBits, hit))
                    return true;
                currentNode = getNextNode(&new_tmax, 5, 3, 8);
                break;
            case 2:
                new_tmin = glm::vec3(tmin->x, tmid.y, tmin->z);
                new_tmax = glm::vec3(tmid.x, tmax->y, tmid.z);
                if(m_children[2 ^ directionBits]->traverseNode(&new_tmin, &new_tmax, directionBits, hit))
                    return true;
                currentNode = getNextNode(&new_tmax, 6, 8, 3);
                break;
            case 3:
                new_tmin = glm::vec3(tmin->x, tmid.y, tmid.z);
                new_tmax = glm::vec3(tmid.x, tmax->y, tmax->z);
                if(m_children[3 ^ directionBits]->traverseNode(&new_tmin, &new_tmax, directionBits, hit))
                    return true;
                currentNode = getNextNode(&new_tmax, 7, 8, 8);
                break;
            case 4:
                new_tmin = glm::vec3(tmid.x, tmin->y, tmin->z);
                new_tmax = glm::vec3(tmax->x, tmid.y, tmid.z);
                if(m_children[4 ^ directionBits]->traverseNode(&new_tmin, &new_tmax, directionBits, hit))
                    return true;
                currentNode = getNextNode(&new_tmax, 8, 6, 5);
                break;
            case 5:
                new_tmin = glm::vec3(tmid.x, tmin->y, tmid.z);
                new_tmax = glm::vec3(tmax->x, tmid.y, tmax->z);
                if(m_children[5 ^ directionBits]->traverseNode(&new_tmin, &new_tmax, directionBits, hit))
                    return true;
                currentNode = getNextNode(&new_tmax, 8, 7, 8);
                break;
            case 6:
                new_tmin = glm::vec3(tmid.x, tmid.y, tmin->z);
                new_tmax = glm::vec3(tmax->x, tmax->y, tmid.z);
                if(m_children[6 ^ directionBits]->traverseNode(&new_tmin, &new_tmax, directionBits, hit))
                    return true;
                currentNode = getNextNode(&new_tmax, 8, 8, 7);
                break;
            case 7:
                new_tmin = glm::vec3(tmid.x, tmid.y, tmid.z);
                new_tmax = glm::vec3(tmax->x, tmax->y, tmax->z);
                if(m_children[7 ^ directionBits]->traverseNode(&new_tmin, &new_tmax, directionBits, hit))
                    return true;
                currentNode = 8;
                break;
        }

    } while(currentNode < 8);

    return false;
}


bool KIRK::CPU::Octree::traverseNode(glm::vec3 *tmin, glm::vec3 *tmax, unsigned char directionBits, float tMax, KIRK::Ray *ray)
{
    if((tmax->x < 0.f) || (tmax->y < 0.f) || (tmax->z < 0.f) || (tmin->x > tMax) || (tmin->y > tMax) ||
       (tmin->z > tMax))
        return false;

    if(m_hasCandidates)
    {
        float tMax = glm::compMin(*tmax);
        if(Container::isIntersectionWithCandidates(ray, tMax))
            return true;
    }

    if(!hasChildren())
        return false;

    glm::vec3 tmid = 0.5f * (*tmin + *tmax);

    unsigned char currentNode = getFirstNode(tmin, &tmid);

    do
    {
        glm::vec3 new_tmin, new_tmax;
        switch(currentNode)
        {
            case 0:
                new_tmin = glm::vec3(tmin->x, tmin->y, tmin->z);
                new_tmax = glm::vec3(tmid.x, tmid.y, tmid.z);
                if(m_children[directionBits]->traverseNode(&new_tmin, &new_tmax, directionBits, tMax, ray))
                    return true;
                currentNode = getNextNode(&new_tmax, 4, 2, 1);
                break;
            case 1:
                new_tmin = glm::vec3(tmin->x, tmin->y, tmid.z);
                new_tmax = glm::vec3(tmid.x, tmid.y, tmax->z);
                if(m_children[1 ^ directionBits]->traverseNode(&new_tmin, &new_tmax, directionBits, tMax, ray))
                    return true;
                currentNode = getNextNode(&new_tmax, 5, 3, 8);
                break;
            case 2:
                new_tmin = glm::vec3(tmin->x, tmid.y, tmin->z);
                new_tmax = glm::vec3(tmid.x, tmax->y, tmid.z);
                if(m_children[2 ^ directionBits]->traverseNode(&new_tmin, &new_tmax, directionBits, tMax, ray))
                    return true;
                currentNode = getNextNode(&new_tmax, 6, 8, 3);
                break;
            case 3:
                new_tmin = glm::vec3(tmin->x, tmid.y, tmid.z);
                new_tmax = glm::vec3(tmid.x, tmax->y, tmax->z);
                if(m_children[3 ^ directionBits]->traverseNode(&new_tmin, &new_tmax, directionBits, tMax, ray))
                    return true;
                currentNode = getNextNode(&new_tmax, 7, 8, 8);
                break;
            case 4:
                new_tmin = glm::vec3(tmid.x, tmin->y, tmin->z);
                new_tmax = glm::vec3(tmax->x, tmid.y, tmid.z);
                if(m_children[4 ^ directionBits]->traverseNode(&new_tmin, &new_tmax, directionBits, tMax, ray))
                    return true;
                currentNode = getNextNode(&new_tmax, 8, 6, 5);
                break;
            case 5:
                new_tmin = glm::vec3(tmid.x, tmin->y, tmid.z);
                new_tmax = glm::vec3(tmax->x, tmid.y, tmax->z);
                if(m_children[5 ^ directionBits]->traverseNode(&new_tmin, &new_tmax, directionBits, tMax, ray))
                    return true;
                currentNode = getNextNode(&new_tmax, 8, 7, 8);
                break;
            case 6:
                new_tmin = glm::vec3(tmid.x, tmid.y, tmin->z);
                new_tmax = glm::vec3(tmax->x, tmax->y, tmid.z);
                if(m_children[6 ^ directionBits]->traverseNode(&new_tmin, &new_tmax, directionBits, tMax, ray))
                    return true;
                currentNode = getNextNode(&new_tmax, 8, 8, 7);
                break;
            case 7:
                new_tmin = glm::vec3(tmid.x, tmid.y, tmid.z);
                new_tmax = glm::vec3(tmax->x, tmax->y, tmax->z);
                if(m_children[7 ^ directionBits]->traverseNode(&new_tmin, &new_tmax, directionBits, tMax, ray))
                    return true;
                currentNode = 8;
                break;
        }

    } while(currentNode < 8);

    return false;
}


unsigned char KIRK::CPU::Octree::getFirstNode(glm::vec3 *tmin, glm::vec3 *tmid)
{
    unsigned char firstNode = 0;

    if(tmin->x > tmin->y)
        if(tmin->x > tmin->z)
        {
            if(tmid->y < tmin->x)
                firstNode |= 2;
            if(tmid->z < tmin->x)
                firstNode |= 1;
        } else
        {
            if(tmid->x < tmin->z)
                firstNode |= 4;
            if(tmid->y < tmin->z)
                firstNode |= 2;
        }
    else if(tmin->y > tmin->z)
    {
        if(tmid->x < tmin->y)
            firstNode |= 4;
        if(tmid->z < tmin->y)
            firstNode |= 1;
    } else
    {
        if(tmid->x < tmin->z)
            firstNode |= 4;
        if(tmid->y < tmin->z)
            firstNode |= 2;
    }
    return firstNode;
}


unsigned char
KIRK::CPU::Octree::getNextNode(glm::vec3 *tmax, unsigned char nodeX, unsigned char nodeY, unsigned char nodeZ)
{
    if(tmax->x < tmax->y)
        if(tmax->x < tmax->z)
            return nodeX;
        else
            return nodeZ;
    else if(tmax->y < tmax->z)
        return nodeY;
    else
        return nodeZ;

    return 0;
}
