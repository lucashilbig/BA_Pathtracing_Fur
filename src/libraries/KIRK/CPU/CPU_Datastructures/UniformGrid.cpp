#include "KIRK/Utils/Log.h"
#include "UniformGrid.h"

//KIRK::UniformGrid::UniformGrid(KIRK::Scene *scene, float sizeOfVoxel)
//{
//	m_voxelSize = glm::vec3(sizeOfVoxel, sizeOfVoxel, sizeOfVoxel);
//	m_numVoxels = glm::u32vec3( glm::ceil( m_size / m_voxelSize ) );
//
//	init();
//	addBaseDataStructure(scene);
//}
//
//KIRK::UniformGrid::UniformGrid(KIRK::Scene *scene, glm::vec3 sizeOfVoxel)
//{
//	m_voxelSize = sizeOfVoxel;
//	m_numVoxels = glm::u32vec3(glm::ceil(m_size / m_voxelSize));
//
//	init();
//	addBaseDataStructure(scene);
//}

KIRK::CPU::UniformGrid::UniformGrid(glm::vec3 minBound, glm::vec3 maxBound, int numVoxelsPerAxis)
        : CPU_DataStructure(minBound, maxBound)
{
    m_numVoxels = glm::u32vec3(numVoxelsPerAxis, numVoxelsPerAxis, numVoxelsPerAxis);
    m_voxelSize = m_size / glm::vec3(m_numVoxels);

    init();
}

KIRK::CPU::UniformGrid::~UniformGrid()
{

}


void KIRK::CPU::UniformGrid::printDebugInfo()
{
	LOG_DEBUG("Uniform Grid of Size: % x % x % ", m_numVoxels.x, m_numVoxels.y, m_numVoxels.z);
	LOG_DEBUG("Memory size: % bytes.", getSizeInBytes());
	LOG_DEBUG("# of Nodes: % ", m_totalVoxels);
}

int KIRK::CPU::UniformGrid::getSizeInBytes()
{
    int result = sizeof(*this);
    result += m_voxelCandidates.size() * sizeof(std::vector<KIRK::Triangle *>);
    for(std::vector<KIRK::Triangle *> v : m_voxelCandidates)
        result += v.size() * sizeof(KIRK::Triangle *);
    return result;
}

void KIRK::CPU::UniformGrid::init()
{
    m_invVoxelSize = glm::vec3(1 / m_voxelSize.x, 1 / m_voxelSize.y, 1 / m_voxelSize.z);
    m_totalVoxels = m_numVoxels.x * m_numVoxels.y * m_numVoxels.z;
    m_maxIndex = m_numVoxels - glm::ivec3(1, 1, 1);

    m_voxelCandidates = std::vector<std::vector<KIRK::Triangle *>>(m_totalVoxels);
}


void KIRK::CPU::UniformGrid::addBaseDataStructure(KIRK::CPU::Scene *scene)
{
    const std::vector<KIRK::Triangle *> obj_list = scene->getSceneObjects();
    setBounds(scene->getBounds());

    for(unsigned int index = 0; index < obj_list.size(); index++)
    {
		KIRK::Triangle *obj = obj_list[index];
        glm::vec3 *obj_bbox = obj->getBounds();
        glm::vec3 voxel_bbox[2];
        glm::vec3 minIndex = (obj_bbox[0] - m_minBound) * m_invVoxelSize;
        glm::vec3 maxIndex = (obj_bbox[1] - m_minBound) * m_invVoxelSize;

        for(int iz = (int)minIndex.z; iz <= std::min((int)maxIndex.z, m_maxIndex.z); iz++)
            for(int iy = (int)minIndex.y; iy <= std::min((int)maxIndex.y, m_maxIndex.y); iy++)
                for(int ix = (int)minIndex.x; ix <= std::min((int)maxIndex.x, m_maxIndex.x); ix++)
                {
                    voxel_bbox[0] = m_minBound + glm::vec3(ix, iy, iz) * m_voxelSize;
                    voxel_bbox[1] = voxel_bbox[0] + m_voxelSize;
                    if(obj->isInAABB(voxel_bbox))
                    {
                        m_voxelCandidates[getIndex(ix, iy, iz)].push_back(obj);
                    }
                }
    }
}


#define d_handleNextNode(param) \
    index. param += int( step. param ); \
    if ((index. param < 0) || (index. param > m_maxIndex. param )) \
        return false; \
    tnext. param += delta. param ;

#define d_handleNode(param) \
    if (!candidates->empty()) \
        if (testClosestIntersectionWithCandidates(candidates, hit, tMin, tnext. param)) \
            return true; \
    tMin = tnext. param; \
    d_handleNextNode(param)

bool KIRK::CPU::UniformGrid::closestIntersection(KIRK::Intersection *hit)
{
    glm::vec3 direction = hit->m_ray.m_direction;
    glm::vec3 origin = hit->m_ray.m_origin;

    //glm::vec3 step((direction.x < 0) ? -1 : 1, (direction.y < 0) ? -1 : 1, (direction.z < 0) ? -1 : 1);
    glm::vec3 step = glm::sign(direction);
    glm::vec3 delta = step / direction * m_voxelSize;
    glm::vec3 posVec = (origin - m_minBound) * m_invVoxelSize;

    glm::ivec3 index;
    index.x = glm::max(glm::min((int)(posVec.x), m_maxIndex.x), 0);
    index.y = glm::max(glm::min((int)(posVec.y), m_maxIndex.y), 0);
    index.z = glm::max(glm::min((int)(posVec.z), m_maxIndex.z), 0);

    glm::vec3 tnext;
    if(step.x < 0)
        tnext.x = (posVec.x - index.x) * delta.x;
    else
        tnext.x = (1.f - posVec.x + index.x) * delta.x;
    if(step.y < 0)
        tnext.y = (posVec.y - index.y) * delta.y;
    else
        tnext.y = (1.f - posVec.y + index.y) * delta.y;
    if(step.z < 0)
        tnext.z = (posVec.z - index.z) * delta.z;
    else
        tnext.z = (1.f - posVec.z + index.z) * delta.z;

    float tMin = 0.f;
    std::vector<KIRK::Triangle *> *candidates;

    for(;;)
    {
        candidates = &m_voxelCandidates[getIndex(index)];

        if(tnext.x < tnext.y)
            if(tnext.x < tnext.z)
            {
                d_handleNode(x)
            } else
            {
                d_handleNode(z)
            }
        else if(tnext.y < tnext.z)
        {
            d_handleNode(y)
        } else
        {
            d_handleNode(z)
        }
    }
    return false;
}

bool KIRK::CPU::UniformGrid::isIntersection(KIRK::Ray *ray, float tMax)
{
    glm::vec3 direction = ray->m_direction;
    glm::vec3 origin = ray->m_origin;

    //glm::vec3 step((direction.x < 0) ? -1 : 1, (direction.y < 0) ? -1 : 1, (direction.z < 0) ? -1 : 1);
    glm::vec3 step = glm::sign(direction);
    glm::vec3 delta = step / direction * m_voxelSize;
    glm::vec3 posVec = (origin - m_minBound) * m_invVoxelSize;

    glm::ivec3 index;
    index.x = glm::max(glm::min((int)(posVec.x), m_maxIndex.x), 0);
    index.y = glm::max(glm::min((int)(posVec.y), m_maxIndex.y), 0);
    index.z = glm::max(glm::min((int)(posVec.z), m_maxIndex.z), 0);

    glm::vec3 tnext;
    if(step.x < 0)
        tnext.x = (posVec.x - index.x) * delta.x;
    else
        tnext.x = (1.f - posVec.x + index.x) * delta.x;
    if(step.y < 0)
        tnext.y = (posVec.y - index.y) * delta.y;
    else
        tnext.y = (1.f - posVec.y + index.y) * delta.y;
    if(step.z < 0)
        tnext.z = (posVec.z - index.z) * delta.z;
    else
        tnext.z = (1.f - posVec.z + index.z) * delta.z;

    std::vector<KIRK::Triangle *> *candidates;

    for(;;)
    {
        candidates = &m_voxelCandidates[getIndex(index)];

        if(!candidates->empty())
            if(testIsIntersectionWithCandidates(candidates, ray, tMax))
                return true;

        if(tnext.x < tnext.y)
            if(tnext.x < tnext.z)
            {
                d_handleNextNode(x)
            } else
            {
                d_handleNextNode(z)
            }
        else if(tnext.y < tnext.z)
        {
            d_handleNextNode(y)
        } else
        {
            d_handleNextNode(z)
        }
    }
    return false;
}


