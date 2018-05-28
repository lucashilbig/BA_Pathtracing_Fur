#include "DS_Visualizer.h"
#include "../CPU/CPU_Datastructures/CPU_BVH.h"
#include "../CPU/CPU_Datastructures/CPU_KD.h"
#include "../CPU/CPU_Datastructures/Octree.h"
#include "../CPU/CPU_Datastructures/UniformGrid.h"

namespace KIRK {

const char *Visualizer::cDrawShaders[2] = {SHADERS_PATH "/OCL_ParticleSystem/ocl_ps.vert",
                                           SHADERS_PATH "/OCL_ParticleSystem/ocl_ps.frag"};
const char *Visualizer::cMixShaders[2] = {SHADERS_PATH "/Examples/screenFill.vert",
                                          SHADERS_PATH "/Mix/mix_textures.frag"};

Visualizer::Visualizer(int width, int height) :
        m_fbo(CVK::FBO(width, height, 1, false)),
        m_draw_shader(CVK::ShaderMinimal(VERTEX_SHADER_BIT | FRAGMENT_SHADER_BIT, cDrawShaders)),
    m_mix_shader(CVK::ShaderMixTextures(VERTEX_SHADER_BIT | FRAGMENT_SHADER_BIT, cMixShaders))
{

    glGenVertexArrays(1, &m_vao);
    glGenBuffers(1, &m_vbo);
}

Visualizer::~Visualizer()
{
    glDeleteBuffers(1, &m_vbo);
    glDeleteVertexArrays(1, &m_vao);
}

void Visualizer::blendBoxes(unsigned int rendered_image)
{
    m_fbo.bind();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    render();
    m_fbo.unbind();
    m_mix_shader.useProgram();
    m_mix_shader.setTextureInput(0, rendered_image);
    m_mix_shader.setTextureInput(1, m_fbo.getColorTexture(0));
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    m_mix_shader.update();
    m_mix_shader.render();
}

void Visualizer::resize(int width, int height)
{
    m_fbo.resize(width, height);
}

template<typename T>
TreeVisualizerBase<T>::TreeVisualizerBase(int width, int height, const T &root) :
        Visualizer(width, height),
        m_root(root)
{
}

template<typename T>
TreeVisualizerBase<T>::~TreeVisualizerBase()
{
}

template<typename T>
void TreeVisualizerBase<T>::render()
{
    CVK::State::getInstance()->setShader(&m_draw_shader); // Simple shader, mvp and red color
    m_draw_shader.update();
    glBindVertexArray(m_vao);
    glDrawArrays(GL_LINES, 0,
                 m_aabb.size()); // Each bounding box has 8 defining coordinates, render lines connecting them, forming a box
}

template<typename T>
void TreeVisualizerBase<T>::addToBuffer(const glm::vec3 &bmin, const glm::vec3 &bmax)
{
    const glm::vec4 a(bmin.x, bmin.y, bmax.z, 1.0f);
    const glm::vec4 b(bmin.x, bmax.y, bmin.z, 1.0f);
    const glm::vec4 c(bmax.x, bmin.y, bmin.z, 1.0f);
    const glm::vec4 d(bmin.x, bmax.y, bmax.z, 1.0f);
    const glm::vec4 e(bmax.x, bmin.y, bmax.z, 1.0f);
    const glm::vec4 f(bmax.x, bmax.y, bmin.z, 1.0f);

    m_aabb.push_back(a);
    m_aabb.push_back(d);
    m_aabb.push_back(a);
    m_aabb.push_back(e);
    m_aabb.push_back(b);
    m_aabb.push_back(d);
    m_aabb.push_back(b);
    m_aabb.push_back(f);
    m_aabb.push_back(c);
    m_aabb.push_back(e);
    m_aabb.push_back(c);
    m_aabb.push_back(f);
    m_aabb.push_back(glm::vec4(bmin, 1));
    m_aabb.push_back(a);
    m_aabb.push_back(glm::vec4(bmin, 1));
    m_aabb.push_back(b);
    m_aabb.push_back(glm::vec4(bmin, 1));
    m_aabb.push_back(c);
    m_aabb.push_back(glm::vec4(bmax, 1));
    m_aabb.push_back(d);
    m_aabb.push_back(glm::vec4(bmax, 1));
    m_aabb.push_back(e);
    m_aabb.push_back(glm::vec4(bmax, 1));
    m_aabb.push_back(f);
}

template<typename T>
TreeVisualizer<T>::TreeVisualizer(int width, int height, const T &root) :
        TreeVisualizerBase<T>(width, height, root)
{
}

/*
 * @brief Class constructor for actual visualizer; works for Octree, N-Tree and such
 * @param width current window width
 * @param height current window height
 * @param root tree to visualize and traverse
 */
template<typename T>
    TreeVisualizer<CPU::Tree<T>>::TreeVisualizer(int width, int height, const T &root) :
        TreeVisualizerBase<T>(width, height, root)
{
}

/*
 * @brief Class constructor for actual visualizer; works for KD-Tree, BVH and such
 * @param width current window width
 * @param height current window height
 * @param root tree to visualize and traverse
 */
template<typename T>
TreeVisualizer<CPU::TreeAccelBase<T>>::TreeVisualizer(int width, int height, const CPU::TreeAccelBase<T> &root) :
        TreeVisualizerBase<CPU::TreeAccelBase<T>>(width, height, root)
{
}

/*
 * @brief Change at which point no more levels should be rendered; this will re-build the vertex buffer. This is not used for KIRK::Tree types
 * @param max_level maximum level to render; should be modified outside of this class, too, hence the reference
 */
template<typename T>
void TreeVisualizer<CPU::TreeAccelBase<T>>::setMaxLevel(int &max_level)
{
    // numeric_limits is the default value, thus the tree depth will be used if no custom level is provided
    if(max_level == std::numeric_limits<int>::max())
    {
        max_level = this->m_root.m_prop.max_depth;
    }

    this->m_max_level = max_level = std::max(0, std::min(max_level,
                                                         this->m_root.m_prop.depth)); // Make sure this->m_max_level isn't negative or above tree depth
    this->m_aabb.clear();

    // Each node has 24 coordinates for 12 lines representing a bounding box
    this->m_aabb.reserve(this->m_root.m_prop.amnt_nodes * 24);
    traverseNode(this->m_root.m_root, 0, *((const_cast < CPU::TreeAccelBase<T> & > ( this->m_root )).getMinBound()),
                 *((const_cast < CPU::TreeAccelBase<T> & > ( this->m_root )).getMaxBound())); // will fill this->m_aabb
    this->m_aabb.shrink_to_fit(); // Since we probably allocated too much memory, free any unused resources

    glBindVertexArray(this->m_vao);
    glBindBuffer(GL_ARRAY_BUFFER, this->m_vbo);
    glBufferData(GL_ARRAY_BUFFER, this->m_aabb.size() * sizeof(glm::vec4), &this->m_aabb[0], GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

// Forward declaration; will be defined later
template<>
void TreeVisualizer<CPU::UniformGrid>::traverseNode(const CPU::UniformGrid &node, int level, const glm::vec3 &bmin,
                                               const glm::vec3 &bmax);

/*
 * @brief This is relevant only when building the grid visualization
 * @param max_level amount of voxels on each axis
 */
template<>
void TreeVisualizer<CPU::UniformGrid>::setMaxLevel(int &max_level)
{
    if(!m_aabb.empty())
    {
        return; // exit if grid has been built already
    }

    // In this case, m_max_level actually denotes the number of voxels in this grid (abusing variables, yay!)
    m_max_level = (max_level == std::numeric_limits<int>::max() || max_level < 0) ? 0 : max_level;

    // Each voxel has 24 coordinates for 12 lines representing it
    // Each dimension spans max_level voxels
    this->m_aabb.reserve(max_level * max_level * max_level * 24);
    traverseNode(this->m_root, 0, *((const_cast < CPU::UniformGrid & > ( this->m_root )).getMinBound()),
                 *((const_cast < CPU::UniformGrid & > ( this->m_root )).getMaxBound())); // will fill this->m_aabb
    this->m_aabb.shrink_to_fit();

    glBindVertexArray(this->m_vao);
    glBindBuffer(GL_ARRAY_BUFFER, this->m_vbo);
    glBufferData(GL_ARRAY_BUFFER, this->m_aabb.size() * sizeof(glm::vec4), &this->m_aabb[0], GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

/*
 * @brief Change at which point no more levels should be rendered; this will re-build the vertex buffer. This is only used for KIRK::Tree types
 * @param max_level maximum level to render; should be modified outside of this class, too, hence the reference
 */
template<typename T>
    void TreeVisualizer<CPU::Tree<T>>::setMaxLevel(int &max_level)
{
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // DISCLAIMER; NOTES ON PROGRAMMING STYLE; MIGHT WANT TO READ THIS
    // Unfortunately, can't make sure that tree depth is exceeded by the given argument. The Octree/N-Tree code does not have a get-method for its internal depth representation.
    // Keep in mind this is actually a good thing; giving away information about internal data of your class is a design flaw. Classes should be defined by their behaviour, not their data.
    // (Consider using structs if you wish to reveal much data to outside code)
    // Normally, you'd integrate the visualization in the datastructure class to prevent this sort of issue.
    // In this case, it is not wanted that the original Tree datastructure we were provided with is changed in any way.
    // It doesn't matter too much after all, since we'll add to the buffer only if the current node has children.
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    this->m_max_level = max_level = std::max(0, max_level); // Can still make sure this->m_max_level isn't negative
    this->m_aabb.clear();

    // Each node has 24 coordinates for 12 lines representing a bounding box
    this->m_aabb.reserve((const_cast < T & > ( this->m_root )).numTotalChildren() * 24);
    const glm::vec3 &bmin = *((const_cast < T & > ( this->m_root )).getMinBound());
    const glm::vec3 &bmax = *((const_cast < T & > ( this->m_root )).getMaxBound());
    traverseNode(this->m_root, 0, bmin, bmax); // will fill this->m_aabb
    this->m_aabb.shrink_to_fit(); // Since we probably allocated too much memory, free any unused resources

    glBindVertexArray(this->m_vao);
    glBindBuffer(GL_ARRAY_BUFFER, this->m_vbo);
    glBufferData(GL_ARRAY_BUFFER, this->m_aabb.size() * sizeof(glm::vec4), &this->m_aabb[0], GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

template<>
void TreeVisualizer<CPU::TreeAccelBase<CPU::KDNode>>::traverseNode(const CPU::KDNode &node, int level, const glm::vec3 &bmin,
                                                         const glm::vec3 &bmax)
{
    // Render only those boxes that are furthest down the tree.
    // If a depth is given that does not exceed the tree depth, that must include any leaf nodes within this range.
    // All parent's boxes are made up of these "leaf" boxes recursively! No need to render all these lines multiple times.
    if(level++ == this->m_max_level || node.m_axis == -1)
    {
        this->addToBuffer(bmin, bmax);

        return;
    }

    glm::vec3 left_clamp = bmax;
    glm::vec3 right_clamp = bmin;
    left_clamp[node.m_axis] = node.m_split_pos;
    right_clamp[node.m_axis] = node.m_split_pos;

    // If this node has no children, axis will be -1 and we will not enter this code block, hence these child pointers are valid
    traverseNode(*(node.m_children[0]), level, bmin, left_clamp);
    traverseNode(*(node.m_children[1]), level, right_clamp, bmax);
}

/*
 * @brief Traverse a BVH node in order to retrieve its bounding volume coordinates
 * @param node current BVH node
 * @param level current recursion depth (level in tree)
 * @param bmin mininum bounds of node's bounding volume
 * @param bmax maximum bounds of node's bounding volume
 */
template<>
    void TreeVisualizer<CPU::TreeAccelBase<CPU::BVHNode>>::traverseNode(const CPU::BVHNode &node, int level, const glm::vec3 &bmin,
                                                          const glm::vec3 &bmax)
{
    // Reached maximum level, stop traversing.
    if(level++ > this->m_max_level)
    {
        return;
    }

    this->addToBuffer(bmin, bmax);

    // If left child exists, so does the right child - continue recursion.
    if(node.LEFT_CHILD)
    {
        traverseNode(*(node.LEFT_CHILD), level, node.LEFT_CHILD->m_bvol.min(), node.LEFT_CHILD->m_bvol.max());
        traverseNode(*(node.RIGHT_CHILD), level, node.RIGHT_CHILD->m_bvol.min(), node.RIGHT_CHILD->m_bvol.max());
    }
}

/*
 * @brief Doesn't make much sense. We abuse this function so we don't have to make a new class. Basically builds the grid visualization.
 * @param node the grid
 * @param level irrelevant
 * @param bmin minimum scene bounds
 * @param bmax maximum scene bounds
 */
template<>
    void TreeVisualizer<CPU::UniformGrid>::traverseNode(const CPU::UniformGrid &node, int level, const glm::vec3 &bmin,
                                               const glm::vec3 &bmax)
{
    if(m_max_level == 0)
    {
        return; // Please no division by zero
    }

    const glm::vec3 size = bmax - bmin;
    const glm::vec3 voxel_size = size / glm::vec3(m_max_level); // use max level as amount of voxels

    // Go over each dimension, increasing the step by voxel size. This will collect all voxels spanning the scene
    for(float x = bmin.x; x < bmax.x; x += voxel_size.x)
    {
        for(float y = bmin.y; y < bmax.y; y += voxel_size.y)
        {
            for(float z = bmin.z; z < bmax.z; z += voxel_size.z)
            {
                glm::vec3 tmin(x, y, z);
                glm::vec3 tmax = tmin +
                                 voxel_size; // Each voxel we collect has a corresponding max bound that we need in order to set the line coordinates
                addToBuffer(tmin, tmax); // Adds 12 vertices
            }
        }
    }
}

/*
 * @brief Traverse a KIRK::Tree node in order to retrieve its bounding coordinates
 * @param node current KIRK::Tree node
 * @param level current recursion depth (level in tree)
 * @param bmin mininum bounds of node's bounding box
 * @param bmax maximum bounds of node's bounding box
 */
template<typename T>
    void TreeVisualizer<CPU::Tree<T>>::traverseNode(const T &node, int level, const glm::vec3 &bmin, const glm::vec3 &bmax)
{
    // Render only those boxes that are furthest down the tree.
    // If a depth is given that does not exceed the tree depth, that must include any leaf nodes within this range.
    // All parent's boxes are made up of these "leaf" boxes recursively! No need to render all these lines multiple times.
    // There is one exception to this; any node further up the tree that does not have
    // children would not be drawn as it's not on the max level. This provides us with the second (optional) condition.
    if(level++ == this->m_max_level || !((const_cast < T & > ( node )).hasChildren()))
    {
        this->addToBuffer(bmin, bmax);

        return;
    }

    auto &children = *((const_cast < T & > ( node )).getChildren());

    for(auto &child : children)
    {
        traverseNode(*child, level, *(child->getMinBound()), *(child->getMaxBound()));
    }
}

template
    class TreeVisualizer<CPU::TreeAccelBase<CPU::KDNode> >;

template
    class TreeVisualizer<CPU::TreeAccelBase<CPU::BVHNode>>;

template
class TreeVisualizer<CPU::Tree<CPU::Octree>>;

template
class TreeVisualizer<CPU::UniformGrid>;

template
class TreeVisualizerBase<CPU::TreeAccelBase<CPU::KDNode>>;

template
class TreeVisualizerBase<CPU::TreeAccelBase<CPU::BVHNode>>;

template
class TreeVisualizerBase<CPU::Octree>;

template
class TreeVisualizerBase<CPU::UniformGrid>;

}
