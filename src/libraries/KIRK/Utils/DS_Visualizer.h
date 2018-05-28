#ifndef __CVK_DS_VISUALIZER_H__
#define __CVK_DS_VISUALIZER_H__

#include "CVK/CVK_2/CVK_ShaderMinimal.h"
#include "CVK/CVK_2/CVK_FBO.h"
#include "CVK/CVK_Utils/CVK_ShaderMixTextures.h"
#include "../CPU/CPU_Datastructures/Tree.h"
#include "../CPU/CPU_Datastructures/TreeAccel.h"

namespace KIRK {

namespace CPU{
    class BVHNode;

    class KDNode;

    class Octree;

    class UniformGrid;
}

/*** WARNING: THIS MIGHT NOT WORK. The Visualizer was not yet modified to work with our new framework version. Use at own risk. */
// PS: if you fix this please create a merge request on gitlab
/*
 * @brief Abstract class for visualizing raytracing acceleration datastructures
 */
class Visualizer
{
public:

	/*
	* @brief Visualizer class destructor; free resources on the GPU
	*/
    virtual ~Visualizer();

	/*
	* @brief Render the visualization and mix the result with a given rendered image
	* @param rendered_image usually a result of raytracing a scene
	*/
    void blendBoxes(unsigned int rendered_image);

	/*
	* @brief Resize the fbo so that it fits the new viewport
	* @param width new window width
	* @param height new window height
	*/
    void resize(int width, int height);
    virtual void setMaxLevel(int &max_level) = 0;
    virtual void render() = 0;
protected:

	/*
	* @brief Visualizer class constructor; allocate resources on the GPU, initialize shaders
	* @param width current window width
	* @param height current window height
	*/
    Visualizer(int width, int height);

    CVK::FBO m_fbo;
    CVK::ShaderMinimal m_draw_shader;
    CVK::ShaderMixTextures m_mix_shader;
    GLuint m_vao;
    GLuint m_vbo;
    static const char *cDrawShaders[2];
    static const char *cMixShaders[2];
};

/*
 * @brief Basic functionality for tree visualization
 */
template<typename T>
class TreeVisualizerBase : public Visualizer
{
public:

	/*
	* @brief Tree visualizer base class destructor
	*/
    virtual ~TreeVisualizerBase();

	/*
	* @brief Render the tree visualization with OpenGL
	*/
    virtual void render() override final;
protected:

	/*
	* @brief Tree visualizer base class constructor; call base class constructor, build buffer
	* @param width current window width
	* @param height current window height
	* @param root tree to visualize and traverse
	*/
    TreeVisualizerBase(int width, int height, const T &root);

	/*
	* @brief Expand vertex buffer (well...for now, it's just an array) to include another bounding box
	* @param bmin mininum bounds of node's bounding box
	* @param bmax maximum bounds of node's bounding box
	*/
    void addToBuffer(const glm::vec3 &bmin, const glm::vec3 &bmax);

    std::vector<glm::vec4> m_aabb;
    const T &m_root;
    int m_max_level;
};

/*
 * @brief Visualization of general trees, e.g. a Uniform Grid (you don't believe that's a tree? Think again!)
 */
template<typename T>
class TreeVisualizer : public TreeVisualizerBase<T>
{
public:

	/*
	* @brief Class constructor for actual visualizer; works for Uniform Grid
	* @param width current window width
	* @param height current window height
	* @param root tree to visualize and traverse
	*/
    TreeVisualizer(int width, int height, const T &root);
    virtual void setMaxLevel(int &max_level) override final;
private:
    void traverseNode(const T &node, int level, const glm::vec3 &bmin, const glm::vec3 &bmax);
};

/*
 * @brief Visualization of e.g. KD-Tree, BVH
 */
template<typename T>
    class TreeVisualizer<CPU::TreeAccelBase<T>> : public TreeVisualizerBase<CPU::TreeAccelBase<T>>
{
public:

	/*
	* @brief Class constructor for actual visualizer; works for KD-Tree, BVH and such
	* @param width current window width
	* @param height current window height
	* @param root tree to visualize and traverse
	*/
    TreeVisualizer(int width, int height, const CPU::TreeAccelBase<T> &root);
    virtual void setMaxLevel(int &max_level) override final;
private:
    void traverseNode(const T &node, int level, const glm::vec3 &bmin, const glm::vec3 &bmax);
};

/*
 * @brief Visualization of e.g. N-Tree, Octree (in fact, any type of KIRK::Tree)
 */
template<typename T>
    class TreeVisualizer<CPU::Tree<T>> : public TreeVisualizerBase<T>
{
public:

	/*
	* @brief Class constructor for actual visualizer; works for Octree, N-Tree and such
	* @param width current window width
	* @param height current window height
	* @param root tree to visualize and traverse
	*/
    TreeVisualizer(int width, int height, const T &root);
    virtual void setMaxLevel(int &max_level) override final;
private:
    void traverseNode(const T &node, int level, const glm::vec3 &bmin, const glm::vec3 &bmax);
};

    typedef TreeVisualizer<CPU::TreeAccelBase<CPU::KDNode>> KDVisualizer;
typedef TreeVisualizer<CPU::TreeAccelBase<CPU::BVHNode>> BVHVisualizer;
    typedef TreeVisualizer<CPU::Tree<CPU::Octree>> OctreeVisualizer;
typedef TreeVisualizer<CPU::UniformGrid> GridVisualizer;
}
#include "../CPU/CPU_Datastructures/CPU_KD.h"

#endif
