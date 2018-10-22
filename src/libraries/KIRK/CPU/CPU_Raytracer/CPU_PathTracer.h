#ifndef __CVK_PT_PATHTRACER_H
#define __CVK_PT_PATHTRACER_H

#include <glm/glm.hpp>

#include "CPU_Raytracer.h"
#include "KIRK/Utils/Threading.h"
#include "KIRK/Utils/Gui/Gui.h"

#include <memory>
#include "KIRK/CPU/CPU_Datastructures/CPU_DataStructure.h"
#include "KIRK/Utils/Clock.h"
#include "KIRK/Utils/BufferSegmentation.h"
#include "KIRK/Utils/Sampler.h"
#include "KIRK/Utils/Tonemapping.h"

namespace KIRK {
namespace CPU {

/**
* \brief Holds the color and the light contribution of a Bounce
*/
struct Bounce
{
	Color::RGBA color;
	glm::vec3 radiance;
	int mat_flags;
};

/**
 * \brief PathTracer class
 */
class PathTracer : public CPU_Raytracer, public virtual KIRK::GuiElement
{
public:

	PathTracer() : KIRK::GuiElement(),
	               m_tonemapper(std::make_unique<Tonemapper>()),
	               m_samples_per_pixel(1),
	               m_inverse_samples_per_pixel(1),
				   m_maxBufferSize(1e+8){};

    PathTracer(std::shared_ptr<KIRK::CPU::Scene> cpuscene, int depth);

    virtual ~PathTracer() { }

	/**
	 * \brief Sets the number of samples per pixel
	 * \param samples number of samples per pixel
	 */
	void setSampleCount(const unsigned int samples);
	/**
	* \brief Get the number of the currently finished samples
	* \return the current number of finished samples
	*/
	int getCurrentSampleCount();
	/**
	* \brief Not being called anywhere. Should provide an interface to notify for window resizes.
	*/
	void windowResizeCallback();

	/**
	 * \brief Resets the rendering progress
	 */
	void reset();
	void onGui() override;

	UniformSampler& getSampler() { return m_sampler; }

	std::unique_ptr<Tonemapper> m_tonemapper;

protected:

	/**
	 * \brief This is the render method and provides the basic render pipeline for segmentated rendering
	 */
	void render() override;

	/**
	 * \brief This holds the actual PathTracing Routine performed everytime the a segment gets rendered
	 */
	virtual void processSegment();

	/**
	 * \brief Resets/Initializes all used buffers and (re)calculates the desired segmentation
	 */
	virtual void resetBuffers();

	/**
	 * \brief Tries to initiate the buffer Segmentation and catches exception if maxBufferSize isn't sufficient.
	 * \param	staticBufferObjectSize	The sum of the size (in bytes) of all objects being used in total Buffers. (e.g. a depth buffer with the total size of a renderTexture)
	 * \param	dynamicBufferObjectSize	The sum of the size (in bytes) of all objects being used in a segmented Buffer. (e.g. a hit Buffer)
	 */
	void initBufferSegmentation(size_t staticBufferObjectSize, size_t dynamicBufferObjectSize);

	/**
	 * \brief Writes sample weighted colors of processed segment in textureBuffer and updates renderTexture
	 */
	virtual void drawTexture();

	void applyToneMapping();

	/**
	 * \brief Sets the bounceBuffers colors to vec4(0) and radiance to vec4(1)
	 */
	virtual void clearBufferWeights();

	/**
	 * \brief Generates all primary rays and writes them to m_rayBuffer.
	 */
	virtual void generatePrimaryRays();

	/**
	 * \brief Traces all rays iteratively for depth from m_rayBuffer and shades their hits.
	*/
	virtual void traceRays();

	/**
	 * \brief Traces a single ray at given buffer index, determines its hit and writes it to hitBuffer.
	 * \param segmentIndex index for accessing segmented buffer
	*/
	virtual void traceRay(int segmentIndex, bool is_primary, int mat_flags);


	// Buffers

	/**
	 * \brief Color Buffer of total image
	 */
	std::vector <glm::vec4> m_textureBuffer;
	/**
	 * \brief Segmented Buffer for Rays
	 */
	std::vector <Ray> m_rayBuffer;
	/**
	 * \brief Segmented Buffer for Bounces
	 */
	std::vector <Bounce> m_bounceBuffer;
	/**
	 * \brief Segmented Buffer for Intersections
	 */
	std::vector <Intersection> m_hitBuffer;

	UniformSampler m_sampler;
	/**
	 * \brief Protection for threaded execution of setPixel in the KIRK::Texture
	 */
	std::mutex m_texture_mutex;
	unsigned int m_samples_per_pixel;
	float m_inverse_samples_per_pixel;
	int c_sample = 0;
	bool m_use_tonemapping = false;
	bool m_has_resized = false; //!< If window was resized

	/**
	 * \brief Utility class for buffer segmentation
	 */
	BufferSegmentation m_bufferSegmentation;
	/**
	 * \brief The maximum size in Bytes all Buffers are allowed to allocate.
	 */
	size_t m_maxBufferSize;

	int m_num_threads = KIRK::ThreadManager::maxThreads();
	bool m_use_progress_bar = false;
};
}
}

#endif
