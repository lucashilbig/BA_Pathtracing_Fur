#pragma once

#include "KIRK/Utils/Log.h"

namespace KIRK {
/**
*	The BufferSegmentation is a utility class to help segmenting 2D Buffers (e.g. a 2D Image) into smaller segments.
*	The main reason to do this is to limit the used memory of a Buffer by working in smaller portions.
*
*	@remarks The term 'Buffer' stands for the whole sized Buffer (e.g. the whole 2D Image). 
*	@remarks The term 'Segment' stands for a rectangular subset of the 'Buffer', which is not necessarily a strict subset.
*	@remarks The term 'total' indicates something refers to a Buffer and its whole size, not a Segment.
*	@remarks The term 'segmented' indicates something refers to a Segment and its segmented size, not a Buffer.
*/
class BufferSegmentation
{
public:
	/**
	*	The constructor is empty because the BufferSegmentation will most likely not happen at the time the constructor is called.
	*	@see BufferSegmentation#init
	*/
	BufferSegmentation() {};

	/**
	*	Performs initial calculations to determine the optimal segmentation.
	*	Whenever the total Buffer Size changes (e.g. changing resolution of a renderTexture) this must be called to reinitiate the Segmentation.
	*
	*	@param bufferWidth	The total width of the Buffer. (e.g. the width of a renderTexture)
	*	@param bufferHeight	The total height of the Buffer.
	*	@param	maxBufferSize	The maximum size in bytes all your buffers (total and segmented) may take. (e.g. the systems free memory)
	*	@param	staticBufferObjectSize	The sum of the size (in bytes) of all objects being used in total Buffers. (e.g. a depth buffer with the total size of a renderTexture)
	*	@param	dynamicBufferObjectSize	The sum of the size (in bytes) of all objects being used in a segmented Buffer. (e.g. a hit Buffer)
	*/
	void init (float bufferWidth, float bufferHeight, size_t maxBufferSize, size_t staticBufferObjectSize, size_t dynamicBufferObjectSize)
	{
		m_bufferWidth = bufferWidth;
		m_bufferHeight = bufferHeight;
		// Multiply staticBufferObjectSize with buffer Width and Height
		staticBufferObjectSize *= m_bufferWidth * m_bufferHeight;

		// Test if maximum Buffer Size is sufficient
		if (maxBufferSize < staticBufferObjectSize) {
			throw std::runtime_error("BufferSegmentation failed. Maximum Buffer Size exceeded.");
		}

		// Calculate Segment Data
		m_segmentLength = std::fmin(m_bufferHeight, std::fmin(m_bufferWidth, std::floor(std::sqrt((maxBufferSize - staticBufferObjectSize) / dynamicBufferObjectSize))));
		m_segmentSize = std::fmin(m_bufferWidth, m_segmentLength) * std::fmin(m_bufferHeight, m_segmentLength);
		m_segmentCount = std::ceil(m_bufferWidth / m_segmentLength) * std::ceil(m_bufferHeight / m_segmentLength);
		m_segmentId = 0;
		m_progressedOverall = 0.0f;

		updateBufferSegmentation();

		LOG_INFO("Performed Buffer Segmentation: # of Segments: %, SegmentEdgeLength: %, SegmentSize: %", m_segmentCount, m_segmentLength, m_segmentSize);
	};

	/**
	*	Increments the active segment Id, updates the indexing positions and determines if the segmentation is still active or not.
	*	This should be called at the end of your render loop.
	*/
	void update ()
	{
		m_progressedLastSegment = (m_curSegmentWidth * m_curSegmentHeight) / (m_bufferWidth * m_bufferHeight);
		m_progressedOverall += m_progressedLastSegment;

		m_segmentId++;
		if (m_segmentId >= m_segmentCount)
		{
			m_isActive = false;
			m_segmentId = 0;
		}
		else { m_isActive = true; }
		updateBufferSegmentation();
	};

	/**
	*	True when the current segment is not the last segment.
	*	@return Returns true when all segments have been processed.
	*/
	bool isActive () const { return m_isActive; };

	/**
	*	Sets isActive to be false and reset segmentCount.
	*/
	void setInactive ()
	{
		m_isActive = false;
		m_segmentId = 0;
		m_segmentCount = 0;
		m_progressedOverall = 0.0f;
	}

	/**
	*	Get the calculated perfect edge length (length = width = height) of the quadratic segment
	*	The actual width or height of the current segment may vary, this is the calculated perfect length.
	*	To get the calculated length for a segmented Buffer see getSegmentSize()
	*	@see BufferSegmentation#getSegmentWidth BufferSegmentation#getSegmentHeight
	*	@return segmentLength
	*/
	unsigned int getSegmentLength () const { return m_segmentLength; };

	/**
	*	Get the actual width of the current segment.
	*	The width may be smaller than the segmentLength when an edge of the total buffer is reached.
	*	@return segmentWidth
	*/
	unsigned int getSegmentWidth () const { return m_curSegmentWidth; };

	/**
	*	Get the actual height of the current segment.
	*	The height may be smaller than the segmentLength when an edge of the total buffer is reached.
	*	@return segmentHeight
	*/
	unsigned int getSegmentHeight () const { return m_curSegmentHeight; };

	/**
	*	Get the calculated size of a segment. This is the size for a segmented buffer.
	*	This may be smaller than segmentLength * segmentLength when the segmentLength is greater than the total buffer width or height.
	*	@return segmentSize The size for a segmented buffer.
	*/
	unsigned int getSegmentSize () const { return m_segmentSize; };

	/**
	*	Get the count of segments.
	*	@return segmentCount 
	*/
	unsigned int getSegmentCount () const { return m_segmentCount; };

	/**
	*	Get the id of the currently active segment.
	*	@return segmentId
	*/
	unsigned int getSegmentId () const { return m_segmentId; };

	/**
	*	Get the size of the total buffer. This is the size for a total buffer.
	*	@return bufferSize The size for a total buffer.
	*/
	unsigned int getBufferSize () const { return m_bufferWidth * m_bufferHeight; };

	/**
	*	Get the X position of the total buffer. This can be used in loops in combination with the offset param.
	*	@param offset Use to add or subtract an offset to the position.
	*	@return bufferPosX The value gets clamped to 0 - bufferWidth
	*/
	unsigned int getBufferPosX (int offset = 0) const { return std::fmax(std::fmin(static_cast <signed>(m_bufferPosX + offset), m_bufferWidth - 1), 0); };

	/**
	*	Get the Y position of the total buffer. This can be used in loops in combination with the offset param.
	*	@param offset Use to add or subtract an offset to the position.
	*	@return bufferPosY The value gets clamped to 0 - bufferHeight
	*/
	unsigned int getBufferPosY (int offset = 0) const { return std::fmax(std::fmin(static_cast <signed>(m_bufferPosY + offset), m_bufferHeight - 1), 0); };

	/**
	*	Get the total buffer width.
	*	@return bufferWidth The total buffer width.
	*/
	unsigned int getBufferWidth () const { return m_bufferWidth; };

	/**
	*	Get the total buffer height.
	*	@return bufferWidth The total buffer height.
	*/
	unsigned int getBufferHeight () const { return m_bufferHeight; };

	/**
	*	Get the 1D index to index a total 1d buffer at the current position. This can be used in loops in combination with the offset params.
	*	@param offsetX Use to add or subtract an offset to the X position.
	*	@param offsetY Use to add or subtract an offset to the Y position.
	*	@return bufferIndex The index gets clamped to 0 - bufferSize
	*/
	unsigned int get1DBufferIndex (int offsetX, int offsetY) const { return getBufferPosY(offsetY) * m_bufferWidth + getBufferPosX(offsetX); };

	/**
	*	Get the 1D index to index a segmented 1d buffer at the current position. This can be used in loops for easier indexing.
	*	@param x The x position.
	*	@param y The y position.
	*	@return segmentIndex (y * segmentLength + x)
	*/
	unsigned int get1DSegmentIndex (unsigned int x, unsigned int y) const { return y * m_segmentLength + x; };

	/**
	 * \brief Get the progress in percentage done last segment in ratio to total sample count. Gets updated when update() gets called.
	 * \param samples The total sample count.
	 * \return progress done last segment in percentage in ratio to total sample count
	 */
	float getProgressLastSegment(int samples = 1) const {
		return m_progressedLastSegment * (1.0f / samples) ;
	}

	/**
	* \brief Get the overall progress in percentage in ratio to total sample count. Gets updated when update() gets called.
	* \param samples Multiplier to be used for passing the current sample / max samples count.
	* \return Overall progress done in percentage in ratio to total sample count.
	*/
	float getProgressOverall(int samples = 1) const {
		return m_progressedOverall * (1.0f / samples);
	}

private:
	/**
	*	Updates current bufferPos and current segment's width and height
	*/
	void updateBufferSegmentation ()
	{
		m_bufferPosX = std::fmod(m_segmentId, std::ceil(m_bufferWidth / m_segmentLength)) * m_segmentLength;
		m_bufferPosY = std::floor(m_segmentId / std::ceil(m_bufferWidth / m_segmentLength)) * m_segmentLength;
		m_curSegmentWidth = std::fmin(m_segmentLength, m_bufferWidth - m_bufferPosX);
		m_curSegmentHeight = std::fmin(m_segmentLength, m_bufferHeight - m_bufferPosY);
		LOG_DEBUG("Next Segment #% (X: %, Y: %, W: %, H: %",
			m_segmentId, m_bufferPosX, m_bufferPosY, m_curSegmentWidth, m_curSegmentHeight);
	}

	bool m_isActive = false; //!< True when all segments have been processed.
	unsigned int m_segmentLength; //!< The edge length of the quadratic segment
	unsigned int m_segmentSize; //!< The segments surface size
	unsigned int m_curSegmentWidth; //!< The current segment width, can be smaller than segmentLength if we are close to the buffers edge
	unsigned int m_curSegmentHeight; //!< The current segment height
	unsigned int m_segmentId; //!< The id of the current segment
	unsigned int m_segmentCount; //!< The total number of segments
	unsigned int m_bufferPosX; //!< The calculated x Position in total Buffer
	unsigned int m_bufferPosY; //!< The calculated y Position in total Buffer
	float m_bufferWidth; //!< The total buffer width, stored as float for division sake
	float m_bufferHeight; //!< The total buffer height
	float m_progressedLastSegment; //!< The progress of last segment
	float m_progressedOverall; //!< The overall progress
};

}
