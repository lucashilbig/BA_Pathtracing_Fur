#include "KIRK/CPU/CPU_Raytracer/CPU_PathTracer.h"
#include "KIRK/Common/Shading/Shader.h"

namespace KIRK {
	namespace CPU {

		PathTracer::PathTracer(std::shared_ptr<KIRK::CPU::Scene> cpuscene, int depth) : PathTracer()
		{
			PathTracer::init(cpuscene);
		}

		// per sample:
		//		first generate primary rays.
		//		then trace them.
		//		then partly clear LightSampleBuffer except the colors.
		// Finally after all samples put the LightSamples onto the texture.
		void KIRK::CPU::PathTracer::render()
		{
			if ((c_sample == 0 && !m_bufferSegmentation.isActive()) || m_has_resized) { resetBuffers(); }


			if (c_sample < m_samples_per_pixel || m_bufferSegmentation.isActive())
			{
				LOG_INFO("Rendering Sample %/% -> Segment %/%",
					c_sample + 1, m_samples_per_pixel,
					m_bufferSegmentation.getSegmentId() + 1, m_bufferSegmentation.getSegmentCount());

				Clock<> clock;

				processSegment();
				drawTexture();
				m_bufferSegmentation.update();

				float cputime = clock.getElapsedTime<Milliseconds>();
				float progressLastSegment = m_bufferSegmentation.getProgressLastSegment(m_samples_per_pixel);
				float progressOverall = m_bufferSegmentation.getProgressOverall(m_samples_per_pixel);
				float estimatedTimeLeft = ((1.0f - progressOverall) / progressLastSegment) * cputime;

				LOG_INFO("Processed % %% in %ms \n\tProcess overall % %% \n\tEstimated time left: %ms\n",
					progressLastSegment * 100, cputime,
					progressOverall * 100, estimatedTimeLeft);

				if (!m_bufferSegmentation.isActive())
				{
					c_sample++;
					if (m_use_tonemapping)
						applyToneMapping();
				}

			}

		}

		void PathTracer::processSegment()
		{
			generatePrimaryRays();
			clearBufferWeights();
			traceRays();
		}

		void PathTracer::drawTexture()
		{
			KIRK::ThreadManager::for_loop_double(0, 0, m_bufferSegmentation.getSegmentWidth(), m_bufferSegmentation.getSegmentHeight(), [&](int x, int y)
			{
				int bufferIndex = m_bufferSegmentation.get1DBufferIndex(x, y);
				int segmentIndex = m_bufferSegmentation.get1DSegmentIndex(x, y);

				if (c_sample == 0)
				{
					m_textureBuffer[bufferIndex] = m_bounceBuffer[segmentIndex].color;
				}
				else
				{
					//m_textureBuffer[bufferIndex] += (m_bounceBuffer[segmentIndex].color - m_textureBuffer[bufferIndex]) / (c_sample + 1);
					//Alternative, yields same results.
					m_textureBuffer[bufferIndex] = m_textureBuffer[bufferIndex] * (c_sample / static_cast <float>(c_sample + 1));
					m_textureBuffer[bufferIndex] += (1 / static_cast <float>(c_sample + 1)) * m_bounceBuffer[segmentIndex].color;
				}




				{
					// Let's protect this with a mutex because it changes state in the texture.
					// Doesn't seem to change anything but this way we can be sure we are not doing things wrong
					std::unique_lock<std::mutex> l(m_texture_mutex);
					m_render_texture->setPixel(m_bufferSegmentation.getBufferPosX(x), m_bufferSegmentation.getBufferPosY(y),
						m_textureBuffer[bufferIndex]);
				}
			}, m_num_threads, m_use_progress_bar);
		}

		void PathTracer::applyToneMapping()
		{
			auto m_tonemapBuffer = std::vector <glm::vec3>(m_bufferSegmentation.getBufferSize());
			for (int x = 0; x < m_bufferSegmentation.getBufferWidth(); ++x) { for (int y = 0; y < m_bufferSegmentation.getBufferHeight(); ++y) { m_tonemapBuffer[y * m_bufferSegmentation.getBufferWidth() + x] = glm::vec3(m_textureBuffer[y * m_bufferSegmentation.getBufferWidth() + x]); } }
			m_tonemapper->map(m_tonemapBuffer, m_bufferSegmentation.getBufferWidth(), m_bufferSegmentation.getBufferHeight());
			for (int x = 0; x < m_bufferSegmentation.getBufferWidth(); ++x)
			{
				for (int y = 0; y < m_bufferSegmentation.getBufferHeight(); ++y)
				{
					m_render_texture->setPixel(x, y,
						glm::vec4(m_tonemapBuffer[y * m_bufferSegmentation.getBufferWidth() + x], 1));
				}
			}
		}

		void PathTracer::clearBufferWeights()
		{
			KIRK::ThreadManager::for_loop_double(0, 0, m_bufferSegmentation.getSegmentWidth(), m_bufferSegmentation.getSegmentHeight(), [&](int x, int y)
			{
				int segmentIndex = m_bufferSegmentation.get1DSegmentIndex(x, y);
				m_bounceBuffer[segmentIndex].radiance = glm::vec3(1.f);
				m_bounceBuffer[segmentIndex].color = glm::vec4(0.f);
				m_bounceBuffer[segmentIndex].mat_flags = 0;
			}, m_num_threads, m_use_progress_bar);
		}

		void PathTracer::generatePrimaryRays()
		{
			KIRK::ThreadManager::for_loop_double(0, 0, m_bufferSegmentation.getSegmentWidth(), m_bufferSegmentation.getSegmentHeight(), [&](int x, int y)
			{
				int segmentIndex = m_bufferSegmentation.get1DSegmentIndex(x, y);
				m_rayBuffer[segmentIndex] = m_cpuscene->getActiveCamera().getRayFromPixel(m_bufferSegmentation.getBufferPosX(x), m_bufferSegmentation.getBufferPosY(y), m_sampler.sample1D(), m_sampler.sample1D());
				//Reset m_lamba
				m_hitBuffer[segmentIndex].m_lambda = FLT_MAX;
			}, m_num_threads, m_use_progress_bar);
		}

		void PathTracer::traceRays()
		{
			for (int bounces = 0; bounces < m_depth; ++bounces)
			{
				KIRK::ThreadManager::for_loop_double(0, 0, m_bufferSegmentation.getSegmentWidth(), m_bufferSegmentation.getSegmentHeight(), [&](int x, int y)
				{
					//Check if bounce is still contributing to color at all
					int segmentIndex = m_bufferSegmentation.get1DSegmentIndex(x, y);
					if (m_bounceBuffer[segmentIndex].radiance != glm::vec3(0))
					{
						traceRay(segmentIndex, bounces == 0, m_bounceBuffer[segmentIndex].mat_flags);

						if (m_hitBuffer[segmentIndex].m_lambda == FLT_MAX)
						{
							//No Intersection found, shade Environment
							m_cpuscene->getEnvironment().m_eshader->shade(*this, m_hitBuffer[segmentIndex], m_bounceBuffer[segmentIndex], m_rayBuffer[segmentIndex]);
						}
						else
						{
							//We found an intersection, let's check if it's one with a light or regular mesh
							if (m_hitBuffer[segmentIndex].m_barycentric_coord.x == -1)
							{
								//We hit a light
								auto light = m_cpuscene->getLights()[m_hitBuffer[segmentIndex].m_barycentric_coord.y].get();
								light->m_lshader->shade(*this, m_hitBuffer[segmentIndex], m_bounceBuffer[segmentIndex], m_rayBuffer[segmentIndex]);
							}
							else
							{
								//We actually found an regular intersection, shade it.
								m_hitBuffer[segmentIndex].m_object->getMaterial()->m_shader->shade(*this, m_hitBuffer[segmentIndex], m_bounceBuffer[segmentIndex], m_rayBuffer[segmentIndex]);
								if ((m_bounceBuffer[segmentIndex].mat_flags & BSDFHelper::MATFLAG_EMISSIVE_BOUNCE) == BSDFHelper::MATFLAG_EMISSIVE_BOUNCE)
								{
									m_hitBuffer[segmentIndex].m_lambda = -1;
								}
							}
						}
					}
				}, m_num_threads, m_use_progress_bar);
			}
		}

		void PathTracer::traceRay(int segmentIndex, bool is_primary, int mat_flags)
		{
			if (m_rayBuffer[segmentIndex].m_direction == glm::vec3(0))
				m_hitBuffer[segmentIndex].m_lambda = -1;
			if (m_hitBuffer[segmentIndex].m_lambda == -1) { return; }
			m_hitBuffer[segmentIndex].m_ray = m_rayBuffer[segmentIndex];
			m_hitBuffer[segmentIndex].m_lambda = FLT_MAX;
			m_hitBuffer[segmentIndex].m_object = nullptr;
			bool is_intersection = m_cpuscene->getDataStructure().closestIntersection(&m_hitBuffer[segmentIndex]);

			if (is_intersection)
			{
				m_hitBuffer[segmentIndex].m_object->calcNormal(&m_hitBuffer[segmentIndex]);
				m_hitBuffer[segmentIndex].m_object->calcTcoord(&m_hitBuffer[segmentIndex]);
			}
			bool is_light_intersection = false;
			float t_lights = FLT_MAX;
			int t_index = -1;
			for (int lightIndex = 0; lightIndex < m_cpuscene->getLights().size(); lightIndex++)
			{
				auto light = m_cpuscene->getLights()[lightIndex].get();
				float t = FLT_MAX;
				bool intsect = light->isIntersection(m_hitBuffer[segmentIndex].m_ray, t);
				is_light_intersection = is_light_intersection || intsect;
				if (intsect)
				{
					t_lights = glm::min(t_lights, t);
					t_index = t == t_lights ? lightIndex : t_index;
				}
			}
			if (t_lights < m_hitBuffer[segmentIndex].m_lambda)
			{
				// Set lamba to distance of nearest light intersection
				m_hitBuffer[segmentIndex].m_lambda = t_lights;
				//Use barycentric coord x as indicator that we hit a light
				m_hitBuffer[segmentIndex].m_barycentric_coord.x = -1;
				//Use barycentric coord y to save lightIndex
				m_hitBuffer[segmentIndex].m_barycentric_coord.y = t_index;
			}
		}

		void PathTracer::resetBuffers()
		{
			initBufferSegmentation(sizeof(glm::vec4), sizeof(Bounce) + sizeof(Ray) + sizeof(Intersection));
			m_rayBuffer = std::vector <Ray>(m_bufferSegmentation.getSegmentSize());
			m_textureBuffer = std::vector <glm::vec4>(m_bufferSegmentation.getBufferSize());
			m_bounceBuffer = std::vector <Bounce>(m_bufferSegmentation.getSegmentSize());
			m_hitBuffer = std::vector <Intersection>(m_bufferSegmentation.getSegmentSize(), Intersection(Ray()));

			m_has_resized = false;

			LOG_DEBUG("Max Buffer Size: %, Actual Buffer Size: %", m_maxBufferSize,
				sizeof(m_rayBuffer[0]) * m_rayBuffer.capacity()
				+ sizeof(m_bounceBuffer[0]) * m_bounceBuffer.capacity()
				+ sizeof(m_textureBuffer[0]) * m_textureBuffer.capacity()
				+ sizeof(m_hitBuffer[0]) * m_hitBuffer.capacity()
			);
		}

		void PathTracer::initBufferSegmentation(size_t staticBufferObjectSize, size_t dynamicBufferObjectSize)
		{
			try
			{
				m_bufferSegmentation.init(m_render_texture->getSize().x, m_render_texture->getSize().y,
					m_maxBufferSize, staticBufferObjectSize, dynamicBufferObjectSize);
			}
			catch (std::runtime_error e)
			{
				LOG_ERROR("%", e.what());
				std::exit(EXIT_FAILURE);
			}
		}

		void PathTracer::windowResizeCallback() { m_has_resized = true; }

		void PathTracer::reset()
		{
			c_sample = 0;
			m_bufferSegmentation.setInactive();
		}

		void PathTracer::setSampleCount(const unsigned int samples)
		{
			m_samples_per_pixel = samples;
			m_inverse_samples_per_pixel = 1 / (float)samples;
		}

		int PathTracer::getCurrentSampleCount()
		{
			return c_sample;
		}

		void PathTracer::onGui()
		{
			int samples = m_samples_per_pixel;
			int maxBufferSize = static_cast <int>(m_maxBufferSize * 1e-6);

			ImGui::PushID("pathtracer_cpu");
			if (ImGui::CollapsingHeader("Pathtracer")) {
				if (ImGui::TreeNode("Tools"))
				{
					ImGui::SliderInt("Thread count", &m_num_threads, 1, ThreadManager::maxThreads());
					ImGui::Checkbox("Progressbar", &m_use_progress_bar);
					ImGui::DragInt("Max Buffer Size (MB)", &maxBufferSize, 10, 10, 48000);
					ImGui::TreePop();
				}

				if (ImGui::TreeNode("Sampling"))
				{
					ImGui::DragInt("Samples", &samples, 1, 1, 10000);
					ImGui::SliderInt("Max Bounces", &m_depth, 1, 20);
					ImGui::TreePop();
				}

				if (m_tonemapper && ImGui::TreeNode("ToneMapping"))
				{
					ImGui::Checkbox("Tonemap every frame", &m_use_tonemapping);
					if (ImGui::Button("Tonemap current frame", ImVec2(ImGui::GetContentRegionAvailWidth(), 32)))
					{
						applyToneMapping();
					}
					m_tonemapper->draw();
					ImGui::TreePop();
				}

			}
			ImGui::PopID();
			m_samples_per_pixel = samples;
			m_maxBufferSize = static_cast <size_t>(maxBufferSize) * 1e+6;
		}

	}
}
