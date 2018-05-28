#ifndef __CLOCK_H__
#define __CLOCK_H__

#include <chrono>
#include <stack>
#include <vector>
#include <GL/glew.h>

/**
 * \brief Starts the GPU Clock.
 */
#define START_GPU_CLOCK KIRK::GPU_Clock::instance().start()
/**
 * \brief Ends the GPU Clock, returning the measured time in nanoseconds.
 */
#define END_GPU_CLOCK KIRK::GPU_Clock::instance().end()
/**
 * \brief Starts the GPU Clock only in Debug-Mode.
 */
#define START_GPU_CLOCK_DEBUG KIRK::GPU_Clock::instance().start()
/**
 * \brief Ends the GPU Clock, returning the measured time in nanoseconds. Only in Debug-Mode.
 */
#define END_GPU_CLOCK_DEBUG KIRK::GPU_Clock::instance().end()
#ifdef NDEBUG
#undef START_GPU_CLOCK_DEBUG
#undef END_GPU_CLOCK_DEBUG
#define START_GPU_CLOCK_DEBUG ;
#define END_GPU_CLOCK_DEBUG 0
#endif

namespace KIRK {
	typedef std::chrono::nanoseconds Nanoseconds;
	typedef std::chrono::microseconds Microseconds;
	typedef std::chrono::milliseconds Milliseconds;
	typedef std::chrono::duration<float> Seconds;
	typedef std::chrono::duration<double, std::ratio<60>> Minutes;
	typedef std::chrono::duration<long double, std::ratio<3600>> Hours;

	/**
	 * @brief A simple timer that can be used for measuring runtime of code fragments.
	 *
	 * @example
	 * CVK::Clock < std::chrono::steady_clock > clock;
	 * ... // do something
	 * LOG_INFO << clock.getElapsedTime < CVK::Milliseconds > () << " ms";
	 */
	template<typename C = std::chrono::steady_clock>
	class Clock
	{
	public:
		/**
		 * @brief Clock class constructor; essentially the same thing as CVK::Clock::restart
		 */
		Clock() : m_start(C::now()) {};

		/**
		 * @brief Get the time elapsed since the clock was last restarted (or constructed)
		 * @return the elapsed time as internal representation from std::chrono
		 */
		template<typename T = Milliseconds>
		typename T::rep getElapsedTime() const
		{
			return std::chrono::duration_cast<T>(C::now() - m_start).count();
		};

		/**
		 * @brief Restart the clock so that next time we call getElapsedTime () less (or equally much) time has elapsed
		 */
		void restart() { m_start = C::now(); };
	private:
		std::chrono::time_point<C> m_start;
	};

	/**
	 * This Clock measures the time of any OpenGL-Calls executed in between of start() and end().
	 * end() also returns the time as an GLuint64 measured in nanoseconds.
	 * This Clock can handle nested start/end blocks in each other of arbitrary depth and nestedness.
	 * The Clock effectively ignores the time it takes to measure time itself when doing nested measuring.
	 * However, time measuring on the GPU consumes more time than on the CPU.
	 * Because of that, you should limit your measurements.
	 * Usage through the macros START_GPU_CLOCK and END_GPU_CLOCK is recommended.
	 * For detailed measurements, use START_GPU_CLOCK_DEBUG and END_GPU_CLOCK_DEBUG to avoid performance issues in Release.
	 * \brief A Clock for measuring time on the GPU.
	 */
	class GPU_Clock {
	public:
		/**
		 * \brief Retrieves/Creates the Singleton instance.
		 * \return the instance.
		 */
		static GPU_Clock& instance();
		~GPU_Clock();

		GPU_Clock(GPU_Clock&&) = delete;
		GPU_Clock& operator=(GPU_Clock&& other) = delete;
		// automatically no copy operations now

		/**
		 * \brief Start the Clock.
		 */
		void start();

		/**
		 * \brief End the Clock.
		 * \return The time measured in Nanoseconds.
		 */
		GLuint64 end();

	private:
		/**
		 * \brief Wrapper-Class for GPU-Queries with GL_TIME_ELAPSED.
		 * Cannot be nested, thus this class is private because GPU_Clock works around that problem.
		 */
		class TimerQuery {
		public:

			TimerQuery() { glGenQueries(1, &m_id); };
			explicit TimerQuery(GLuint id) : m_id(id) { };
			~TimerQuery() { glDeleteQueries(1, &m_id); };

			TimerQuery(TimerQuery&&) = delete;
			TimerQuery& operator=(TimerQuery&& other) = delete;
			// automatically no copy operations now

			/**
			 * \brief Begins the Query GL_TIME_ELAPSED.
			 */
			void issue() const { glBeginQuery(GL_TIME_ELAPSED, m_id); };

			/**
			 * \brief Ends the Query GL_TIME_ELAPSED, waits for the result and returns it.
			 * \return The elapsed time since issue(), measured in Nanoseconds.
			 */
			GLuint64 get() const {
				glEndQuery(GL_TIME_ELAPSED);
				GLuint64 result;
				glGetQueryObjectui64v(m_id, GL_QUERY_RESULT, &result);
				return result;
			};

		private:
			GLuint m_id;
		};

		GPU_Clock();

		TimerQuery m_query;
		int m_current_index = 0;
		std::vector<GLuint64> m_times;
		std::stack<int> m_start_indices;
	};
}

#endif
