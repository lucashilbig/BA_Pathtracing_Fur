#pragma once

#include "KIRK/Utils/Log.h"
#include "externals/CTPL/ctpl_stl.h"
#include <thread>
#include <atomic>
#include <mutex>
#include <map>

namespace KIRK {
/**
*	A class realising the multithreading concept of a barrier.
*	The constructor is called with the amount of threads that need to be synchronised.
*	Each Thread can then call the member function count_down_and_wait() to synchronise.
*	Internally, a Busy-Wait is realised using yield.
*/
class SpinlockBarrier
{
public:
	SpinlockBarrier(const SpinlockBarrier&) = delete;
	SpinlockBarrier& operator=(const SpinlockBarrier&) = delete;

	/**
	*	The constructor. Initialises the Barrier with the amount of threads that are to be synchronised.
	*	@param count The number of threads.
	*/
	explicit SpinlockBarrier(unsigned int count) :
		m_count(check_counter(count)), m_generation(0),
		m_count_reset_value(count) {}

	/**
	*	Check whether count is zero. If so, throws an exception. Otherwise, count is returned.
	*	@param count The number to check.
	*	@return count
	*	@throws std::exception
	*/
	static unsigned int check_counter(unsigned int count)
	{
		if(count == 0)
			throw new std::invalid_argument("SpinlockBarrier constructor: count cannot be zero.");
		return count;
	}

	/**
	*	Yield thread until this function has been called as many times as specified in constructor.
	*	Resets the counter afterwards.
	*/
	void count_down_and_wait()
	{
		unsigned int gen = m_generation.load();
		if(--m_count == 0)
		{
			if(m_generation.compare_exchange_weak(gen, gen + 1))
			{
				m_count = m_count_reset_value;
			}
			return;
		}
		while((gen == m_generation) && (m_count != 0))
			std::this_thread::yield();
	}

private:
	std::atomic<unsigned int> m_count;
	std::atomic<unsigned int> m_generation;

	unsigned int m_count_reset_value;
};

/**
 * Here we can collect some functions dedicated to threading.
 */
class ThreadManager
{
public:
	ThreadManager() = delete;
	~ThreadManager() = delete;

	using thread_func_double = std::function<void(int, int)>;
	using thread_func_single = std::function<void(int)>;

	/**
	 * @return The maximum count of threads to be used on this Computer.
	 */
	static int maxThreads() { return std::thread::hardware_concurrency(); }

	/**
	 * Equivalent to two nested for loops like

	 \code{.cpp}
		for(int y=start_y; y<end_excl_y; y++)
			for(int x=start_x; x<end_excl_x; x++)
	 \endcode

	 * but where the workload is worked on by several threads line by line (y).
	 * @param start_x Start value of the inner for loop
	 * @param start_y Start value of the outer for loop
	 * @param end_excl_x Exclusive end value of the inner for loop
	 * @param end_excl_y Exclusive end value of the outer for loop
	 * @param func function which will be executed for each value pair int x and int y ( as in starting to raytrace for one pixel with coordinates (x,y) )
	 * @param thread_count Thread count to use. Defaults to maxThreads().
	   * @param use_progress_bar If true, a progress bar will be printed in the console. Defaults to false.
	 */
	static void for_loop_double(int start_x, int start_y, int end_excl_x, int end_excl_y, thread_func_double func, int thread_count = maxThreads(), bool use_progress_bar = false);

	/**
	 * \brief Calls func with values in [start, end_excl), spreading the interval into chunks that are assigned as jobs to threads.
	 * \param start The beginning of the interval to work on
	 * \param end_excl The end of the interval to work on, exclusive
	 * \param func The function to call
	 * \param thread_count How many threads to use
	 * \param use_progress_bar Whether to output a progress bar to std::cout, running in an additional thread
	 */
	static void for_loop_single(int start, int end_excl, thread_func_single func, int thread_count = maxThreads(), bool use_progress_bar = false);

private:
	static std::atomic_int m_counter; /*!< Helper counter for the loops */
	static ctpl::thread_pool m_threads; /*!< Thread pool */

	static void addProgressBar(int num_items, int divisor);
};

/**
 * A struct for wrapping thread-dependant variables.
 * Internally, there will be one value stored for each thread.
 * From the outside, this can be used just as if it actually was a variable of the type it is wrapping:
 * \code{.cpp}
 * ThreadDependant<bool> b1;
 * ThreadDependant<bool> b2;
 * b1 = false;
 * b2 = true;
 * bool b3 = b1 && b2; // b3 will be false
 *
 * ThreadDependant<int> a;
 * ThreadDependant<int> b;
 * a = 5;
 * b = 6;
 * int c = a + b; // c will be 11
 * std::cout << a << std::endl // will output 5
 * \endcode
 * The struct will automatically handle this and use/assign the value belonging to the thread it is called from.
 *
 * Important notice: if you wrap something and want to use the .-operator to access a member or member function, you need to change:
 * \code{.cpp}
 * variable.member
 * \endcode
 * to
 * \code{.cpp}
 * variable().member
 * \endcode
 * That is because the .-operator cannot be overloaded.
 *
 * If you start parallelising an existing Raytracer and need a (member) variable thread-dependant, all you need to do is change it to ThreadDependant<WhateverType> in the declaration, all the other code should still work.
 */
template<typename T>
struct ThreadDependant
{
public:
	ThreadDependant() : m_default_value() {}

	/**
	* Constructs a ThreadDependant where every entry in the internal array is initialised to the same value.
	* @param val The value to initialise all entrys.
	*
	*/
	ThreadDependant(T val) : m_default_value(val) {}
	~ThreadDependant() {}

	operator T&()
	{
		if(!value.count(std::this_thread::get_id())) { return value[std::this_thread::get_id()] = m_default_value; }
		return value[std::this_thread::get_id()];
	}

	T& operator()()
	{
		if(!value.count(std::this_thread::get_id())) { return value[std::this_thread::get_id()] = m_default_value; }
		return value[std::this_thread::get_id()];
	}

	ThreadDependant<T>& operator=(const T& arg)
	{
		value[std::this_thread::get_id()] = arg;
		return *this;
	}

	ThreadDependant<T>& operator=(ThreadDependant<T>& arg)
	{
		value[std::this_thread::get_id()] = arg();
		return *this;
	}

private:
	std::map<std::thread::id, T> value;
	T m_default_value;
};
}
