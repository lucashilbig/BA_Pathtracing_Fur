#include "Threading.h"

std::atomic_int KIRK::ThreadManager::m_counter;
ctpl::thread_pool KIRK::ThreadManager::m_threads(maxThreads());

void KIRK::ThreadManager::for_loop_single(int start, int end_excl, thread_func_single func, int thread_count, bool use_progress_bar)
{
	// Calculate number of work items
	const int chunk_size =  ((end_excl - start) % thread_count > 0 ? 1 : 0) + (end_excl - start) / thread_count;
	const int final_chunk_size = (end_excl - start) % chunk_size;
	const int num_full_chunks = (end_excl - start) / chunk_size;
	const int num_chunks = (final_chunk_size > 0 ? 1 : 0) + num_full_chunks;

	// Resize thread pool if needed
	if(m_threads.size() != thread_count + use_progress_bar)
	{
		m_threads.resize(thread_count + use_progress_bar);
	}

	m_counter = num_chunks;

	// Add progress bar if needed
	if(use_progress_bar)
	{
		addProgressBar(num_chunks, 30);
	}

	// Add tasks to the thread pool
	for(int i = 0; i < num_full_chunks; ++i)
	{
		m_threads.push([i, chunk_size, func, start](int id)
			{
				for(int x = 0; x < chunk_size; ++x)
				{
					func(start + x + i * chunk_size);
				}
				--m_counter;
			});
	}

	// Add the last, smaller chunk, in case it exists
	if(num_full_chunks < num_chunks)
	{
		m_threads.push([func, chunk_size, num_full_chunks, final_chunk_size, start](int id)
			{
				for(int x = 0; x < final_chunk_size; ++x)
				{
					func(start + x + num_full_chunks * chunk_size);
				}
				--m_counter;
			});
	}

	// Wait until all tasks are finished
	while(m_counter > 0)
	{
		std::this_thread::yield();
	}
}

void KIRK::ThreadManager::addProgressBar(int num_items, int divisor)
{
	int progress_rate = num_items / divisor;
	m_threads.push([progress_rate, divisor](int id)
		{
			int old_progress = -1;
			while(true)
			{
				int current_counter = m_counter;
				if(current_counter < 0)
				{
					break;
				}
				auto div_r = std::div(current_counter, progress_rate);
				if(old_progress < divisor - div_r.quot)
				{
					old_progress = divisor - div_r.quot;
					std::string output = "| " + std::string(old_progress, '=') + std::string(divisor - old_progress, ' ') + " |\r";
					std::cout << output << std::flush;
				}
				//std::this_thread::sleep_for(std::chrono::milliseconds(100));
				/* How yield works depends on the scheduler of the OS */
				std::this_thread::yield();
			}
		});
}

void KIRK::ThreadManager::for_loop_double(int start_x, int start_y, int end_excl_x, int end_excl_y, thread_func_double func, int thread_count, bool use_progress_bar)
{
	// Resize thread pool if needed
	if(m_threads.size() != thread_count + use_progress_bar)
	{
		m_threads.resize(thread_count + use_progress_bar);
	}

	m_counter = end_excl_y - start_y - 1;

	// Add progress bar as task to thread pool if needed
	if(use_progress_bar)
	{
		addProgressBar(end_excl_y - start_y, 30);
	}

	// Add tasks to the thread pool
	for(int i = start_y; i < end_excl_y; ++i)
	{
		m_threads.push([i, func, start_x, end_excl_x](int id)
			{
				for(int x = start_x; x < end_excl_x; ++x)
				{
					func(x, i);
				}
				--m_counter;
			});
	}

	// Wait until all tasks are finished
	while(m_counter >= 0)
	{
		std::this_thread::yield();
	}
}
