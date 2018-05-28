#include "Clock.h"
#include "Log.h"
#include <numeric>

KIRK::GPU_Clock::GPU_Clock()
{
	LOG_DEBUG("Constructing GPU_Clock instance...");
}

KIRK::GPU_Clock& KIRK::GPU_Clock::instance()
{
	static KIRK::GPU_Clock _instance;
	return _instance;
}

void KIRK::GPU_Clock::start()
{
	if(!m_start_indices.empty())
	{
		m_times.push_back(m_query.get());
		m_current_index++;
	}
	m_start_indices.push(m_current_index);
	m_query.issue();
}

GLuint64 KIRK::GPU_Clock::end()
{
	m_current_index++;
	m_times.push_back(m_query.get());
	int index = m_start_indices.top();
	m_start_indices.pop();
	GLuint64 elapsed_time = 0;
	elapsed_time = std::accumulate(m_times.begin() + index, m_times.begin() + m_current_index, elapsed_time);
	if(m_start_indices.empty())
	{
		m_current_index = 0;
		m_times.clear();
	}
	else { m_query.issue(); }
	return elapsed_time;
}

KIRK::GPU_Clock::~GPU_Clock() {}
